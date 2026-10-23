//! Scoring matrices, gap costs and Gumbel statistics.
//!
//! The standard tables and their published ALP parameters are generated from
//! DIAMOND's own `src/stats/matrices/*.h` files by `scripts/gen_matrices.py`,
//! so E-values here are calibrated exactly like the C++ implementation.

use crate::alp::AlpParams;
use crate::bail;
use crate::error::{Error, Result};
use crate::matrix_data::{StdMatrix, MATRICES, SYMBOLS};
use std::path::Path;

pub const LN_2: f64 = std::f64::consts::LN_2;

/// BLOSUM62 background frequencies (NCBI), in [`SYMBOLS`] order.
const BLOSUM62_BACKGROUND: [f64; 20] = [
    0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025,
    0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073,
];

#[derive(Clone, Debug)]
pub struct ScoreMatrix {
    pub name: String,
    /// Row-major `25 x 25` scores in NCBI order.
    scores: Vec<i8>,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// Gapped Karlin-Altschul parameters.
    pub lambda: f64,
    pub k: f64,
    ln_k: f64,
    /// ALP parameters used for DIAMOND-compatible E-values; `None` for custom
    /// matrices, which fall back to the textbook formula.
    alp: Option<AlpParams>,
}

impl ScoreMatrix {
    /// Load one of the built-in matrices (`BLOSUM62`, `BLOSUM80`, ...).
    ///
    /// `gap_open`/`gap_extend` of `-1` select the matrix defaults, matching
    /// DIAMOND's command line.
    pub fn standard(name: &str, gap_open: i32, gap_extend: i32) -> Result<ScoreMatrix> {
        let upper = name.to_ascii_uppercase();
        let def: &StdMatrix = MATRICES
            .iter()
            .find(|m| m.name == upper)
            .ok_or_else(|| Error::msg(format!("Unknown scoring matrix: {name}")))?;
        let gap_open = if gap_open == -1 {
            def.default_gap.0
        } else {
            gap_open
        };
        let gap_extend = if gap_extend == -1 {
            def.default_gap.1
        } else {
            gap_extend
        };
        if gap_open <= 0 || gap_extend <= 0 {
            bail!("gap penalties must be positive (got {gap_open}/{gap_extend})");
        }
        let params = def
            .params
            .iter()
            .find(|p| p.0 == gap_open && p.1 == gap_extend)
            .ok_or_else(|| {
                Error::msg(format!(
                    "Gap penalty settings are outside the supported range for {} (--gapopen {} --gapextend {}).",
                    def.name, gap_open, gap_extend
                ))
            })?;
        let mut scores = Vec::with_capacity(25 * 25);
        for row in def.scores.iter() {
            scores.extend_from_slice(row);
        }
        let lambda = params.2;
        let k = params.3;
        let alp = Some(AlpParams::new(
            lambda,
            k,
            params.4,
            params.5,
            params.6,
            gap_open,
            gap_extend,
            def.ungapped.0,
            def.ungapped.1,
        ));
        Ok(ScoreMatrix {
            name: def.name.to_string(),
            scores,
            gap_open,
            gap_extend,
            lambda,
            k,
            ln_k: k.ln(),
            alp,
        })
    }

    /// Read an NCBI-format matrix file (`--custom-matrix`).
    ///
    /// Statistics are estimated: `lambda` solves the ungapped fixed-point
    /// equation with BLOSUM62 background frequencies, `K` uses a conservative
    /// default. DIAMOND runs the full ALP fit here; exact statistics for a
    /// custom matrix should therefore be taken as approximate.
    pub fn from_file(path: &Path, gap_open: i32, gap_extend: i32) -> Result<ScoreMatrix> {
        if gap_open <= 0 || gap_extend <= 0 {
            bail!("gap penalties must be positive (got {gap_open}/{gap_extend})");
        }
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::msg(format!("cannot read {}: {e}", path.display())))?;
        let mut columns: Option<Vec<u8>> = None;
        let mut scores = vec![-gap_extend as i8; 25 * 25];
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut tokens = line.split_whitespace();
            let Some(first) = tokens.next() else { continue };
            if columns.is_none() {
                let mut cols = Vec::new();
                for token in std::iter::once(first).chain(tokens) {
                    cols.push(symbol_index(token)?);
                }
                columns = Some(cols);
                continue;
            }
            let row = symbol_index(first)?;
            let cols = columns.as_ref().unwrap();
            for (j, token) in tokens.enumerate() {
                if j >= cols.len() {
                    break;
                }
                let value: i32 = token
                    .parse()
                    .map_err(|_| Error::msg(format!("bad matrix value {token:?}")))?;
                scores[row as usize * 25 + cols[j] as usize] = value as i8;
            }
        }
        if columns.is_none() {
            bail!("{}: empty matrix", path.display());
        }
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("custom")
            .to_string();
        let lambda = solve_lambda(&scores, &BLOSUM62_BACKGROUND).unwrap_or(0.267);
        let k = 0.1;
        Ok(ScoreMatrix {
            name,
            scores,
            gap_open,
            gap_extend,
            lambda,
            k,
            ln_k: k.ln(),
            alp: None,
        })
    }

    #[inline]
    pub fn score(&self, a: u8, b: u8) -> i32 {
        let a = (a as usize).min(24);
        let b = (b as usize).min(24);
        self.scores[a * 25 + b] as i32
    }

    #[inline]
    pub fn score_i8(&self, a: u8, b: u8) -> i8 {
        self.scores[(a as usize).min(24) * 25 + (b as usize).min(24)]
    }

    /// Cost of a gap of `len` residues (BLAST convention: open + len * extend).
    #[inline]
    pub fn gap_cost(&self, len: i32) -> i32 {
        self.gap_open + len * self.gap_extend
    }

    /// Raw score of an alignment that scores `bits` bit score.
    pub fn raw_score(&self, bits: f64) -> i32 {
        ((bits * LN_2 + self.ln_k) / self.lambda).ceil() as i32
    }

    /// BLAST bit score: `(lambda * S - ln K) / ln 2`.
    pub fn bit_score(&self, raw: i32) -> f64 {
        (self.lambda * raw as f64 - self.ln_k) / LN_2
    }

    /// DIAMOND's E-value, including the ALP finite-size correction.
    ///
    /// `query_len`/`subject_len` are the full sequence lengths (not the HSP
    /// spans); `db_letters` is the number of residues in the database.
    pub fn evalue(&self, raw: i32, query_len: usize, subject_len: usize, db_letters: u64) -> f64 {
        match &self.alp {
            Some(alp) => alp.evalue(raw, query_len, subject_len, db_letters),
            None => self.evalue_blast(raw, query_len, db_letters),
        }
    }

    /// Textbook Karlin-Altschul E-value, used for custom matrices and as a
    /// fallback when the ALP area degenerates.
    pub fn evalue_blast(&self, raw: i32, query_len: usize, db_letters: u64) -> f64 {
        let exponent = -self.lambda * raw as f64 + self.ln_k;
        if exponent > 700.0 {
            f64::INFINITY
        } else {
            query_len as f64 * db_letters as f64 * exponent.exp()
        }
    }

    /// Minimum raw score that keeps the ungapped E-value below `threshold`.
    pub fn raw_score_for_evalue(&self, threshold: f64, query_len: usize, db_letters: u64) -> i32 {
        if threshold <= 0.0 {
            return 0;
        }
        let space = query_len.max(1) as f64 * db_letters.max(1) as f64;
        (((space * self.k / threshold).ln()) / self.lambda).ceil() as i32
    }

    /// Expected score under the background distribution; must be negative for
    /// local alignment statistics to be valid.
    pub fn expected_score(&self) -> f64 {
        let mut e = 0.0;
        for (i, pi) in BLOSUM62_BACKGROUND.iter().enumerate() {
            for (j, pj) in BLOSUM62_BACKGROUND.iter().enumerate() {
                e += *pi * *pj * self.score(i as u8, j as u8).max(-20) as f64;
            }
        }
        e
    }
}

fn symbol_index(token: &str) -> Result<u8> {
    let c = token.as_bytes().first().copied().unwrap_or(b'?');
    let idx = SYMBOLS.iter().position(|&s| s == c.to_ascii_uppercase());
    match idx {
        Some(i) => Ok(i as u8),
        None => bail!("invalid symbol {token:?} in custom matrix"),
    }
}

/// Solve `sum_ij p_i p_j exp(lambda * s_ij) = 1` for `lambda > 0`.
fn solve_lambda(scores: &[i8], background: &[f64; 20]) -> Option<f64> {
    let f = |lambda: f64| {
        let mut total = 0.0;
        for i in 0..20 {
            for j in 0..20 {
                total += background[i] * background[j] * (lambda * scores[i * 25 + j] as f64).exp();
            }
        }
        total - 1.0
    };
    let mut hi = 1.0;
    for _ in 0..32 {
        if f(hi) > 0.0 {
            break;
        }
        hi *= 2.0;
    }
    let mut lo = 1e-6;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blosum62_defaults() {
        let m = ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap();
        assert_eq!((m.gap_open, m.gap_extend), (11, 1));
        assert_eq!(m.score(0, 0), 4); // A/A
        assert_eq!(m.score(17, 17), 11); // W/W
        assert_eq!(m.score(0, 17), -3); // A/W
        assert!((m.lambda - 0.267).abs() < 1e-9);
        assert!((m.k - 0.041).abs() < 1e-9);
    }

    #[test]
    fn bit_scores_are_blast_compatible() {
        let m = ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap();
        // (0.267 * 100 - ln 0.041) / ln 2 = 43.1 bits.
        let b = m.bit_score(100);
        assert!((b - 43.1).abs() < 0.1, "{b}");
        assert_eq!(m.raw_score(b), 100);
    }

    #[test]
    fn evalues_match_diamond_alp_statistics() {
        let m = ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap();
        let got = m.evalue(198, 497, 601, 2_263_996);
        assert!((got / 7.60e-17 - 1.0).abs() < 0.01, "{got}");
    }

    #[test]
    fn blast_evalues_match_the_formula() {
        let m = ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap();
        let expected = 0.041 * 100.0 * 1e6 * (-0.267f64 * 100.0).exp();
        let got = m.evalue_blast(100, 100, 1_000_000);
        assert!((got / expected - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unsupported_gap_costs_are_rejected() {
        assert!(ScoreMatrix::standard("BLOSUM62", 42, 7).is_err());
    }

    #[test]
    fn unknown_matrix_is_rejected() {
        assert!(ScoreMatrix::standard("NOPE", -1, -1).is_err());
    }
}
