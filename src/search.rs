//! The search pipeline: seeds → ungapped X-drop → gapped Smith-Waterman.

use crate::bail;
use crate::db::Database;
use crate::error::{Error, Result};
use crate::gapped;
use crate::hsp::Hsp;
use crate::matrix::ScoreMatrix;
use crate::reduction::Reduction;
use crate::seed::{self, SeedHit};
use crate::shape;
use crate::ungapped;
use rayon::prelude::*;
use std::collections::HashMap;

/// Sensitivity presets. The seed shapes, reduction and seed-complexity cut
/// mirror `Search::sensitivity_traits` in the C++ source.
#[derive(Clone, Debug)]
pub struct SearchConfig {
    pub matrix: ScoreMatrix,
    pub evalue: f64,
    pub max_target_seqs: usize,
    pub max_hsps: usize,
    pub sensitivity: String,
    pub threads: usize,
    /// `--index-chunks`; `None` keeps the sensitivity preset's value.
    pub index_chunks: Option<usize>,
    /// `--ungapped-evalue`; only used by DIAMOND's diagonal re-seeking mode,
    /// so the value has no effect here.
    pub ungapped_evalue: f64,
    /// `--seed-cut` override; `None` keeps the sensitivity preset's value.
    pub seed_cut: Option<f64>,
}

impl SearchConfig {
    pub fn new(matrix: ScoreMatrix) -> SearchConfig {
        SearchConfig {
            matrix,
            evalue: 0.001,
            max_target_seqs: 25,
            max_hsps: 1,
            sensitivity: "default".to_string(),
            threads: 1,
            index_chunks: None,
            ungapped_evalue: 10000.0,
            seed_cut: None,
        }
    }

    /// `--xdrop 12.3` bits, converted to raw score units.
    pub fn ungapped_xdrop(&self) -> i32 {
        self.matrix.raw_score(12.3)
    }

    fn seed_cut_for(&self) -> f64 {
        // `Search::sensitivity_traits`: --fast/--faster use 0.9, --default 0.8
        // and every more sensitive preset 1.0 (only homopolymer seeds masked).
        if let Some(cut) = self.seed_cut {
            return cut;
        }
        match self.sensitivity.as_str() {
            "faster" | "fast" => 0.9,
            "mid-sensitive" | "sensitive" | "more-sensitive" | "very-sensitive"
            | "ultra-sensitive" => 1.0,
            _ => 0.8,
        }
    }

    fn index_chunks_for(&self) -> usize {
        // `Search::sensitivity_traits`: the most sensitive presets index the
        // whole database in one slice so the join sees every seed at once.
        if let Some(n) = self.index_chunks {
            return n.max(1);
        }
        match self.sensitivity.as_str() {
            "very-sensitive" | "ultra-sensitive" => 1,
            _ => 4,
        }
    }
}

/// One translated frame of one query sequence.
#[derive(Clone, Debug)]
pub struct Target {
    pub query: usize,
    pub frame: u8,
    pub letters: Vec<u8>,
}

/// Align every target frame against the database.
///
/// The returned vector is parallel to `targets`; each entry holds the HSPs of
/// that frame sorted by score. `--max-target-seqs`/`--max-hsps` are applied by
/// the caller after frames have been merged.
pub fn run(db: &Database, targets: &[Target], cfg: &SearchConfig) -> Result<Vec<Vec<Hsp>>> {
    if targets.is_empty() || db.is_empty() {
        return Ok(vec![Vec::new(); targets.len()]);
    }
    let shapes = shape::shapes_for(&cfg.sensitivity, 0);
    let reduction = Reduction::for_sensitivity(&cfg.sensitivity);
    let seed_cut = cfg.seed_cut_for();
    // `setup.cpp`: `cfg.seed_complexity_cut = seed_cut * ln(2) * shapes[0].weight`.
    // `seed_cut` is expressed in bits per seed position while [`seed::seed_is_complex`]
    // compares a log multinomial coefficient, so the two must be converted.
    let complexity_cut = seed_cut * std::f64::consts::LN_2 * shapes[0].weight as f64;
    let index_chunks = cfg.index_chunks_for();

    // Query seed arrays are independent of the database chunk.
    let query_seqs: Vec<&[u8]> = targets.iter().map(|t| t.letters.as_slice()).collect();
    let mut query_seeds = Vec::with_capacity(shapes.len());
    for shape in &shapes {
        query_seeds.push(seed::build(
            &query_seqs,
            shape,
            &reduction,
            0,
            complexity_cut,
        ));
    }

    let refs: Vec<&[u8]> = (0..db.len()).map(|i| db.seq(i)).collect();

    // Double-indexed join, chunked so only one slice of the database is
    // indexed at a time (`--index-chunks`).
    let mut hits: Vec<HashMap<u32, Vec<SeedHit>>> =
        (0..targets.len()).map(|_| HashMap::new()).collect();
    for (c0, c1) in seed::chunks(db.len(), index_chunks) {
        for (si, shape) in shapes.iter().enumerate() {
            let ref_seeds =
                seed::build(&refs[c0..c1], shape, &reduction, c0 as u32, complexity_cut);
            seed::join_into(&query_seeds[si], &ref_seeds, &mut hits, 4096);
        }
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads.max(1))
        .build()
        .map_err(|e| Error::msg(format!("cannot start thread pool: {e}")))?;

    Ok(pool.install(|| {
        (0..targets.len())
            .into_par_iter()
            .map(|ti| align_target(db, &targets[ti], &hits[ti], cfg))
            .collect()
    }))
}

fn align_target(
    db: &Database,
    target: &Target,
    hits: &HashMap<u32, Vec<SeedHit>>,
    cfg: &SearchConfig,
) -> Vec<Hsp> {
    let query = &target.letters;
    if query.is_empty() || hits.is_empty() {
        return Vec::new();
    }
    let db_letters = db.total_letters();
    let matrix = &cfg.matrix;
    // DIAMOND's chaining backtrace only emits HSPs scoring at least 19 raw
    // score units (`Chaining::run(..., 19, band)`); applying the same floor to
    // the segments keeps weak diagonal noise out of the gapped stage, which
    // matters because every segment costs a banded DP here.
    let cutoff = 19;
    let xdrop = cfg.ungapped_xdrop();

    let mut subjects: Vec<u32> = hits.keys().copied().collect();
    subjects.sort_unstable();

    // Per subject: ungapped segments, then one gapped extension each.
    let mut found: Vec<Vec<Hsp>> = Vec::new();
    for sid in subjects {
        let subject = db.seq(sid as usize);
        if subject.is_empty() {
            continue;
        }
        let mut seeds: Vec<&SeedHit> = hits[&sid].iter().collect();
        seeds.sort_unstable_by_key(|h| (h.qpos as isize - h.spos as isize, h.spos));

        // Diagonal segments; hits that fall inside an already extended
        // segment on the same diagonal are sub-optimal (DIAMOND's adjacency
        // filter).
        let mut segments: Vec<ungapped::Segment> = Vec::new();
        for seed_hit in seeds {
            let diag = seed_hit.qpos as isize - seed_hit.spos as isize;
            if let Some(last) = segments.last() {
                if last.diag() == diag && last.s_end() >= seed_hit.spos as usize {
                    continue;
                }
            }
            let seg = ungapped::xdrop(
                query,
                subject,
                seed_hit.qpos as usize,
                seed_hit.spos as usize,
                matrix,
                xdrop,
            );
            if seg.score >= cutoff {
                segments.push(seg);
            }
        }
        if segments.is_empty() {
            continue;
        }
        segments.sort_unstable_by_key(|s| std::cmp::Reverse(s.score));

        // DIAMOND sizes the band from the query length alone, then merges the
        // bands of a target while they overlap (`min-band-overlap` defaults to
        // 0, so touching intervals merge). One DP per band replaces one DP per
        // segment, which is what lets a single HSP span several diagonals.
        let band = gapped::band_for(query.len());
        let mut order: Vec<usize> = (0..segments.len()).collect();
        order.sort_unstable_by_key(|&i| segments[i].diag());
        let mut bands: Vec<(isize, isize, usize)> = Vec::new();
        for &idx in &order {
            let diag = segments[idx].diag();
            let (b0, b1) = (diag - band as isize, diag + band as isize);
            match bands.last_mut() {
                Some(last) if b0 <= last.1 => {
                    last.1 = last.1.max(b1);
                    if segments[idx].score > segments[last.2].score {
                        last.2 = idx;
                    }
                }
                _ => bands.push((b0, b1, idx)),
            }
        }
        // Give the best-scoring band the first shot at --max-hsps.
        bands.sort_unstable_by(|a, b| segments[b.2].score.cmp(&segments[a.2].score));

        let mut hsps: Vec<Hsp> = Vec::new();
        for (d0, d1, best) in bands {
            if hsps.len() >= cfg.max_hsps {
                break;
            }
            let seg = &segments[best];
            let anchor = (seg.q_start, seg.s_start);
            let Some(hsp) = gapped_hsp(
                query,
                subject,
                anchor,
                (d0, d1),
                band,
                cfg,
                db_letters,
                target,
                sid,
            ) else {
                continue;
            };
            if hsp.evalue <= cfg.evalue {
                hsps.push(hsp);
            }
        }
        if !hsps.is_empty() {
            hsps.sort_by(cmp_hsp);
            found.push(hsps);
        }
    }

    // Order targets by their best HSP and apply --max-target-seqs.
    found.sort_by(|a, b| cmp_hsp(&a[0], &b[0]));
    found.truncate(cfg.max_target_seqs);
    let mut out = Vec::new();
    for mut hsps in found {
        hsps.truncate(cfg.max_hsps);
        out.extend(hsps);
    }
    out
}

fn cmp_hsp(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .then_with(|| a.subject.cmp(&b.subject))
}

// Nine parameters is a lot, but they are the natural inputs of one band: both
// sequences, the window it covers and the bookkeeping needed for E-values.
#[allow(clippy::too_many_arguments)]
fn gapped_hsp(
    query: &[u8],
    subject: &[u8],
    anchor: (usize, usize),
    span: (isize, isize),
    band: usize,
    cfg: &SearchConfig,
    db_letters: u64,
    target: &Target,
    subject_id: u32,
) -> Option<Hsp> {
    // The band already restricts the DP to O(query_len * span) cells, so the
    // window is the whole sequence pair: DIAMOND's banded swipe walks the
    // subject to its end instead of clipping a box around the seed, and a
    // clipped window would truncate alignments that trail off the seed.
    // Clipping is kept only as a guard against absurdly large tracebacks.
    let row_cells = (span.1 - span.0).max(0) as usize + 1;
    let budget = 16_000_000usize;
    let max_rows = (budget / row_cells).max(2 * band + 16);
    let mut q0 = 0usize;
    let mut q1 = query.len();
    if q1 - q0 > max_rows {
        q0 = anchor.0.saturating_sub(max_rows / 2);
        q1 = (q0 + max_rows).min(query.len());
        q0 = q1.saturating_sub(max_rows);
    }
    let s0 = 0usize;
    let s1 = subject.len();

    let params = gapped::Params {
        gap_open: cfg.matrix.gap_open,
        gap_extend: cfg.matrix.gap_extend,
        band,
    };
    // The DP uses window-local diagonals, so the span shifts with `q0`.
    let shift = q0 as isize - s0 as isize;
    let aligned = gapped::extend_span(
        &query[q0..q1],
        &subject[s0..s1],
        span.0 - shift,
        span.1 - shift,
        &cfg.matrix,
        &params,
    )?;

    let q_begin = q0 + aligned.q_begin;
    let q_end = q0 + aligned.q_end;
    let s_begin = s0 + aligned.s_begin;
    let s_end = s0 + aligned.s_end;
    // DIAMOND's E-value uses the full sequence lengths, not the HSP spans.
    let evalue = cfg
        .matrix
        .evalue(aligned.score, query.len(), subject.len(), db_letters);
    Some(Hsp {
        query: target.query,
        subject: subject_id as usize,
        frame: target.frame,
        q_begin,
        q_end,
        s_begin,
        s_end,
        score: aligned.score,
        bits: cfg.matrix.bit_score(aligned.score),
        evalue,
        stats: aligned.stats,
        ops: aligned.ops,
    })
}

/// Apply `--max-hsps`/`--max-target-seqs` after merging frames of a query.
pub fn merge_frames(
    results: Vec<Vec<Hsp>>,
    targets: &[Target],
    queries: usize,
    max_hsps: usize,
    max_targets: usize,
) -> Vec<Vec<Hsp>> {
    let mut per_query: Vec<Vec<Hsp>> = (0..queries).map(|_| Vec::new()).collect();
    for (ti, hsps) in results.into_iter().enumerate() {
        if let Some(slot) = per_query.get_mut(targets[ti].query) {
            slot.extend(hsps);
        }
    }
    for hsps in per_query.iter_mut() {
        hsps.sort_by(cmp_hsp);
        let mut order: Vec<usize> = Vec::new();
        let mut by_subject: HashMap<usize, Vec<Hsp>> = HashMap::new();
        for h in hsps.drain(..) {
            match by_subject.get_mut(&h.subject) {
                Some(list) => {
                    if list.len() < max_hsps {
                        list.push(h);
                    }
                }
                None => {
                    if order.len() < max_targets {
                        order.push(h.subject);
                        by_subject.insert(h.subject, vec![h]);
                    }
                }
            }
        }
        for sid in order {
            if let Some(list) = by_subject.remove(&sid) {
                hsps.extend(list);
            }
        }
    }
    per_query
}

/// Validate that the selected matrix exists before the search starts.
pub fn check_matrix(name: &str) -> Result<()> {
    if ScoreMatrix::standard(name, -1, -1).is_err() {
        bail!("Unknown scoring matrix: {name}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_cut_follows_the_sensitivity_preset() {
        let mut cfg = SearchConfig::new(ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap());
        for (name, cut) in [
            ("faster", 0.9),
            ("fast", 0.9),
            ("default", 0.8),
            ("mid-sensitive", 1.0),
            ("sensitive", 1.0),
            ("more-sensitive", 1.0),
            ("very-sensitive", 1.0),
            ("ultra-sensitive", 1.0),
        ] {
            cfg.sensitivity = name.to_string();
            assert_eq!(cfg.seed_cut_for(), cut, "wrong seed cut for --{name}");
        }
    }

    #[test]
    fn explicit_seed_cut_overrides_the_preset() {
        let mut cfg = SearchConfig::new(ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap());
        cfg.sensitivity = "very-sensitive".to_string();
        cfg.seed_cut = Some(0.5);
        assert_eq!(cfg.seed_cut_for(), 0.5);
    }

    #[test]
    fn index_chunks_follow_the_sensitivity_preset() {
        let mut cfg = SearchConfig::new(ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap());
        for (name, n) in [
            ("faster", 4),
            ("fast", 4),
            ("default", 4),
            ("mid-sensitive", 4),
            ("sensitive", 4),
            ("more-sensitive", 4),
            ("very-sensitive", 1),
            ("ultra-sensitive", 1),
        ] {
            cfg.sensitivity = name.to_string();
            assert_eq!(cfg.index_chunks_for(), n, "wrong index chunks for --{name}");
        }
        cfg.index_chunks = Some(7);
        assert_eq!(cfg.index_chunks_for(), 7);
    }
}
