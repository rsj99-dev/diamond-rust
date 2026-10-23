//! HSP records shared by the alignment stages and the output layer.

use crate::alphabet;
use crate::matrix::ScoreMatrix;

/// One column of a pairwise alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// Aligned residues (match or substitution).
    Match,
    /// Gap in the subject: the query consumes a residue, the subject does not.
    Ins,
    /// Gap in the query: the subject consumes a residue, the query does not.
    Del,
}

impl Op {
    pub fn code(self) -> char {
        match self {
            Op::Match => 'M',
            Op::Ins => 'I',
            Op::Del => 'D',
        }
    }
}

/// Column counts of an alignment.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HspStats {
    /// Identical residue pairs.
    pub identities: u32,
    /// Residue pairs with a positive substitution score.
    pub positives: u32,
    /// Aligned columns with different residues (gaps excluded).
    pub mismatches: u32,
    /// Number of gap runs.
    pub gap_openings: u32,
    /// Total number of gap residues.
    pub gaps: u32,
    /// Number of alignment columns (residues + gaps).
    pub length: u32,
}

impl HspStats {
    #[inline]
    pub fn add_match(&mut self, a: u8, b: u8, matrix: &ScoreMatrix) {
        self.length += 1;
        if a == b {
            self.identities += 1;
        } else {
            self.mismatches += 1;
        }
        if matrix.score(a, b) > 0 {
            self.positives += 1;
        }
    }
}

/// A high-scoring segment pair.
///
/// `q_begin..q_end` and `s_begin..s_end` are 0-based half-open coordinates.
/// For translated searches they refer to the *translated* frame; the
/// nucleotide coordinates are derived with [`Hsp::nt_range`].
#[derive(Clone, Debug)]
pub struct Hsp {
    /// Index of the query record (not the frame).
    pub query: usize,
    /// Database sequence index.
    pub subject: usize,
    /// Frame index: 0 for `blastp`, 0..5 for `blastx`.
    pub frame: u8,
    pub q_begin: usize,
    pub q_end: usize,
    pub s_begin: usize,
    pub s_end: usize,
    pub score: i32,
    pub bits: f64,
    pub evalue: f64,
    pub stats: HspStats,
    /// Run-length encoded alignment, left to right.
    pub ops: Vec<(Op, u32)>,
}

impl Hsp {
    /// `blastx` qframe column: +1..+3 for forward frames, -1..-3 for reverse.
    pub fn q_frame(&self) -> i32 {
        if self.frame < 3 {
            self.frame as i32 + 1
        } else {
            2 - self.frame as i32
        }
    }

    pub fn is_reverse(&self) -> bool {
        self.frame >= 3
    }

    /// Nucleotide range on the forward strand (0-based, half-open).
    pub fn nt_range(&self, dna_len: usize) -> (usize, usize) {
        let f = self.frame as usize;
        if f < 3 {
            (f + 3 * self.q_begin, f + 3 * self.q_end)
        } else {
            let o = f - 3;
            (
                dna_len - (o + 3 * self.q_end),
                dna_len - (o + 3 * self.q_begin),
            )
        }
    }

    /// Percentage of identical residues over the alignment length.
    pub fn pident(&self) -> f64 {
        if self.stats.length == 0 {
            0.0
        } else {
            self.stats.identities as f64 * 100.0 / self.stats.length as f64
        }
    }

    /// Percentage of positive-scoring residue pairs.
    pub fn ppos(&self) -> f64 {
        if self.stats.length == 0 {
            0.0
        } else {
            self.stats.positives as f64 * 100.0 / self.stats.length as f64
        }
    }

    /// Query coverage of this HSP in percent.
    ///
    /// For translated searches `query_len` is the length of the translated
    /// frame; the caller passes the amino-acid length.
    pub fn qcovhsp(&self, query_len: usize) -> f64 {
        if query_len == 0 {
            0.0
        } else {
            (self.q_end - self.q_begin) as f64 * 100.0 / query_len as f64
        }
    }

    /// Subject coverage of this HSP in percent.
    pub fn scovhsp(&self, subject_len: usize) -> f64 {
        if subject_len == 0 {
            0.0
        } else {
            (self.s_end - self.s_begin) as f64 * 100.0 / subject_len as f64
        }
    }

    /// `M12I3M40` style CIGAR string (SAM conventions, gaps are runs).
    pub fn cigar(&self) -> String {
        let mut out = String::new();
        for &(op, n) in &self.ops {
            out.push_str(&n.to_string());
            out.push(op.code());
        }
        out
    }

    /// Gapped query and subject strings for the aligned region.
    ///
    /// `query`/`subject` are the full encoded sequences the coordinates refer to.
    pub fn aligned_strings(&self, query: &[u8], subject: &[u8]) -> (String, String) {
        let mut q = String::with_capacity(self.stats.length as usize);
        let mut s = String::with_capacity(self.stats.length as usize);
        let mut qi = self.q_begin;
        let mut si = self.s_begin;
        for &(op, n) in &self.ops {
            for _ in 0..n {
                match op {
                    Op::Match => {
                        q.push(alphabet::decode_byte(query[qi]) as char);
                        s.push(alphabet::decode_byte(subject[si]) as char);
                        qi += 1;
                        si += 1;
                    }
                    Op::Ins => {
                        q.push(alphabet::decode_byte(query[qi]) as char);
                        s.push('-');
                        qi += 1;
                    }
                    Op::Del => {
                        q.push('-');
                        s.push(alphabet::decode_byte(subject[si]) as char);
                        si += 1;
                    }
                }
            }
        }
        (q, s)
    }

    /// BLAST traceback operations (`BTOP`) for the aligned region.
    pub fn btop(&self, query: &[u8], subject: &[u8]) -> String {
        let mut out = String::new();
        let mut matches = 0u32;
        let mut qi = self.q_begin;
        let mut si = self.s_begin;
        let flush = |out: &mut String, matches: &mut u32| {
            if *matches > 0 {
                out.push_str(&matches.to_string());
                *matches = 0;
            }
        };
        for &(op, n) in &self.ops {
            for _ in 0..n {
                match op {
                    Op::Match => {
                        let a = alphabet::decode_byte(query[qi]);
                        let b = alphabet::decode_byte(subject[si]);
                        if a == b {
                            matches += 1;
                        } else {
                            flush(&mut out, &mut matches);
                            out.push(a as char);
                            out.push(b as char);
                        }
                        qi += 1;
                        si += 1;
                    }
                    Op::Ins => {
                        // An insertion is an extra query residue: query char
                        // followed by a gap in the subject.
                        flush(&mut out, &mut matches);
                        out.push(alphabet::decode_byte(query[qi]) as char);
                        out.push('-');
                        qi += 1;
                    }
                    Op::Del => {
                        flush(&mut out, &mut matches);
                        out.push('-');
                        out.push(alphabet::decode_byte(subject[si]) as char);
                        si += 1;
                    }
                }
            }
        }
        flush(&mut out, &mut matches);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_follow_diamond_numbering() {
        for (frame, expect) in [(0, 1), (1, 2), (2, 3), (3, -1), (4, -2), (5, -3)] {
            let h = Hsp {
                query: 0,
                subject: 0,
                frame,
                q_begin: 0,
                q_end: 0,
                s_begin: 0,
                s_end: 0,
                score: 0,
                bits: 0.0,
                evalue: 0.0,
                stats: HspStats::default(),
                ops: Vec::new(),
            };
            assert_eq!(h.q_frame(), expect);
        }
    }

    #[test]
    fn reverse_frame_maps_back_to_forward_strand() {
        let h = Hsp {
            query: 0,
            subject: 0,
            frame: 3, // -1 frame, offset 0
            q_begin: 2,
            q_end: 5,
            s_begin: 0,
            s_end: 0,
            score: 0,
            bits: 0.0,
            evalue: 0.0,
            stats: HspStats::default(),
            ops: Vec::new(),
        };
        // 30 nt -> aa 2..5 of the reverse complement covers nt 15..24.
        assert_eq!(h.nt_range(30), (30 - 15, 30 - 6));
    }
}
