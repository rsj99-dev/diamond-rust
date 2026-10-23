//! Banded affine-gap Smith-Waterman with traceback.
//!
//! The DP runs over a diagonal band centred on the seed hit, like DIAMOND's
//! `Extension` stage: a full quadratic matrix would be far too expensive once
//! a few hundred seeds per query are extended. Gap costs follow BLAST
//! conventions, i.e. a gap of length `L` costs `open + L * extend`.

use crate::hsp::{HspStats, Op};
use crate::matrix::ScoreMatrix;

/// A very negative value that still survives a couple of subtractions.
const NEG: i32 = i32::MIN / 4;

#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub gap_open: i32,
    pub gap_extend: i32,
    /// Half-width of the diagonal band in residues. Used by [`extend`] only;
    /// [`extend_span`] takes the diagonal range directly.
    pub band: usize,
}

/// DIAMOND's band width by query length (`Extension::band`).
///
/// The default sensitivity sets `--extend-mode` to `banded-fast`, so the
/// `BANDED_FAST` row of DIAMOND's table applies, not the wider default row.
pub fn band_for(len: usize) -> usize {
    const THRESHOLD: [usize; 7] = [50, 100, 150, 200, 250, 350, 500];
    const BAND: [usize; 8] = [12, 16, 30, 30, 30, 40, 64, 64];
    let bucket = THRESHOLD.iter().filter(|&&t| len >= t).count();
    BAND[bucket]
}

#[derive(Clone, Debug)]
pub struct Aligned {
    pub score: i32,
    pub q_begin: usize,
    pub q_end: usize,
    pub s_begin: usize,
    pub s_end: usize,
    pub stats: HspStats,
    pub ops: Vec<(Op, u32)>,
}

/// Extend a seed at `anchor` (0-based positions in `query`/`subject`).
///
/// The band is the symmetric interval around the seed diagonal; use
/// [`extend_span`] when several diagonals must share one DP.
/// Returns `None` when no positive-scoring alignment exists in the band.
pub fn extend(
    query: &[u8],
    subject: &[u8],
    anchor: (usize, usize),
    matrix: &ScoreMatrix,
    params: &Params,
) -> Option<Aligned> {
    let diag = anchor.0 as isize - anchor.1 as isize;
    let band = params.band.max(1) as isize;
    extend_span(query, subject, diag - band, diag + band, matrix, params)
}

/// Banded SW restricted to the diagonals `d0..=d1` (`d = query_pos - subject_pos`).
///
/// DIAMOND widens the band of every approximate HSP to
/// `[d_min - band, d_max + band]` and merges the bands while they overlap
/// (`gapped_score.cpp`), so one DP can cover several diagonal segments.
pub fn extend_span(
    query: &[u8],
    subject: &[u8],
    d0: isize,
    d1: isize,
    matrix: &ScoreMatrix,
    params: &Params,
) -> Option<Aligned> {
    let n = query.len();
    let m = subject.len();
    if n == 0 || m == 0 {
        return None;
    }
    let (d0, d1) = (d0.min(d1), d0.max(d1));

    // Per-row subject range [j0, j1] (inclusive).
    let mut j0 = vec![0usize; n + 1];
    let mut j1 = vec![0usize; n + 1];
    let mut off = vec![0usize; n + 2];
    let mut max_w = 1usize;
    for i in 0..=n {
        let lo = (i as isize - d1).clamp(0, m as isize) as usize;
        let hi = (i as isize - d0).clamp(0, m as isize) as usize;
        j0[i] = lo;
        j1[i] = hi.max(lo);
        max_w = max_w.max(j1[i] - j0[i] + 1);
    }
    for i in 0..=n {
        off[i + 1] = off[i] + (j1[i] - j0[i] + 1);
    }
    let mut trace = vec![0u8; off[n + 1]];

    let open_cost = params.gap_open + params.gap_extend;
    let ge = params.gap_extend;
    let mut prev_h = vec![0i32; max_w];
    let mut prev_e = vec![NEG; max_w];
    let mut cur_h = vec![0i32; max_w];
    let mut cur_e = vec![NEG; max_w];

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=n {
        let qi = query[i - 1];
        let a = j0[i];
        let b = j1[i];
        let pa = j0[i - 1];
        let pb = j1[i - 1];
        let mut f_left = NEG;
        for j in a..=b {
            // Diagonal step; H outside the band counts as 0 so alignments can
            // start on the band border.
            let mut d = NEG;
            if j >= 1 {
                let hd = if j > pa && j - 1 <= pb {
                    prev_h[j - 1 - pa]
                } else {
                    0
                };
                d = hd + matrix.score(qi, subject[j - 1]);
            }

            // E: gap in the subject (query consumes a residue).
            let mut e = NEG;
            let mut e_open = 0u8;
            if j >= pa && j <= pb {
                let hp = prev_h[j - pa];
                if hp > NEG / 2 {
                    let v = hp - open_cost;
                    if v > e {
                        e = v;
                    }
                }
                let pe = prev_e[j - pa];
                if pe > NEG / 2 {
                    let v = pe - ge;
                    if v > e {
                        e = v;
                        e_open = 1;
                    }
                }
            }

            // F: gap in the query (subject consumes a residue).
            let mut f = NEG;
            let mut f_open = 0u8;
            if j > a {
                let hl = cur_h[j - 1 - a];
                if hl > NEG / 2 {
                    let v = hl - open_cost;
                    if v > f {
                        f = v;
                    }
                }
                if f_left > NEG / 2 {
                    let v = f_left - ge;
                    if v > f {
                        f = v;
                        f_open = 1;
                    }
                }
            }

            let mut hv = 0i32;
            let mut choice = 0u8;
            if d > hv {
                hv = d;
                choice = 1;
            }
            if e > hv {
                hv = e;
                choice = 2;
            }
            if f > hv {
                hv = f;
                choice = 3;
            }

            cur_h[j - a] = hv;
            cur_e[j - a] = e;
            trace[off[i] + j - a] = choice | (e_open << 2) | (f_open << 3);
            f_left = f;

            if hv > best {
                best = hv;
                best_i = i;
                best_j = j;
            }
        }
        std::mem::swap(&mut prev_h, &mut cur_h);
        std::mem::swap(&mut prev_e, &mut cur_e);
    }

    if best <= 0 {
        return None;
    }

    // Traceback.
    let mut rev: Vec<Op> = Vec::with_capacity(best as usize + 1);
    let mut stats = HspStats::default();
    let mut state = 0u8; // 0 = H, 1 = E (gap in subject), 2 = F (gap in query)
    let mut i = best_i;
    let mut j = best_j;
    let mut in_gap = false;
    while i > 0 && j > 0 {
        if j < j0[i] || j > j1[i] {
            break;
        }
        let t = trace[off[i] + j - j0[i]];
        match state {
            0 => match t & 3 {
                0 => break,
                1 => {
                    stats.add_match(query[i - 1], subject[j - 1], matrix);
                    rev.push(Op::Match);
                    in_gap = false;
                    i -= 1;
                    j -= 1;
                }
                2 => state = 1,
                _ => state = 2,
            },
            1 => {
                stats.length += 1;
                stats.gaps += 1;
                if !in_gap {
                    stats.gap_openings += 1;
                    in_gap = true;
                }
                rev.push(Op::Ins);
                if (t >> 2) & 1 == 0 {
                    state = 0;
                }
                i -= 1;
            }
            _ => {
                stats.length += 1;
                stats.gaps += 1;
                if !in_gap {
                    stats.gap_openings += 1;
                    in_gap = true;
                }
                rev.push(Op::Del);
                if (t >> 3) & 1 == 0 {
                    state = 0;
                }
                j -= 1;
            }
        }
    }

    rev.reverse();
    let mut ops: Vec<(Op, u32)> = Vec::with_capacity(rev.len());
    for op in rev {
        match ops.last_mut() {
            Some(last) if last.0 == op => last.1 += 1,
            _ => ops.push((op, 1)),
        }
    }
    if ops.is_empty() {
        return None;
    }

    Some(Aligned {
        score: best,
        q_begin: i,
        q_end: best_i,
        s_begin: j,
        s_end: best_j,
        stats,
        ops,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap()
    }

    fn params(band: usize) -> Params {
        Params {
            gap_open: 11,
            gap_extend: 1,
            band,
        }
    }

    #[test]
    fn identical_sequences_align_end_to_end() {
        let q = alphabet::encode(b"MKWVTFISLLFLFSSAYSR");
        let score: i32 = (0..q.len()).map(|i| matrix().score(q[i], q[i])).sum();
        let a = extend(&q, &q, (0, 0), &matrix(), &params(30)).unwrap();
        assert_eq!(a.score, score);
        assert_eq!((a.q_begin, a.s_begin), (0, 0));
        assert_eq!((a.q_end, a.s_end), (q.len(), q.len()));
        assert_eq!(a.stats.identities as usize, q.len());
        assert_eq!(a.stats.mismatches, 0);
        assert_eq!(a.stats.gaps, 0);
    }

    #[test]
    fn gaps_are_placed_in_the_shorter_sequence() {
        let q = alphabet::encode(b"MKWVTFISLLFLFSSAYSR");
        // The subject is missing one residue (the second F).
        let s = alphabet::encode(b"MKWVTFISLLLFSSAYSR");
        let a = extend(&q, &s, (0, 0), &matrix(), &params(30)).unwrap();
        assert!(a.stats.gaps > 0, "expected a gap: {:?}", a.stats);
        assert_eq!(a.stats.identities as usize, s.len());
        let perfect: i32 = (0..s.len()).map(|i| matrix().score(s[i], s[i])).sum();
        assert!(a.score < perfect);
    }

    #[test]
    fn unrelated_sequences_give_no_alignment() {
        let q = alphabet::encode(b"WWWWWWWWWW");
        let s = alphabet::encode(b"AAAAAAAAAA");
        assert!(extend(&q, &s, (0, 0), &matrix(), &params(15)).is_none());
    }

    #[test]
    fn band_widths_match_diamond() {
        assert_eq!(band_for(10), 12);
        assert_eq!(band_for(50), 16);
        assert_eq!(band_for(120), 30);
        assert_eq!(band_for(400), 64);
        assert_eq!(band_for(900), 64);
    }
}
