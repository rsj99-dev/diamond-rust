//! Ungapped X-drop extension on a seed diagonal.
//!
//! This is a direct port of DIAMOND's `xdrop_ungapped()`
//! (`src/dp/ungapped_align.cpp`): the extension starts at the seed, walks left
//! and right as long as the score stays within `xdrop` of the best score seen
//! so far, and returns the best-scoring segment.

use crate::matrix::ScoreMatrix;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Segment {
    pub q_start: usize,
    pub s_start: usize,
    pub len: usize,
    pub score: i32,
}

impl Segment {
    #[inline]
    pub fn diag(&self) -> isize {
        self.q_start as isize - self.s_start as isize
    }

    #[inline]
    pub fn q_end(&self) -> usize {
        self.q_start + self.len
    }

    #[inline]
    pub fn s_end(&self) -> usize {
        self.s_start + self.len
    }
}

/// Extend the seed whose first residue is `query[qa]`/`subject[sa]`.
///
/// `xdrop` is in raw score units (`raw_ungapped_xdrop` in DIAMOND).
pub fn xdrop(
    query: &[u8],
    subject: &[u8],
    qa: usize,
    sa: usize,
    matrix: &ScoreMatrix,
    xdrop: i32,
) -> Segment {
    let mut score = 0i32;
    let mut st = 0i32;
    let mut n = 1usize;
    let mut delta = 0usize;
    let mut len = 0usize;

    // Left extension.
    let mut q = qa as isize - 1;
    let mut s = sa as isize - 1;
    while score - st < xdrop && q >= 0 && s >= 0 {
        st += matrix.score(query[q as usize], subject[s as usize]);
        if st > score {
            score = st;
            delta = n;
        }
        q -= 1;
        s -= 1;
        n += 1;
    }

    // Right extension, continuing from the left-hand best score.
    q = qa as isize;
    s = sa as isize;
    st = score;
    n = 1;
    while score - st < xdrop && (q as usize) < query.len() && (s as usize) < subject.len() {
        st += matrix.score(query[q as usize], subject[s as usize]);
        if st > score {
            score = st;
            len = n;
        }
        q += 1;
        s += 1;
        n += 1;
    }

    Segment {
        q_start: qa - delta,
        s_start: sa - delta,
        len: len + delta,
        score,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::standard("BLOSUM62", -1, -1).unwrap()
    }

    #[test]
    fn identical_region_extends_over_the_whole_seed() {
        let q = alphabet::encode(b"MKWVTFISLLFLFSSAYSR");
        let s = alphabet::encode(b"MKWVTFISLLFLFSSAYSR");
        let seg = xdrop(&q, &s, 0, 0, &matrix(), 40);
        assert_eq!(seg.q_start, 0);
        assert_eq!(seg.len, q.len());
        let expect: i32 = (0..q.len()).map(|i| matrix().score(q[i], q[i])).sum();
        assert_eq!(seg.score, expect);
    }

    #[test]
    fn stops_at_the_xdrop_boundary() {
        let q = alphabet::encode(b"WWWWWWWWWWCCCCCCCCCC");
        let s = alphabet::encode(b"WWWWWWWWWWDDDDDDDDDD");
        let seg = xdrop(&q, &s, 0, 0, &matrix(), 5);
        // Each C/D pair costs 3; with xdrop 5 the extension is cut two
        // columns after the last improvement, so the reported segment ends
        // at the poly-W run.
        assert_eq!(seg.len, 10);
        assert_eq!(seg.score, 110);
    }

    #[test]
    fn coordinates_are_clipped_at_the_sequence_start() {
        let q = alphabet::encode(b"AAAAAA");
        let s = alphabet::encode(b"AAAAAA");
        let seg = xdrop(&q, &s, 0, 0, &matrix(), 30);
        assert_eq!((seg.q_start, seg.s_start), (0, 0));
        assert_eq!(seg.len, 6);
    }
}
