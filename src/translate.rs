//! Six-frame translation for `blastx`.
//!
//! Frames 0..2 are the forward strand with offsets 0..2, frames 3..5 the
//! reverse strand, matching DIAMOND's numbering (`frame < 3` = forward).

use crate::alphabet;

/// Codon table in TCAG order (index = 16*b0 + 4*b1 + b2).
const CODON: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const STOP: u8 = alphabet::STOP_LETTER;

#[inline]
fn base(b: u8) -> u8 {
    match b {
        b'T' | b't' | b'U' | b'u' => 0,
        b'C' | b'c' => 1,
        b'A' | b'a' => 2,
        b'G' | b'g' => 3,
        _ => 4,
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

/// Reverse complement of a nucleotide sequence (encoded as ASCII).
pub fn reverse_complement(nt: &[u8]) -> Vec<u8> {
    nt.iter().rev().map(|&b| complement(b)).collect()
}

/// Translate one frame starting at `offset` (0, 1 or 2) of `nt`.
///
/// Codons containing ambiguous bases become `X`; stop codons become `*`.
pub fn translate_frame(nt: &[u8], offset: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(nt.len().saturating_sub(offset) / 3);
    let mut i = offset;
    while i + 3 <= nt.len() {
        let a = base(nt[i]);
        let b = base(nt[i + 1]);
        let c = base(nt[i + 2]);
        if a > 3 || b > 3 || c > 3 {
            out.push(alphabet::X_LETTER);
        } else {
            out.push(alphabet::encode_byte(CODON[(a * 16 + b * 4 + c) as usize]));
        }
        i += 3;
    }
    out
}

/// All six reading frames, in DIAMOND's frame order.
pub fn six_frames(nt: &[u8]) -> Vec<Vec<u8>> {
    let rc = reverse_complement(nt);
    let mut frames = Vec::with_capacity(6);
    for offset in 0..3 {
        frames.push(translate_frame(nt, offset));
    }
    for offset in 0..3 {
        frames.push(translate_frame(&rc, offset));
    }
    frames
}

/// Length of the longest stop-free stretch (DIAMOND's `--min-orf` filter).
pub fn longest_orf(frame: &[u8]) -> usize {
    let mut best = 0usize;
    let mut run = 0usize;
    for &l in frame {
        if l == STOP {
            run = 0;
        } else {
            run += 1;
            best = best.max(run);
        }
    }
    best
}

/// Mask every stop-free stretch shorter than `min_len` with `X`, mirroring
/// DIAMOND's `find_orfs()` (short ORFs cannot seed or score).
pub fn mask_short_orfs(frame: &mut [u8], min_len: usize) {
    if min_len <= 1 {
        return;
    }
    let mut begin = 0usize;
    for i in 0..frame.len() {
        if frame[i] == STOP {
            if i - begin < min_len {
                for letter in &mut frame[begin..i] {
                    *letter = alphabet::X_LETTER;
                }
            }
            begin = i + 1;
        }
    }
    if frame.len() - begin < min_len {
        for letter in &mut frame[begin..] {
            *letter = alphabet::X_LETTER;
        }
    }
}

/// DIAMOND's default minimum ORF length for `blastx`.
pub fn min_orf_len(dna_len: usize, user: Option<usize>) -> usize {
    match user {
        Some(n) => n,
        None => {
            if dna_len < 30 {
                1
            } else if dna_len < 100 {
                20
            } else {
                40
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aa(s: &str) -> Vec<u8> {
        alphabet::encode(s.as_bytes())
    }

    fn nt(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    #[test]
    fn standard_code_is_used() {
        assert_eq!(translate_frame(&nt("ATGGCCTAA"), 0), aa("MA*"));
        assert_eq!(translate_frame(&nt("TTTTTCCTT"), 0), aa("FFL"));
        assert_eq!(translate_frame(&nt("atgaaa"), 0), aa("MK"));
    }

    #[test]
    fn ambiguous_codons_become_x() {
        assert_eq!(translate_frame(&nt("ATGNNTAAA"), 0), aa("MXK"));
        assert_eq!(translate_frame(&nt("ATGTAA"), 0), aa("M*"));
    }

    #[test]
    fn reverse_frames_read_the_complement() {
        let frames = six_frames(&nt("ATGAAATTT"));
        assert_eq!(frames[0], aa("MKF"));
        // Reverse complement of ATGAAATTT is AAATTTCAT -> KFH.
        assert_eq!(frames[3], aa("KFH"));
    }

    #[test]
    fn orf_lengths_are_measured_between_stops() {
        assert_eq!(longest_orf(&aa("MK*MKF*")), 3);
        assert_eq!(longest_orf(&aa("MKF")), 3);
    }

    #[test]
    fn short_orfs_are_masked() {
        let mut frame = aa("MK*MKF*WWWW");
        mask_short_orfs(&mut frame, 3);
        assert_eq!(frame, aa("XX*MKF*WWWW"));
    }
}
