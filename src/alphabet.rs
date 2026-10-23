//! Residue encoding.
//!
//! DIAMOND numbers amino acids in NCBI order so that a scoring matrix can be
//! indexed with `a * 25 + b` without a translation step:
//!
//! ```text
//! A R N D C Q E G H I L K M F P S T W Y V  B  J  Z  X  *
//! 0 1 2 3 4 5 6 7 8 9 ...                20 21 22 23 24
//! ```
//!
//! Only `A`..`V` (0..20) may take part in seeds; `B`/`J`/`Z` are accepted by
//! the dynamic programming but excluded from seeding, and `X`/`*` are the
//! ambiguous / stop symbols.

pub const AA_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";
pub const AA_COUNT: usize = 25;
/// Number of unambiguous amino acids.
pub const TRUE_AA: usize = 20;
/// Ambiguous residue (`X`).
pub const X_LETTER: u8 = 23;
/// Stop codon (`*`).
pub const STOP_LETTER: u8 = 24;

const fn build_table() -> [u8; 256] {
    let mut t = [X_LETTER; 256];
    let symbols = b"ARNDCQEGHILKMFPSTWYVBJZX*";
    let mut i = 0;
    while i < 25 {
        let c = symbols[i];
        t[c as usize] = i as u8;
        // Lower-case mirror.
        if c.is_ascii_uppercase() {
            t[c.to_ascii_lowercase() as usize] = i as u8;
        }
        i += 1;
    }
    // Selenocysteine / pyrrolysine are scored as their closest relatives,
    // matching DIAMOND's `UO-` ignore list (which maps them to X otherwise).
    t[b'U' as usize] = 4; // C
    t[b'u' as usize] = 4;
    t[b'O' as usize] = 11; // K
    t[b'o' as usize] = 11;
    t
}

static TABLE: [u8; 256] = build_table();

#[inline]
pub fn encode_byte(c: u8) -> u8 {
    TABLE[c as usize]
}

/// Upper-case and strip layout characters (`-`, `.`, whitespace).
pub fn encode(raw: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len());
    for &c in raw {
        if c.is_ascii_whitespace() || c == b'-' || c == b'.' {
            continue;
        }
        out.push(TABLE[c as usize]);
    }
    out
}

#[inline]
pub fn decode_byte(l: u8) -> u8 {
    AA_SYMBOLS[(l as usize).min(AA_COUNT - 1)]
}

pub fn decode(letters: &[u8]) -> String {
    letters.iter().map(|&l| decode_byte(l) as char).collect()
}

/// An unambiguous amino acid (`A`..`V`).
#[inline]
pub fn is_canonical(l: u8) -> bool {
    l < TRUE_AA as u8
}

/// A letter that can be part of a seed (`A`..`Z`, i.e. everything but `X`/`*`).
#[inline]
pub fn is_amino_acid(l: u8) -> bool {
    l < X_LETTER
}

/// First whitespace-delimited token of a FASTA header.
pub fn seq_id(header: &str) -> &str {
    header
        .trim_start_matches(['>', ' '])
        .split(|c: char| c.is_ascii_whitespace())
        .next()
        .unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ncbi_order_roundtrip() {
        for (i, &c) in AA_SYMBOLS.iter().enumerate() {
            assert_eq!(encode_byte(c), i as u8, "{}", c as char);
        }
        assert_eq!(encode_byte(b'a'), 0);
        assert_eq!(encode_byte(b'U'), 4);
        assert_eq!(encode_byte(b'?'), X_LETTER);
    }

    #[test]
    fn strips_layout_characters() {
        assert_eq!(encode(b"AC-DE\nFG"), vec![0, 4, 3, 6, 13, 7]);
    }

    #[test]
    fn header_id() {
        assert_eq!(seq_id(">sp|P1|ABC description here"), "sp|P1|ABC");
        assert_eq!(seq_id(">  lone"), "lone");
    }
}
