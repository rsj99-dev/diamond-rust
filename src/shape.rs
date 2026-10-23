//! Seed shapes and seed hashing.
//!
//! A shape is a string of `0`/`1`; the `1` positions are read from the
//! sequence and folded through the reduction into a base-`size` number. This
//! is exactly DIAMOND's `Shape::set_seed`, so a shape such as `111101110111`
//! (weight 10, span 12) has the same seeds here as in the C++ implementation.

use crate::reduction::Reduction;

#[derive(Clone, Debug)]
pub struct Shape {
    /// The `0`/`1` pattern as given on the command line or in a preset.
    pub code: String,
    /// Number of `1`s (residues that make up the seed).
    pub weight: u32,
    /// Number of residues spanned by the shape.
    pub span: u32,
    /// Offsets of the `1` positions.
    pub positions: Vec<u32>,
}

impl Shape {
    pub fn parse(code: &str) -> Shape {
        let mut positions = Vec::new();
        let mut span = 0u32;
        for (i, c) in code.bytes().enumerate() {
            match c {
                b'1' => {
                    positions.push(i as u32);
                    span = i as u32 + 1;
                }
                b'0' => span = i as u32 + 1,
                _ => panic!("invalid seed shape {:?}: expected only 0/1", code),
            }
        }
        assert!(!positions.is_empty(), "seed shape {:?} has weight 0", code);
        Shape {
            code: code.to_string(),
            weight: positions.len() as u32,
            span,
            positions,
        }
    }

    /// Contiguous shape of the given weight, e.g. `k=5` -> `11111`.
    pub fn contiguous(weight: u32) -> Shape {
        Shape::parse(&"1".repeat(weight as usize))
    }

    /// Maximum key value + 1 for this shape under `reduction`.
    pub fn key_space(&self, reduction: &Reduction) -> u64 {
        reduction.size().saturating_pow(self.weight)
    }

    /// Hash the seed starting at `pos`.
    ///
    /// Returns `None` when the seed would leave the sequence or contains a
    /// residue that is not a canonical amino acid.
    #[inline]
    pub fn packed(&self, seq: &[u8], pos: usize, reduction: &Reduction) -> Option<u64> {
        if pos + self.span as usize > seq.len() {
            return None;
        }
        let size = reduction.size();
        let mut key = 0u64;
        for &off in &self.positions {
            let letter = seq[pos + off as usize];
            if letter >= 20 {
                return None;
            }
            key = key * size + reduction.bucket(letter) as u64;
        }
        Some(key)
    }
}

/// DIAMOND's preset shapes for each sensitivity level (`Search::shape_codes`).
pub fn shapes_for(sensitivity: &str, limit: usize) -> Vec<Shape> {
    let codes: &[&str] = match sensitivity {
        "faster" | "fast" => &["1101110101101111"],
        "default" => &["111101110111", "111011010010111"],
        "mid-sensitive" => &[
            "11110110111",
            "1101100111101",
            "1110010101111",
            "11010101100111",
            "11101110001011",
            "1110100100010111",
            "1101000011010111",
            "1110011000011011",
        ],
        "sensitive" => &[
            "1011110111",
            "110100100010111",
            "11001011111",
            "101110001111",
            "11011101100001",
            "1111010010101",
            "111001001001011",
            "10101001101011",
            "111101010011",
            "1111000010000111",
            "1100011011011",
            "1101010000011011",
            "1110001010101001",
            "110011000110011",
            "11011010001101",
            "1101001100010011",
        ],
        "more-sensitive" => &[
            "1011110111",
            "110100100010111",
            "11001011111",
            "101110001111",
            "11011101100001",
            "1111010010101",
            "111001001001011",
            "10101001101011",
            "111101010011",
            "1111000010000111",
            "1100011011011",
            "1101010000011011",
            "1110001010101001",
            "110011000110011",
            "11011010001101",
            "1101001100010011",
        ],
        "very-sensitive" => &[
            "11110111",
            "11100100111",
            "110010101011",
            "11010001001011",
            "10101100001101",
            "110100100100011",
            "1010010100010011",
            "1100101000001011",
            "11100000100010101",
            "11000100010010011",
            "11010000001000111",
            "110001001000010011",
            "1010001000100001011",
            "1100010100000010011",
            "1100100000101000011",
            "1101000010000001011",
        ],
        "ultra-sensitive" => &[
            "1111111",
            "11101111",
            "110011111",
            "110110111",
            "111111001",
            "1010111011",
            "1011110101",
            "1111000111",
            "10011110011",
            "10101101101",
            "10111010101",
            "11001010111",
            "11001100111",
            "11010101101",
            "11110001011",
            "100111010011",
            "101100110101",
            "101110000111",
            "110100101011",
            "110110001101",
            "111000110011",
            "1010001011011",
            "1010101000111",
            "1010110100011",
            "1100100110011",
            "1100101001011",
            "1101001100101",
            "1101010101001",
            "1110001010101",
            "1110010010011",
            "10100001101101",
            "11000100010111",
            "11010000100111",
            "11010100110001",
            "11101000011001",
            "11110000001101",
            "11110100000011",
            "101001000001111",
            "110000100101011",
            "110010010000111",
            "110101100001001",
            "110110000010011",
            "111001000100011",
            "111100000100101",
            "1000110010010101",
            "1001000100101101",
            "1001000110011001",
            "1010001001001011",
            "1010001010010011",
            "1010010001010101",
            "1010010100010011",
            "1010010101001001",
            "1010100000101011",
            "1010100011000101",
            "1011000010001011",
            "1100010000111001",
            "1100010010001011",
            "1100100001001011",
            "1100100100100011",
            "1100110000001101",
            "1101000100010011",
            "1101000110000101",
            "1110000001010011",
            "1110100000010101",
        ],
        _ => &["111101110111", "111011010010111"],
    };
    let n = if limit == 0 {
        codes.len()
    } else {
        limit.min(codes.len())
    };
    codes[..n].iter().map(|c| Shape::parse(c)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet;

    fn encoded(s: &str) -> Vec<u8> {
        alphabet::encode(s.as_bytes())
    }

    #[test]
    fn contiguous_seed_key() {
        let shape = Shape::contiguous(3);
        let red = Reduction::none();
        let seq = encoded("MKWV");
        // M=K index 12, K=11, W=17 in the 20-bucket reduction.
        let k = shape.packed(&seq, 0, &red).unwrap();
        assert_eq!(k, 12 * 400 + 11 * 20 + 17);
        assert!(
            shape.packed(&seq, 2, &red).is_none(),
            "seed leaves the sequence"
        );
    }

    #[test]
    fn spaced_seed_uses_only_ones() {
        let shape = Shape::parse("101");
        assert_eq!(shape.weight, 2);
        assert_eq!(shape.span, 3);
        let red = Reduction::none();
        let seq = encoded("ACD");
        // A and D only, C is skipped: A fills the high digit, D the low one.
        assert_eq!(shape.packed(&seq, 0, &red).unwrap(), 3);
    }

    #[test]
    fn reduced_seed_merges_similar_residues() {
        let shape = Shape::contiguous(2);
        let red = Reduction::murphy10();
        // L and I are in the same bucket, so LI and LL hash identically.
        assert_eq!(
            shape.packed(&encoded("LI"), 0, &red),
            shape.packed(&encoded("LL"), 0, &red)
        );
    }

    #[test]
    fn presets_match_diamond() {
        let s = shapes_for("default", 0);
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].code, "111101110111");
        assert_eq!(s[1].code, "111011010010111");
        assert_eq!(s[0].weight, 10);
        assert_eq!(s[1].weight, 10);
        assert_eq!(shapes_for("sensitive", 4).len(), 4);
    }

    #[test]
    fn more_sensitive_reuses_the_sensitive_shapes() {
        // DIAMOND's `shape_codes` lists the same 16x8 iedera set for both.
        let sensitive: Vec<String> = shapes_for("sensitive", 0)
            .iter()
            .map(|s| s.code.clone())
            .collect();
        let more: Vec<String> = shapes_for("more-sensitive", 0)
            .iter()
            .map(|s| s.code.clone())
            .collect();
        assert_eq!(sensitive, more);
        assert_eq!(sensitive.len(), 16);
    }

    #[test]
    fn very_sensitive_uses_the_14x7_shapes() {
        let s = shapes_for("very-sensitive", 0);
        assert_eq!(s.len(), 16);
        assert_eq!(s[0].code, "11110111");
        assert_eq!(s[0].weight, 7);
    }
}
