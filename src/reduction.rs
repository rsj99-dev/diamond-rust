//! Reduced amino acid alphabets used to hash seeds.
//!
//! DIAMOND does not hash 20-letter seeds: it first folds residues into a
//! smaller alphabet (`murphy10` by default), which shrinks the seed space and
//! makes a seed match tolerant to conservative substitutions. The groups are
//! written exactly like DIAMOND's `Reduction` definition strings, e.g.
//! `"A KR EDNQ C G H ILVM FYW P ST"` (10 groups).

use crate::alphabet;

#[derive(Clone, Debug)]
pub struct Reduction {
    map: [u8; alphabet::AA_COUNT],
    size: u64,
}

impl Reduction {
    /// Build a reduction from a DIAMOND-style group definition.
    pub fn parse(def: &str) -> Reduction {
        let mut map = [u8::MAX; alphabet::AA_COUNT];
        let mut bucket = 0u8;
        for group in def.split_whitespace() {
            for &c in group.as_bytes() {
                let letter = alphabet::encode_byte(c);
                if alphabet::is_canonical(letter) {
                    map[letter as usize] = bucket;
                }
            }
            bucket += 1;
        }
        assert!(bucket > 0, "empty reduction definition");
        Reduction {
            map,
            size: bucket as u64,
        }
    }

    #[inline]
    pub fn size(&self) -> u64 {
        self.size
    }

    #[inline]
    pub fn bucket(&self, letter: u8) -> u8 {
        self.map[(letter as usize).min(alphabet::AA_COUNT - 1)]
    }

    pub fn murphy10() -> Reduction {
        Reduction::parse("A KR EDNQ C G H ILVM FYW P ST")
    }

    pub fn diamond9() -> Reduction {
        Reduction::parse("A G P CH ST NDE RQK FWY ILMV")
    }

    pub fn diamond14() -> Reduction {
        Reduction::parse("A C Q E G H M P S T ND RK FWY ILV")
    }

    pub fn diamond18() -> Reduction {
        Reduction::parse("A R N D C Q E G H L K M F P S T WY IV")
    }

    /// No reduction: one bucket per canonical amino acid.
    pub fn none() -> Reduction {
        Reduction::parse("A R N D C Q E G H I L K M F P S T W Y V")
    }

    /// Look up the reduction used by DIAMOND for a given sensitivity name.
    ///
    /// Only `--very-sensitive` folds to a smaller alphabet (`diamond9`); every
    /// other preset, including `--ultra-sensitive`, hashes `murphy10`
    /// (`Search::sensitivity_traits`).
    pub fn for_sensitivity(sensitivity: &str) -> Reduction {
        match sensitivity {
            "very-sensitive" => Reduction::diamond9(),
            _ => Reduction::murphy10(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn murphy10_groups() {
        let r = Reduction::murphy10();
        assert_eq!(r.size(), 10);
        assert_eq!(r.bucket(alphabet::encode_byte(b'A')), 0);
        assert_eq!(r.bucket(alphabet::encode_byte(b'K')), 1);
        assert_eq!(r.bucket(alphabet::encode_byte(b'R')), 1);
        // L, I, V, M share the last group.
        assert_eq!(
            r.bucket(alphabet::encode_byte(b'L')),
            r.bucket(alphabet::encode_byte(b'M'))
        );
    }

    #[test]
    fn no_reduction_is_bijective() {
        let r = Reduction::none();
        assert_eq!(r.size(), 20);
        for i in 0..20u8 {
            assert_eq!(r.bucket(i), i);
        }
    }

    #[test]
    fn sensitivity_presets_pick_diamonds_reduction() {
        // `Search::sensitivity_traits`: only --very-sensitive folds to 9 groups.
        assert_eq!(Reduction::for_sensitivity("default").size(), 10);
        assert_eq!(Reduction::for_sensitivity("sensitive").size(), 10);
        assert_eq!(Reduction::for_sensitivity("very-sensitive").size(), 9);
        assert_eq!(Reduction::for_sensitivity("ultra-sensitive").size(), 10);
    }
}
