//! Seed enumeration and the double-index hash join.
//!
//! Seeds are hashed from a reduced alphabet and a spaced shape (see
//! [`crate::shape`]); a hit is a pair of equal seed keys on the query and the
//! reference side. Both arrays are sorted by key and merged, which is the
//! `--algo 0` (double-indexed) strategy of DIAMOND.

use crate::reduction::Reduction;
use crate::shape::Shape;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug)]
pub struct SeedEntry {
    pub key: u64,
    pub seq: u32,
    pub pos: u32,
}

#[derive(Clone, Copy, Debug)]
pub struct SeedHit {
    pub query: u32,
    pub qpos: u32,
    pub subject: u32,
    pub spos: u32,
}

fn ln_factorial() -> &'static [f64; 33] {
    static TABLE: OnceLock<[f64; 33]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [0.0f64; 33];
        for n in 2..33 {
            t[n] = t[n - 1] + (n as f64).ln();
        }
        t
    })
}

/// DIAMOND's `Search::seed_is_complex`: the log multinomial coefficient of the
/// reduced seed must reach `cut`.
///
/// `cut` is `cfg.seed_complexity_cut`, i.e. `seed_cut * ln2 * shapes[0].weight`
/// (see `setup.cpp`), not the raw `--seed-cut` value: the preset's 0.8/1.0 is
/// "extra bits per seed position" and has to be scaled to the shape weight.
pub fn seed_is_complex(seq: &[u8], shape: &Shape, reduction: &Reduction, cut: f64) -> bool {
    if cut <= 0.0 {
        return true;
    }
    let mut count = [0u32; 32];
    for &off in &shape.positions {
        let letter = seq[off as usize];
        if letter >= 20 {
            return false;
        }
        let bucket = reduction.bucket(letter) as usize;
        if bucket >= count.len() {
            return false;
        }
        count[bucket] += 1;
    }
    let lnfact = ln_factorial();
    let mut entropy = lnfact[shape.weight as usize];
    for &c in count.iter() {
        if c > 0 {
            entropy -= lnfact[c as usize];
        }
    }
    entropy >= cut
}

/// Hash every seed of every sequence and return them sorted by key.
///
/// `base` is added to the sequence index so that seeds of a database chunk
/// keep global numbering. `complexity_cut` is `cfg.seed_complexity_cut`.
pub fn build(
    seqs: &[&[u8]],
    shape: &Shape,
    reduction: &Reduction,
    base: u32,
    complexity_cut: f64,
) -> Vec<SeedEntry> {
    let mut out = Vec::new();
    for (i, seq) in seqs.iter().enumerate() {
        if seq.len() < shape.span as usize {
            continue;
        }
        for pos in 0..=seq.len() - shape.span as usize {
            let Some(key) = shape.packed(seq, pos, reduction) else {
                continue;
            };
            if !seed_is_complex(&seq[pos..], shape, reduction, complexity_cut) {
                continue;
            }
            out.push(SeedEntry {
                key,
                seq: base + i as u32,
                pos: pos as u32,
            });
        }
    }
    out.sort_unstable_by_key(|e| e.key);
    out
}

/// Merge-join the query and reference seed arrays.
///
/// `max_expand` caps the number of hits produced by a single key, which keeps
/// low-complexity repeats from exploding the hit list.
pub fn join_into(
    query: &[SeedEntry],
    reference: &[SeedEntry],
    out: &mut [HashMap<u32, Vec<SeedHit>>],
    max_expand: usize,
) {
    let (mut i, mut j) = (0usize, 0usize);
    while i < query.len() && j < reference.len() {
        match query[i].key.cmp(&reference[j].key) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let key = query[i].key;
                let qi = i + query[i..].iter().take_while(|e| e.key == key).count();
                let rj = j + reference[j..].iter().take_while(|e| e.key == key).count();
                if (qi - i) * (rj - j) <= max_expand {
                    for q in &query[i..qi] {
                        let bucket = &mut out[q.seq as usize];
                        for r in &reference[j..rj] {
                            bucket.entry(r.seq).or_default().push(SeedHit {
                                query: q.seq,
                                qpos: q.pos,
                                subject: r.seq,
                                spos: r.pos,
                            });
                        }
                    }
                }
                i = qi;
                j = rj;
            }
        }
    }
}

/// Split `n` items into at most `parts` balanced, contiguous chunks.
pub fn chunks(n: usize, parts: usize) -> Vec<(usize, usize)> {
    let parts = parts.max(1);
    if n == 0 {
        return Vec::new();
    }
    let parts = parts.min(n);
    let mut out = Vec::with_capacity(parts);
    let mut start = 0usize;
    for p in 0..parts {
        let end = if p + 1 == parts {
            n
        } else {
            start + n / parts + usize::from(p < n % parts)
        };
        out.push((start, end));
        start = end;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet;

    fn entries(s: &str) -> Vec<u8> {
        alphabet::encode(s.as_bytes())
    }

    #[test]
    fn identical_sequences_share_seeds() {
        let shape = Shape::contiguous(4);
        let red = Reduction::murphy10();
        // `--default` uses 0.8, scaled to the shape weight like setup.cpp does.
        let cut = 0.8 * std::f64::consts::LN_2 * shape.weight as f64;
        // Gly, Cys, Pro and Asp keep every 4-mer at four distinct murphy10
        // buckets, so the complexity filter keeps all seven seeds.
        let a = entries("MKWVTFGCPD");
        let q = build(&[&a], &shape, &red, 0, cut);
        let r = build(&[&a], &shape, &red, 7, cut);
        assert_eq!(q.len(), 7);
        let mut hits = vec![HashMap::new(); 1];
        join_into(&q, &r, &mut hits, 4096);
        let mut h = hits[0]
            .get(&7)
            .expect("hit on the reference sequence")
            .clone();
        assert_eq!(h.len(), 7);
        h.sort_by_key(|x| x.qpos);
        assert_eq!((h[0].qpos, h[0].spos), (0, 0));
    }

    #[test]
    fn low_complexity_seeds_are_filtered() {
        let shape = Shape::contiguous(8);
        let red = Reduction::murphy10();
        let cut = 0.8 * std::f64::consts::LN_2 * shape.weight as f64;
        assert!(!seed_is_complex(&entries("AAAAAAAAAA"), &shape, &red, cut));
        assert!(seed_is_complex(&entries("MKWVTFGCPD"), &shape, &red, cut));
        // Six Ala plus two Cys score ln(28) = 3.33: the scaled cut rejects the
        // seed while the unscaled `--seed-cut` value would accept it.
        assert!(!seed_is_complex(&entries("AAAAAACC"), &shape, &red, cut));
        assert!(seed_is_complex(&entries("AAAAAACC"), &shape, &red, 0.8));
    }

    #[test]
    fn chunking_covers_everything() {
        assert_eq!(chunks(10, 4).iter().map(|c| c.1 - c.0).sum::<usize>(), 10);
        assert_eq!(chunks(10, 4), vec![(0, 3), (3, 6), (6, 8), (8, 10)]);
        assert_eq!(chunks(0, 4), vec![]);
    }
}
