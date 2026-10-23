//! Gumbel statistics for gapped local alignments.
//!
//! DIAMOND does not use the textbook `K * m * n * exp(-lambda S)` formula.
//! It delegates to the ALP library (Sheetlin, Park, Frith & Spouge), which
//! replaces the search space `m * n` by an *effective* area that shrinks as
//! the score grows. That finite-size correction is what makes DIAMOND's
//! E-values differ from BLAST's by many orders of magnitude for strong hits,
//! so the formula is reproduced here bit-for-bit.
//!
//! The ALP parameters (`alpha`, `alpha_v`, `sigma`) are read from DIAMOND's
//! `src/stats/matrices/*.h` by `scripts/gen_matrices.py`; the intercepts are
//! then derived from the gap costs exactly like `alp_params()` in
//! `stats/score_matrix.cpp`.

// The coefficients below are transcribed verbatim from GSL and ALP and carry
// more digits than an f64 can represent. Truncating them would silently change
// the low bits of every E-value, which is exactly what this module reproduces.
#![allow(clippy::excessive_precision)]

use crate::matrix::LN_2;

/// `nat_cut_off_in_max` from `sls_pvalues.cpp`: finite-size correction
/// cut-off, expressed in nats.
const NAT_CUT_OFF_IN_MAX: f64 = 2.0;

/// `GSL_ROOT6_DBL_EPSILON`.
const ROOT6_EPS: f64 = 2.460_783_300_575_925_1e-3;

const SQRT_PI: f64 = 1.772_453_850_905_516;
const INV_SQRT_2PI: f64 = 0.398_942_280_401_432_7;

/// ALP Gumbel parameters plus the derived finite-size correction thresholds.
#[derive(Clone, Debug)]
pub struct AlpParams {
    pub lambda: f64,
    pub k: f64,
    /// `a_I`/`b_I`: length intercept and slope for the subject.
    a: f64,
    b: f64,
    alpha: f64,
    beta: f64,
    /// `a_J`/`b_J`: the same for the query. DIAMOND passes identical values.
    a_j: f64,
    b_j: f64,
    alpha_j: f64,
    beta_j: f64,
    sigma: f64,
    tau: f64,
    vi_thr: f64,
    vj_thr: f64,
    c_thr: f64,
}

impl AlpParams {
    /// Build the parameters the way `alp_params()` in `score_matrix.cpp` does.
    ///
    /// `alpha`/`alpha_v`/`sigma` come from the gapped row of the matrix table,
    /// `u_alpha`/`u_alpha_v` from its ungapped row.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        lambda: f64,
        k: f64,
        alpha: f64,
        alpha_v: f64,
        sigma: f64,
        gap_open: i32,
        gap_extend: i32,
        u_alpha: f64,
        u_alpha_v: f64,
    ) -> AlpParams {
        let g = (gap_open + gap_extend) as f64;
        let b = 2.0 * g * (u_alpha - alpha);
        let beta = 2.0 * g * (u_alpha_v - alpha_v);
        let tau = 2.0 * g * (u_alpha_v - sigma);
        AlpParams {
            lambda,
            k,
            a: alpha,
            b,
            alpha: alpha_v,
            beta,
            a_j: alpha,
            b_j: b,
            alpha_j: alpha_v,
            beta_j: beta,
            sigma,
            tau,
            vi_thr: (NAT_CUT_OFF_IN_MAX * alpha_v / lambda).max(0.0),
            vj_thr: (NAT_CUT_OFF_IN_MAX * alpha_v / lambda).max(0.0),
            c_thr: (NAT_CUT_OFF_IN_MAX * sigma / lambda).max(0.0),
        }
    }

    /// Textbook Karlin-Altschul E-value; used when the ALP area degenerates.
    pub fn blast_evalue(&self, raw: i32, query_len: usize, db_letters: u64) -> f64 {
        let exponent = -self.lambda * raw as f64 + self.k.ln();
        if exponent > 700.0 {
            f64::INFINITY
        } else {
            query_len as f64 * db_letters as f64 * exponent.exp()
        }
    }

    /// DIAMOND's E-value: `area * K * exp(-lambda * S) * db_letters / subject_len`.
    ///
    /// `m` is the subject length, `n` the query length, matching the argument
    /// order of `ScoreMatrix::evalue()` in the C++ code.
    pub fn evalue(&self, raw: i32, query_len: usize, subject_len: usize, db_letters: u64) -> f64 {
        let y = raw as f64;
        let m = subject_len.max(1) as f64;
        let n = query_len.max(1) as f64;
        let log_area = self.log_area(y, m, n);
        if !log_area.is_finite() {
            return self.blast_evalue(raw, query_len, db_letters);
        }
        let log_e = log_area + self.k.ln() - self.lambda * y;
        if log_e > 700.0 {
            f64::INFINITY
        } else {
            log_e.exp() * db_letters as f64 / m
        }
    }

    /// `pvalues::log_area()` from `sls_pvalues.cpp`.
    pub fn log_area(&self, y: f64, m: f64, n: f64) -> f64 {
        let m_li_y = m - (self.a * y + self.b);
        let vi_y = self.vi_thr.max(self.alpha * y + self.beta);
        let sqrt_vi = vi_y.sqrt();
        let m_f = if sqrt_vi == 0.0 {
            1e100
        } else {
            m_li_y / sqrt_vi
        };
        let log_p_m = 0.5f64.ln() + log_erfc(-0.5f64.sqrt() * m_f);
        let log_e_m = INV_SQRT_2PI.ln() - 0.5 * m_f * m_f;
        let log_sqrt_e_m = ln_or_neg_inf(sqrt_vi) + log_e_m;
        let log_p1 = if m_li_y < 0.0 {
            log_diff(log_sqrt_e_m, ln_or_neg_inf(-m_li_y) + log_p_m)
        } else {
            log_sum(log_sqrt_e_m, ln_or_neg_inf(m_li_y) + log_p_m)
        };

        let n_lj_y = n - (self.a_j * y + self.b_j);
        let vj_y = self.vj_thr.max(self.alpha_j * y + self.beta_j);
        let sqrt_vj = vj_y.sqrt();
        let n_f = if sqrt_vj == 0.0 {
            1e100
        } else {
            n_lj_y / sqrt_vj
        };
        let log_p_n = 0.5f64.ln() + log_erfc(-0.5f64.sqrt() * n_f);
        let log_e_n = INV_SQRT_2PI.ln() - 0.5 * n_f * n_f;
        let log_sqrt_e_n = ln_or_neg_inf(sqrt_vj) + log_e_n;
        let log_p2 = if n_lj_y < 0.0 {
            log_diff(log_sqrt_e_n, ln_or_neg_inf(-n_lj_y) + log_p_n)
        } else {
            log_sum(log_sqrt_e_n, ln_or_neg_inf(n_lj_y) + log_p_n)
        };

        let c_y = self.c_thr.max(self.sigma * y + self.tau);
        let log_c_y = ln_or_neg_inf(c_y);
        log_sum(log_p1 + log_p2, log_c_y + log_p_m + log_p_n)
    }
}

#[inline]
fn ln_or_neg_inf(x: f64) -> f64 {
    if x > 0.0 {
        x.ln()
    } else {
        f64::NEG_INFINITY
    }
}

fn log_sum(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if !hi.is_finite() {
        return hi;
    }
    if lo == f64::NEG_INFINITY {
        return hi;
    }
    hi + (1.0 + (lo - hi).exp()).ln()
}

fn log_diff(a: f64, b: f64) -> f64 {
    if a < b {
        return f64::NEG_INFINITY;
    }
    if !a.is_finite() {
        return a;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let t = 1.0 - (b - a).exp();
    if t <= 0.0 {
        f64::NEG_INFINITY
    } else {
        a + t.ln()
    }
}

// ---------------------------------------------------------------------------
// log(erfc(x)) — port of GSL's `gsl_sf_log_erfc_e` as bundled with DIAMOND.
// ---------------------------------------------------------------------------

/// `log(erfc(x))` for the whole real line.
pub fn log_erfc(x: f64) -> f64 {
    if x < -8.0 {
        // erfc(-8) differs from 2 by less than 1e-29, i.e. it rounds to 2.
        return LN_2;
    }
    if x * x < 10.0 * ROOT6_EPS {
        let y = x / SQRT_PI;
        const C3: f64 = (4.0 - std::f64::consts::PI) / 3.0;
        const C4: f64 = 2.0 * (1.0 - std::f64::consts::PI / 3.0);
        const C5: f64 = -0.001_829_764_677_455_021;
        const C6: f64 = 0.026_296_515_210_574_65;
        const C7: f64 = -0.016_215_753_788_354_04;
        const C8: f64 = 0.001_259_939_617_621_16;
        const C9: f64 = 0.005_569_646_491_38;
        const C10: f64 = -0.004_556_333_980_2;
        const C11: f64 = 0.000_946_158_903_2;
        const C12: f64 = 0.001_320_024_317_4;
        const C13: f64 = -0.001_429_06;
        const C14: f64 = 0.000_482_04;
        let mut series = C8 + y * (C9 + y * (C10 + y * (C11 + y * (C12 + y * (C13 + C14 * y)))));
        series = y
            * (1.0 + y * (1.0 + y * (C3 + y * (C4 + y * (C5 + y * (C6 + y * (C7 + y * series)))))));
        return -2.0 * series;
    }
    if x > 8.0 {
        return erfc8_sum(x).ln() - x * x;
    }
    erfc(x).ln()
}

fn erfc(x: f64) -> f64 {
    let ax = x.abs();
    let value = if ax <= 1.0 {
        cheb_eval(&ERFC_XLT1, 2.0 * ax - 1.0)
    } else if ax <= 5.0 {
        (-x * x).exp() * cheb_eval(&ERFC_X15, 0.5 * (ax - 3.0))
    } else if ax < 10.0 {
        (-x * x).exp() / ax * cheb_eval(&ERFC_X510, (2.0 * ax - 15.0) / 5.0)
    } else {
        erfc8_sum(ax) * (-x * x).exp()
    };
    if x < 0.0 {
        2.0 - value
    } else {
        value
    }
}

/// Chebyshev evaluation on [-1, 1] (GSL's `cheb_eval_e` without the error).
fn cheb_eval(coeffs: &[f64], x: f64) -> f64 {
    let n = coeffs.len() - 1;
    let y2 = 2.0 * x;
    let mut d = 0.0;
    let mut dd = 0.0;
    for j in (1..=n).rev() {
        let temp = d;
        d = y2 * d - dd + coeffs[j];
        dd = temp;
    }
    x * d - dd + 0.5 * coeffs[0]
}

/// `erfc8_sum()`: rational approximation of `erfc(x) * exp(x^2)` for 8 < x < 100.
fn erfc8_sum(x: f64) -> f64 {
    const P: [f64; 6] = [
        2.978_865_626_393_992_888_62,
        7.409_740_605_964_741_794_425,
        6.160_209_853_109_630_544_090_6,
        5.019_049_726_784_267_463_450_058,
        1.275_366_644_729_965_952_479_585_264,
        0.564_189_583_547_755_074_125_320_170_4,
    ];
    const Q: [f64; 7] = [
        3.369_075_206_982_752_767_7,
        9.608_965_327_192_787_870_698,
        17.081_440_747_466_004_315_710_95,
        12.048_951_927_855_129_036_034_049_1,
        9.396_034_016_235_054_150_430_579_648,
        2.260_528_520_767_326_969_591_866_945,
        1.0,
    ];
    let mut num = P[5];
    for &p in P.iter().rev().skip(1) {
        num = x * num + p;
    }
    let mut den = Q[6];
    for &q in Q.iter().rev().skip(1) {
        den = x * den + q;
    }
    num / den
}

/// Chebyshev fit for `erfc((t+1)/2)`, -1 < t < 1.
const ERFC_XLT1: [f64; 20] = [
    1.060_734_164_217_699_803_451_741_550_56,
    -0.425_824_458_043_810_435_692_047_352_91,
    0.049_552_626_796_204_340_403_576_830_80,
    0.004_492_934_887_683_827_495_580_012_42,
    -0.001_291_941_046_584_969_534_942_247_61,
    -0.000_018_363_892_921_493_962_704_169_79,
    0.000_022_111_147_040_995_262_915_385_56,
    -5.233_374_852_342_571_346_736_931_790_20e-7,
    -2.781_847_888_335_378_853_825_309_895_78e-7,
    1.411_580_927_488_131_145_603_166_842_49e-8,
    2.725_712_963_305_616_999_845_391_418_65e-9,
    -2.063_439_048_720_706_294_064_014_924_76e-10,
    -2.142_739_919_967_853_679_242_014_018_12e-11,
    2.229_902_555_393_582_045_802_850_981_19e-12,
    1.362_500_746_506_982_805_758_079_341_55e-13,
    -1.951_440_109_222_930_918_989_959_130_38e-14,
    -6.856_271_692_317_045_994_428_063_706_90e-16,
    1.445_064_928_696_999_382_395_216_074_93e-16,
    2.459_353_064_605_364_880_375_762_000_30e-18,
    -9.295_995_612_205_233_960_073_593_285_40e-19,
];

/// Chebyshev fit for `erfc(x) exp(x^2)`, 1 < x < 5, x = 2t + 3.
const ERFC_X15: [f64; 25] = [
    0.440_458_320_243_381_110_776_374_666_16,
    -0.143_958_836_762_168_335_790_826_895_326,
    0.044_786_499_817_939_267_247_056_666_937,
    -0.013_343_124_200_271_211_203_618_353_102,
    0.003_824_682_739_750_469_767_692_372_556,
    -0.001_058_699_227_195_126_547_306_482_530,
    0.000_283_859_419_210_073_742_736_310_108,
    -0.000_073_906_170_662_206_760_483_959_432,
    0.000_018_725_312_521_489_179_015_872_934,
    -4.625_309_811_649_194_451_312_972_644_30e-6,
    1.115_586_572_444_328_574_878_840_064_22e-6,
    -2.630_986_626_508_341_300_678_088_327_25e-7,
    6.074_621_227_245_517_773_721_194_087_10e-8,
    -1.374_608_655_398_654_447_772_510_117_93e-8,
    3.051_570_519_054_751_455_200_967_172_10e-9,
    -6.651_747_897_203_107_137_573_077_247_90e-10,
    1.424_833_462_732_077_844_897_929_997_06e-10,
    -3.001_411_273_953_239_020_920_187_445_45e-11,
    6.221_717_926_453_480_914_729_140_012_50e-12,
    -1.269_946_392_256_684_968_761_528_365_55e-12,
    2.553_858_830_332_575_754_026_818_453_85e-13,
    -5.062_582_375_070_386_983_922_654_997_70e-14,
    9.897_054_094_783_273_216_412_642_271_10e-15,
    -1.906_859_787_891_921_810_519_610_249_95e-15,
    3.508_266_480_327_378_492_451_137_573_40e-16,
];

/// Chebyshev fit for `erfc(x) x exp(x^2)`, 5 < x < 10, x = (5t + 15)/2.
const ERFC_X510: [f64; 20] = [
    1.116_849_901_235_456_986_842_978_658_08,
    0.003_736_240_359_381_998_520_654_927_536,
    -0.000_916_623_948_045_470_238_763_619_870,
    0.000_199_094_325_044_940_833_965_078_819,
    -0.000_040_276_384_918_650_072_591_781_859,
    7.765_152_646_970_610_494_771_276_057_90e-6,
    -1.444_647_942_066_890_704_020_992_253_01e-6,
    2.613_119_303_434_639_583_934_852_419_47e-7,
    -4.618_330_266_348_441_523_453_040_955_60e-8,
    8.002_531_115_129_436_015_987_321_443_40e-9,
    -1.362_911_148_627_930_313_957_121_220_89e-9,
    2.285_704_830_901_608_696_076_830_877_22e-10,
    -3.780_225_215_632_518_050_440_569_745_60e-11,
    6.172_536_838_745_282_857_299_104_621_30e-12,
    -9.960_192_909_553_168_884_458_305_974_30e-13,
    1.589_531_437_069_807_702_695_067_260_00e-13,
    -2.510_459_710_471_625_099_995_274_283_16e-14,
    3.926_078_289_891_258_100_135_812_875_60e-15,
    -6.079_706_193_841_603_743_925_354_534_20e-16,
    9.126_006_072_647_947_173_155_074_776_70e-17,
];

#[cfg(test)]
mod tests {
    use super::*;

    fn blosum62() -> AlpParams {
        // BLOSUM62, gap costs 11/1, from src/stats/matrices/blosum62.h.
        AlpParams::new(0.267, 0.041, 1.9, 42.6028, 43.6362, 11, 1, 0.7916, 4.964660)
    }

    #[test]
    fn log_erfc_matches_reference_values() {
        // Reference values computed with a correctly-rounded erfc.
        assert_eq!(log_erfc(0.0), 0.0);
        assert!((log_erfc(1.0) - (-1.849_605_509_933_248_5)).abs() < 1e-12);
        assert!((log_erfc(4.0) - (-17.987_778_312_103_007)).abs() < 1e-11);
        assert!((log_erfc(5.0) - (-27.200_889_545_537_436)).abs() < 1e-11);
        assert!((log_erfc(-5.0) - 0.693_147_180_559_176_6).abs() < 1e-12);
        assert!((log_erfc(-20.0) - LN_2).abs() < 1e-15);
    }

    #[test]
    fn alp_evalues_match_diamond() {
        // Three HSPs from a DIAMOND 2.2.8 blastp run against 7860 sequences:
        // (raw score, subject length, query length, db letters, E-value).
        let p = blosum62();
        let db = 2_263_996u64;
        let q = 497usize;
        let cases = [
            (818, 454, 3.21e-105),
            (262, 456, 3.89e-25),
            (198, 601, 7.60e-17),
        ];
        for (raw, slen, expected) in cases {
            let got = p.evalue(raw, q, slen, db);
            let ratio = got / expected;
            assert!(
                (ratio - 1.0).abs() < 0.01,
                "raw={raw}: got {got:.3e}, expected {expected:.3e} (ratio {ratio})"
            );
        }
    }

    #[test]
    fn weak_hits_stay_blast_like() {
        // For a weak alignment the finite-size correction is close to m*n.
        let p = blosum62();
        let e = p.evalue(60, 300, 300, 1_000_000);
        let blast = p.blast_evalue(60, 300, 1_000_000);
        assert!(e > 0.0 && e.is_finite());
        assert!((e / blast).ln().abs() < 3.0, "ratio {}", e / blast);
    }
}
