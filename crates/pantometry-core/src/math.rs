//! `exp` and `ln` of the workspace's own, built from `+ − × ÷` alone, so that they give the same
//! bits on every platform.
//!
//! # Why the kernel has these
//!
//! The determinism convention says a result is bit-for-bit the same across platforms,
//! optimisation levels, WebAssembly and thread counts. IEEE-754 makes that true of `+ − × ÷` and
//! `sqrt`, which it requires to be correctly rounded. It does not make it true of `exp` or `ln`:
//! `f64::exp` and `f64::ln` call the platform's C library, which is not required to round
//! correctly and is not the same code on two operating systems. A crate that calls them is
//! reproducible on one machine and not across machines, and a domain crate that calls them has
//! to say so of itself.
//!
//! [`exp`] and [`ln`] here are a fixed sequence of IEEE operations — integer manipulation of the
//! bit pattern, two tables of constants, and `+ − × ÷` in a written order, with no fused
//! multiply–add and no call out — so for every non-NaN argument they are the same function on
//! every target Rust compiles to. A NaN comes back NaN, but Rust does not fix its payload bits.
//! Rust does not contract `a * b + c` into an FMA, and its `f64` arithmetic on every tier-1 target
//! and on `wasm32` is IEEE binary64 rounding to nearest, which is all the evaluation below
//! assumes. `the_bits_are_pinned` holds a digest of both at 8192 points each, and CI runs it on
//! Linux, macOS on arm64, Windows, release with LTO, and `wasm32-wasip1`.
//!
//! They are in the kernel because they know nothing about any physics: they are the maths no
//! domain owns, as [`crate::vector`] and [`crate::transform`] are, and a second domain that wanted
//! them should not have to depend on the first. Nothing in the kernel calls them yet.
//!
//! They are also faster than the platform's on `x86_64-pc-windows-gnu`, which is why they were
//! written: in release, over the arguments of a pair-sum hot path, 1.96–2.03 ns a call for [`exp`]
//! against `f64::exp`'s 33.0–33.7, and 2.86–2.92 ns for [`ln`] against `f64::ln`'s 17.4–17.6
//! (`tests/exp_and_ln_against_mpmath.rs`, `the_cost`).
//!
//! # Where they come from
//!
//! **The code is original to this workspace**, written from published mathematics rather than
//! from any library's source, and is under the crate's own `MIT OR Apache-2.0`. What it follows:
//!
//! - **The table-driven method**, P. T. P. Tang, "Table-driven implementation of the exponential
//!   function in IEEE floating-point arithmetic", *ACM Trans. Math. Softw.* 15(2), 144–157
//!   (1989), and "Table-driven implementation of the logarithm function in IEEE floating-point
//!   arithmetic", *ACM Trans. Math. Softw.* 16(4), 378–400 (1990): reduce the argument against a
//!   table of exactly known values so that what is left is small enough for a short series.
//! - **Cody–Waite reduction**, W. J. Cody and W. Waite, *Software Manual for the Elementary
//!   Functions* (Prentice-Hall, 1980): `ln 2` split into a head with trailing zeros, so that
//!   multiples of it are exact, and a tail. J.-M. Muller, *Elementary Functions: Algorithms and
//!   Implementation* (Birkhäuser), covers both, and rounding to an integer by adding `1.5 · 2⁵²`.
//! - **Exact transformations**, T. J. Dekker, "A floating-point technique for extending the
//!   available precision", *Numer. Math.* 18, 224–242 (1971) — the error of a sum recovered
//!   exactly when one addend dominates — and Sterbenz's lemma, that `a − b` is exact when
//!   `b/2 ≤ a ≤ 2b`.
//! - **The series are Taylor's**: `eʳ − 1 = r + r²/2! + … + r⁵/5!` and
//!   `ln(1 + r) − r = −r²/2 + r³/3 − … + r⁷/7`, their coefficients `1/n!` and `±1/n` as Rust
//!   writes them. No coefficient is fitted.
//!
//! The constants — a table of `2^(j/128)` as a double-double, a table of `c ≈ 1/F` and `−ln c`,
//! and the split of `ln 2` — are computed by `tools/math-reference/generate.py` with mpmath at 50
//! digits and are what its `--exp-table`, `--ln-table` and `--constants` print, word for word.
//! The tests here check each against a closed form computed in double-double arithmetic in the
//! test itself, not against the script.
//!
//! **[`exp`]**: `x = (128k + j)·(ln 2/128) + r`, with `k + j/128` the nearest multiple of
//! `1/128` to `x/ln 2` and `|r| ≤ ln 2/256`; then `eˣ = 2ᵏ · 2^(j/128) · eʳ`, the middle factor
//! from the table and the last from the series. **[`ln`]**: `x = 2ᵏ m` with `m` in
//! `[0.703125, 1.40625)`, so that `k = 0` all about 1; `F` is `m` rounded to eight bits after its
//! leading one, `c` a nine-bit approximation of `1/F` from the table, and
//! `ln x = k ln 2 − ln c + ln(1 + r)` with `r = m c − 1`. Written as `(F c − 1) + (m − F) c`,
//! every operation in `r` is exact — the table's cells are narrow enough and `c` short enough for
//! that, which its test checks cell by cell — so nothing is divided and nothing is lost. About 1,
//! `c = 1` and `r = m − 1`. What each step costs in accuracy is derived where it is done, and the
//! bounds assembled from them are written beside the tests that hold them.
//!
//! # How accurate, measured rather than quoted
//!
//! The claim is **less than one ulp**, and it is stated tighter because it can be derived:
//! **0.520 ulp for [`exp`]** everywhere; for [`ln`], **0.508** in the binade of 1 (`k = 0`, `x`
//! in `[0.703125, 1.40625)`) and **0.501** everywhere else. The
//! derivations are written beside the bounds in `tests/exp_and_ln_against_mpmath.rs`, which holds
//! both functions to them against [mpmath](https://mpmath.org/) at 50 significant digits on the
//! points `tools/math-reference/generate.py` draws: the whole domain of each, every reduction
//! boundary `(k + ½) ln 2` of `exp`, its table's boundaries `(n + ½) ln 2/128` for `|n| ≤ 512`,
//! the thresholds of its branches, every binade's boundary of `ln` at `√2` and at `0.703125`,
//! every edge of its 256 cells at `k = −1, 0, 1`, every power of two, every subnormal, and densely
//! about 1. That is 64 993 and 76 308 committed points; the same recipe at twenty times the random
//! points, 804 287 and 737 812, was run once outside the repository:
//!
//! | | bound | worst, committed | worst, ×20 | where (×20) | not correctly rounded (×20) |
//! | --- | --- | --- | --- | --- | --- |
//! | [`exp`] | 0.520 | 0.5106 ulp | 0.5106 ulp | `x = 0.55506`, at a boundary of `j` | 0.074% |
//! | [`ln`], `c = 1` (`x` in `[1 − 3·2⁻¹⁰, 1 + 2⁻⁹)`) | 0.508 | 0.5000 ulp | 0.5010 ulp | `x = 1.0017` | |
//! | [`ln`], `k = 0`, `c ≠ 1` | 0.508 | 0.5000 ulp | 0.5002 ulp | `x = 0.97946` | |
//! | [`ln`], `k ≠ 0` | 0.501 | 0.5000 ulp | 0.5000 ulp | | |
//! | [`ln`], all | | | | | 0.001% (6 points) |
//!
//! An ulp here is the spacing of doubles at the exact value; below the smallest normal it is the
//! subnormal spacing, `2⁻¹⁰⁷⁴`. Neither function is correctly rounded: for the fraction shown the
//! farther of the two neighbouring doubles is returned.
//!
//! # Special values
//!
//! As IEEE-754 and `f64::exp` / `f64::ln` have them, each tested by name:
//!
//! | `x` | `exp(x)` | `ln(x)` |
//! | --- | --- | --- |
//! | NaN | NaN | NaN |
//! | `+∞` | `+∞` | `+∞` |
//! | `−∞` | `+0` | NaN |
//! | `±0` | `1`, exactly | `−∞` |
//! | `1` | | `+0`, exactly |
//! | `< 0` | | NaN |
//! | `> 709.782712893383973096` | `+∞` (overflow) | |
//! | `< −745.13321910194110842` | `+0` (underflow) | |
//! | subnormal | `1` | finite, by scaling into the normal range first |
//!
//! Between `−745.13…` and `−708.39…` the result of [`exp`] is subnormal and is rounded once onto
//! the subnormal grid.

// The constants below are bit patterns, so that each is the double it says without a decimal
// conversion in between; `generate.py --constants` prints them.
//
// `LN2_HI` is ln 2 rounded to a multiple of 2^-35 — 35 significant bits — and `LN2_LO` the double
// nearest ln 2 − LN2_HI. `INV_STEP` is the double nearest 128/ln 2.
const LN2_HI: u64 = 0x3fe6_2e42_fefc_0000;
const LN2_LO: u64 = 0xbdac_610c_a86c_3899;
const INV_STEP: u64 = 0x4067_1547_652b_82fe;

/// `1.5 · 2⁵²`. Added to a value under `2⁵¹` in magnitude, the sum lies in `[2⁵², 2⁵³)`, where
/// the doubles are the integers, so the addition rounds to the nearest integer, ties to even.
const SHIFT: f64 = 6_755_399_441_055_744.0;

/// `|x|` at and above which [`exp`] leaves its main path: `708.0`. Below it `k` is within
/// `[−1022, 1021]` and `2^(j/128) eʳ ≥ 2^(−0.43)` when `k = −1022`, so the result is normal and
/// scaling by `2ᵏ` is exact.
const EXP_MAIN_LIMIT: u64 = 0x4086_2000_0000_0000;

/// Above this [`exp`] is `+∞` without computing; between the true threshold and this the scaling
/// overflows by itself. `709.8` is past `ln(MAX + ½ ulp) = 709.78271289338399…`.
const EXP_OVERFLOW_BEYOND: f64 = 709.8;

/// Below this [`exp`] is `+0` without computing: past `ln 2⁻¹⁰⁷⁵ = −745.13321910194121…`, where
/// half the smallest subnormal is.
const EXP_UNDERFLOW_BEYOND: f64 = -745.2;

/// `2⁻¹⁰²²`, the smallest normal double.
const TWO_M1022: f64 = f64::MIN_POSITIVE;

/// The series' coefficients, `1/n!`.
const E2: f64 = 0.5;
const E3: f64 = 1.0 / 6.0;
const E4: f64 = 1.0 / 24.0;
const E5: f64 = 1.0 / 120.0;

/// Per `j` in `0..128`, the bits of `T_hi` and `T_lo`: `2^(j/128) = T_hi + T_lo` to `2⁻¹⁰⁶`
/// relatively, `T_hi` the nearest double and `T_lo` the nearest double to the rest.
#[rustfmt::skip]
const EXP_TABLE: [u64; 256] = [
    0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
    0x3ff0_163d_a9fb_3335, 0x3c9b_6129_9ab8_cdb7,
    0x3ff0_2c9a_3e77_8061, 0xbc71_9083_535b_085d,
    0x3ff0_4315_e86e_7f85, 0xbc90_a31c_1977_c96e,
    0x3ff0_59b0_d315_8574, 0x3c8d_73e2_a475_b465,
    0x3ff0_706b_29dd_f6de, 0xbc8c_91df_e2b1_3c27,
    0x3ff0_8745_1875_9bc8, 0x3c61_86be_4bb2_84ff,
    0x3ff0_9e3e_cac6_f383, 0x3c91_4878_1831_6136,
    0x3ff0_b558_6cf9_890f, 0x3c98_a62e_4adc_610b,
    0x3ff0_cc92_2b72_47f7, 0x3c90_1edc_16e2_4f71,
    0x3ff0_e3ec_32d3_d1a2, 0x3c40_3a17_27c5_7b53,
    0x3ff0_fb66_affe_d31b, 0xbc6b_9bed_c44e_bd7b,
    0x3ff1_1301_d012_5b51, 0xbc96_c510_3944_9b3a,
    0x3ff1_2abd_c06c_31cc, 0xbc51_b514_b36c_a5c7,
    0x3ff1_429a_aea9_2de0, 0xbc93_2fbf_9af1_369e,
    0x3ff1_5a98_c8a5_8e51, 0x3c82_406a_b9ee_ab0a,
    0x3ff1_72b8_3c7d_517b, 0xbc81_9041_b9d7_8a76,
    0x3ff1_8af9_388c_8dea, 0xbc91_1023_d197_0f6c,
    0x3ff1_a35b_eb6f_cb75, 0x3c8e_5b4c_7b49_68e4,
    0x3ff1_bbe0_8404_5cd4, 0xbc99_5386_352e_f607,
    0x3ff1_d487_3168_b9aa, 0x3c9e_016e_00a2_643c,
    0x3ff1_ed50_22fc_d91d, 0xbc91_df98_027b_b78c,
    0x3ff2_063b_8862_8cd6, 0x3c8d_c775_814a_8495,
    0x3ff2_1f49_917d_dc96, 0x3c82_a97e_9494_a5ee,
    0x3ff2_387a_6e75_6238, 0x3c99_b07e_b6c7_0573,
    0x3ff2_51ce_4fb2_a63f, 0x3c8a_c155_bef4_f4a4,
    0x3ff2_6b45_65e2_7cdd, 0x3c82_bd33_9940_e9d9,
    0x3ff2_84df_e1f5_6381, 0xbc9a_4c3a_8c3f_0d7e,
    0x3ff2_9e9d_f51f_dee1, 0x3c86_12e8_afad_1255,
    0x3ff2_b87f_d0da_d990, 0xbc41_0adc_d638_1aa4,
    0x3ff2_d285_a6e4_030b, 0x3c90_0247_54db_41d5,
    0x3ff2_ecaf_a93e_2f56, 0x3c71_ca0f_45d5_2383,
    0x3ff3_06fe_0a31_b715, 0x3c86_f46a_d231_82e4,
    0x3ff3_2170_fc4c_d831, 0x3c8a_9ce7_8e18_047c,
    0x3ff3_3c08_b264_16ff, 0x3c93_2721_8436_59a6,
    0x3ff3_56c5_5f92_9ff1, 0xbc8b_5cee_5c4e_4628,
    0x3ff3_71a7_373a_a9cb, 0xbc96_3aea_bf42_eae2,
    0x3ff3_8cae_6d05_d866, 0xbc9e_958d_3c99_04bd,
    0x3ff3_a7db_34e5_9ff7, 0xbc75_e436_d661_f5e3,
    0x3ff3_c32d_c313_a8e5, 0xbc9e_fff8_375d_29c3,
    0x3ff3_dea6_4c12_3422, 0x3c8a_da09_11f0_9ebc,
    0x3ff3_fa45_04ac_801c, 0xbc97_d023_f956_f9f3,
    0x3ff4_160a_21f7_2e2a, 0xbc5e_f369_1c30_9278,
    0x3ff4_31f5_d950_a897, 0xbc81_c7dd_e35f_7999,
    0x3ff4_4e08_6061_892d, 0x3c48_9b7a_04ef_80d0,
    0x3ff4_6a41_ed1d_0057, 0x3c9c_944b_d164_8a76,
    0x3ff4_86a2_b5c1_3cd0, 0x3c73_c1a3_b690_62f0,
    0x3ff4_a32a_f0d7_d3de, 0x3c99_cb62_f3d1_be56,
    0x3ff4_bfda_d536_2a27, 0x3c7d_4397_afec_42e2,
    0x3ff4_dcb2_99fd_dd0d, 0x3c98_ecdb_bc6a_7833,
    0x3ff4_f9b2_769d_2ca7, 0xbc94_b309_d259_57e3,
    0x3ff5_16da_a2cf_6642, 0xbc8f_7685_69bd_93ef,
    0x3ff5_342b_569d_4f82, 0xbc80_7abe_1db1_3cad,
    0x3ff5_51a4_ca5d_920f, 0xbc8d_689c_efed_e59b,
    0x3ff5_6f47_36b5_27da, 0x3c99_bb2c_011d_93ad,
    0x3ff5_8d12_d497_c7fd, 0x3c82_95e1_5b9a_1de8,
    0x3ff5_ab07_dd48_5429, 0x3c96_324c_0546_47ad,
    0x3ff5_c926_8a59_46b7, 0x3c3c_4b1b_8169_86a2,
    0x3ff5_e76f_15ad_2148, 0x3c9b_a6f9_3080_e65e,
    0x3ff6_05e1_b976_dc09, 0xbc93_e242_9b56_de47,
    0x3ff6_247e_b03a_5585, 0xbc93_83c1_7e40_b497,
    0x3ff6_4346_34cc_c320, 0xbc8c_483c_759d_8933,
    0x3ff6_6238_8255_2225, 0xbc9b_b609_8759_1c34,
    0x3ff6_8155_d44c_a973, 0x3c60_38ae_44f7_3e65,
    0x3ff6_a09e_667f_3bcd, 0xbc9b_dd34_13b2_6456,
    0x3ff6_c012_750b_dabf, 0xbc72_8956_67ff_0b0d,
    0x3ff6_dfb2_3c65_1a2f, 0xbc6b_be3a_683c_88ab,
    0x3ff6_ff7d_f951_9484, 0xbc88_3c0f_2586_0ef6,
    0x3ff7_1f75_e8ec_5f74, 0xbc81_6e47_8688_7a99,
    0x3ff7_3f9a_48a5_8174, 0xbc90_a8d9_6c65_d53c,
    0x3ff7_5feb_5642_67c9, 0xbc90_2459_5731_6dd3,
    0x3ff7_8069_4fde_5d3f, 0x3c98_66b8_0a02_162d,
    0x3ff7_a114_73eb_0187, 0xbc84_1577_ee04_992f,
    0x3ff7_c1ed_0130_c132, 0x3c9f_124c_d116_4dd6,
    0x3ff7_e2f3_36cf_4e62, 0x3c70_5d02_ba15_797e,
    0x3ff8_0427_543e_1a12, 0xbc92_7c86_626d_972b,
    0x3ff8_2589_994c_ce13, 0xbc9d_4c1d_d415_32d8,
    0x3ff8_471a_4623_c7ad, 0xbc88_d684_a341_cdfb,
    0x3ff8_68d9_9b44_92ed, 0xbc9f_c6f8_9bd4_f6ba,
    0x3ff8_8ac7_d98a_6699, 0x3c99_94c2_f37c_b53a,
    0x3ff8_ace5_422a_a0db, 0x3c96_e9f1_5686_4b27,
    0x3ff8_cf32_16b5_448c, 0xbc70_d55e_32e9_e3aa,
    0x3ff8_f1ae_9915_7736, 0x3c85_cc13_a2e3_976c,
    0x3ff9_145b_0b91_ffc6, 0xbc9d_d679_2e58_2524,
    0x3ff9_3737_b0cd_c5e5, 0xbc67_5fc7_81b5_7ebc,
    0x3ff9_5a44_cbc8_520f, 0xbc76_4b7c_96a5_f039,
    0x3ff9_7d82_9fde_4e50, 0xbc9d_185b_7c1b_85d1,
    0x3ff9_a0f1_70ca_07ba, 0xbc91_73bd_91ce_e632,
    0x3ff9_c491_82a3_f090, 0x3c7c_7c46_b071_f2be,
    0x3ff9_e863_19e3_2323, 0x3c78_24ca_78e6_4c6e,
    0x3ffa_0c66_7b5d_e565, 0xbc93_5949_5d1c_d533,
    0x3ffa_309b_ec4a_2d33, 0x3c96_305c_7ddc_36ab,
    0x3ffa_5503_b23e_255d, 0xbc9d_2f6e_db8d_41e1,
    0x3ffa_799e_1330_b358, 0x3c9b_cb7e_cac5_63c7,
    0x3ffa_9e6b_5579_fdbf, 0x3c90_fac9_0ef7_fd31,
    0x3ffa_c36b_bfd3_f37a, 0xbc8f_9234_cae7_6cd0,
    0x3ffa_e89f_995a_d3ad, 0x3c97_a1cd_345d_cc81,
    0x3ffb_0e07_298d_b666, 0xbc9b_def5_4c80_e425,
    0x3ffb_33a2_b84f_15fb, 0xbc62_805e_3084_d708,
    0x3ffb_5972_8de5_593a, 0xbc9c_71df_bbba_6de3,
    0x3ffb_7f76_f2fb_5e47, 0xbc75_584f_7e54_ac3b,
    0x3ffb_a5b0_30a1_064a, 0xbc9e_fcd3_0e54_292e,
    0x3ffb_cc1e_904b_c1d2, 0x3c82_3dd0_7a2d_9e84,
    0x3ffb_f2c2_5bd7_1e09, 0xbc9e_fdca_3f6b_9c73,
    0x3ffc_199b_dd85_529c, 0x3c81_1065_8950_48dd,
    0x3ffc_40ab_5fff_d07a, 0x3c9b_4537_e083_c60a,
    0x3ffc_67f1_2e57_d14b, 0x3c92_884d_ff48_3cad,
    0x3ffc_8f6d_9406_e7b5, 0x3c71_acbc_4880_5c44,
    0x3ffc_b720_dcef_9069, 0x3c75_03cb_d1e9_49db,
    0x3ffc_df0b_555d_c3fa, 0xbc8d_d83b_5382_9d72,
    0x3ffd_072d_4a07_897c, 0xbc9c_bc37_4379_7a9c,
    0x3ffd_2f87_080d_89f2, 0xbc9d_487b_719d_8578,
    0x3ffd_5818_dcfb_a487, 0x3c82_ed02_d75b_3707,
    0x3ffd_80e3_16c9_8398, 0xbc91_1ec1_8bed_dfe8,
    0x3ffd_a9e6_03db_3285, 0x3c9c_2300_696d_b532,
    0x3ffd_d321_f301_b460, 0x3c92_da57_78f0_18c3,
    0x3ffd_fc97_337b_9b5f, 0xbc91_a5cd_4f18_4b5c,
    0x3ffe_2646_14f5_a129, 0xbc97_b627_817a_1496,
    0x3ffe_502e_e78b_3ff6, 0x3c83_9e89_80a9_cc8f,
    0x3ffe_7a51_fbc7_4c83, 0x3c92_d522_ca0c_8de2,
    0x3ffe_a4af_a2a4_90da, 0xbc9e_9c23_179c_2893,
    0x3ffe_cf48_2d8e_67f1, 0xbc9c_93f3_b411_ad8c,
    0x3ffe_fa1b_ee61_5a27, 0x3c9d_c7f4_86a4_b6b0,
    0x3fff_252b_376b_ba97, 0x3c93_a1a5_bf0d_8e43,
    0x3fff_5076_5b6e_4540, 0x3c99_d3e1_2dd8_a18b,
    0x3fff_7bfd_ad9c_be14, 0xbc9d_bb12_d006_350a,
    0x3fff_a7c1_819e_90d8, 0x3c87_4853_f3a5_931e,
    0x3fff_d3c2_2b8f_71f1, 0x3c62_eb74_9665_79e7,
];

/// `2ᵏ` for `k` in `[−1022, 1023]`, by its bits.
fn two_to(k: i64) -> f64 {
    f64::from_bits(((k + 1023) as u64) << 52)
}

/// The reduction and the series, shared by both paths of [`exp`]: `(k, T_hi, c)` with
/// `eˣ ≈ 2ᵏ (T_hi + c)` before the one rounding of the sum.
#[inline(always)]
fn exp_parts(x: f64) -> (i64, f64, f64) {
    let ln2_hi = f64::from_bits(LN2_HI);
    let ln2_lo = f64::from_bits(LN2_LO);
    // n, the nearest integer to x · 128/ln 2. |x · INV_STEP| < 2^17.1 here, so the sum is in
    // [2^52, 2^53), and its bits less SHIFT's are n: the doubles there are consecutive integers.
    let shifted = x * f64::from_bits(INV_STEP) + SHIFT;
    let n = shifted.to_bits() as i64 - SHIFT.to_bits() as i64;
    let nf = shifted - SHIFT;
    // r = x − n ln2/128, Cody–Waite. ln2_hi/128 has 35 significant bits and |n| < 2^18, so the
    // product is exact. The difference is exact too: x and the product are multiples of
    // 2^-61 or coarser wherever n ≠ 0 (|x| > 2^-9), and the difference is under 2^-8, so it has
    // at most 53 bits at that spacing. The tail's product and the last subtraction round:
    // |r − r*| ≤ 2^-62 + 2^-79.
    let r_hi = x - nf * (ln2_hi / 128.0);
    let r = r_hi - nf * (ln2_lo / 128.0);
    // eʳ − 1 to r⁵, by Estrin's scheme. |r| ≤ ln 2/256 · (1 + 2^-34): the truncation is at most
    // r⁶/720 < 2^-60.6, and the rounding at most 2^-62 for the last addition and 2^-69.5 inside.
    let r2 = r * r;
    let p = r + r2 * ((E2 + r * E3) + r2 * (E4 + r * E5));
    let j = (n & 127) as usize;
    let t_hi = f64::from_bits(EXP_TABLE[2 * j]);
    let t_lo = f64::from_bits(EXP_TABLE[2 * j + 1]);
    // 2^(j/128) eʳ = (T_hi + T_lo)(1 + p) = T_hi + (T_hi p + T_lo) + T_lo p, the last below 2^-61.5
    // and left out; the two roundings here are each at most 2^-61.
    let c = t_hi * p + t_lo;
    (n >> 7, t_hi, c)
}

/// `eˣ`, to within 0.520 ulp, the same bits on every platform. See the [module](self).
///
/// ```
/// use pantometry_core::math::exp;
/// assert_eq!(exp(0.0), 1.0);
/// assert!((exp(1.0) - std::f64::consts::E).abs() <= f64::EPSILON * std::f64::consts::E);
/// assert_eq!(exp(-746.0), 0.0);
/// assert_eq!(exp(710.0), f64::INFINITY);
/// ```
#[inline]
#[must_use]
pub fn exp(x: f64) -> f64 {
    if x.to_bits() & !(1u64 << 63) >= EXP_MAIN_LIMIT {
        return exp_far(x);
    }
    let (k, t_hi, c) = exp_parts(x);
    // The one rounding of the result, then an exact scaling: k is in [−1022, 1021] and the sum,
    // over 2^-0.43 when k = −1022, keeps the product normal.
    (t_hi + c) * two_to(k)
}

/// [`exp`] for `|x| ≥ 708` and NaN: the thresholds, and the results that are subnormal, in the
/// first normal binade, or within a factor of two of overflow.
#[cold]
#[inline(never)]
fn exp_far(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > EXP_OVERFLOW_BEYOND {
        return f64::INFINITY;
    }
    if x < EXP_UNDERFLOW_BEYOND {
        return 0.0;
    }
    // Here k is in [−1076, 1024].
    let (k, t_hi, c) = exp_parts(x);
    if k > 1023 {
        // 2^1024 is not a double: scale by 2^1023, exactly, and then by 2, which rounds to +∞
        // exactly when the sum rounded to 1 or more, i.e. when the result is MAX + ½ ulp or more.
        return ((t_hi + c) * two_to(1023)) * 2.0;
    }
    if k > -1022 {
        return (t_hi + c) * two_to(k);
    }
    // The result is below 2^-1021. Scale both parts by 2^(k + 1022), exactly — they stay normal
    // — so that the result is y′ · 2^-1022 with y′ < 2.
    let s = two_to(k + 1022);
    let (a, b) = (t_hi * s, c * s);
    let y = a + b;
    if y >= 1.0 {
        // The first normal binade: one rounding, then an exact scaling.
        return y * TWO_M1022;
    }
    // Subnormal. Its grid, 2^-1074, is ulp(1) · 2^-1022, so 1 + y′ rounded once to a double is y′
    // rounded once to the grid. Dekker: z = 1 + a and its error e exactly (1 ≥ a), then
    // z + (e + b), whose one rounding is the only one that matters — the error of e + b is
    // 2^-61 at most, 2^-9 of the grid. z + (e + b) is in [1, 2], so subtracting 1 is exact, and
    // so is the scaling onto the subnormal grid.
    let z = 1.0 + a;
    let e = a - (z - 1.0);
    let w = z + (e + b);
    (w - 1.0) * TWO_M1022
}

/// The bits of `0.703125`, the first `F` and the bottom of the range `m` is reduced into.
const LN_F_FIRST: u64 = 0x3FE6_8000_0000_0000;

/// The bits below `F`'s eight fraction bits: an `F` is a double with these zero.
const LN_CELL: u64 = (1 << 44) - 1;

/// `2⁵²`, which takes a subnormal argument into the normal range exactly.
const TWO_52: f64 = 4_503_599_627_370_496.0;

/// The series' coefficients, `(−1)ⁿ⁺¹/n`.
const L2: f64 = -0.5;
const L3: f64 = 1.0 / 3.0;
const L4: f64 = -0.25;
const L5: f64 = 0.2;
const L6: f64 = -1.0 / 6.0;
const L7: f64 = 1.0 / 7.0;

/// Per `F`, from `0.703125` to `1.40625` in steps of `2⁻⁹` below 1 and `2⁻⁸` above — `F`'s bits
/// are `LN_F_FIRST + (i << 44)` — the bits of `c`, of `(−ln c)_hi`, the multiple of `2⁻³⁵`
/// nearest `−ln c`, and of `(−ln c)_lo`, the double nearest the rest. `c` is `1/F` rounded to
/// nine significant bits — a multiple of `2⁻⁹` below 1, of `2⁻⁸` from 1 up — except in the two
/// cells about 1, `F = 511/512` and `F = 1`, where it is 1.
#[rustfmt::skip]
const LN_TABLE: [[u64; 3]; 257] = [
    [0x3ff6_c000_0000_0000, 0xbfd6_86c8_1e98_0000, 0xbda8_a576_2215_f081],
    [0x3ff6_b000_0000_0000, 0xbfd6_59b5_7300_0000, 0xbdaf_0f94_0ed8_57c7],
    [0x3ff6_a000_0000_0000, 0xbfd6_2c82_f2b8_0000, 0xbd9c_7952_f6f5_f22a],
    [0x3ff6_9000_0000_0000, 0xbfd5_ff30_70a8_0000, 0x3d7b_0b0d_e307_7d7e],
    [0x3ff6_8000_0000_0000, 0xbfd5_d1bd_bf58_0000, 0xbd43_94a1_1b1c_1ee4],
    [0x3ff6_7000_0000_0000, 0xbfd5_a42a_b0f8_0000, 0x3da9_80f3_1d79_6fbe],
    [0x3ff6_6000_0000_0000, 0xbfd5_7677_1748_0000, 0x3da5_2c9d_5b2a_49b0],
    [0x3ff6_5000_0000_0000, 0xbfd5_48a2_c3b0_0000, 0x3da1_6ce9_819c_f7e3],
    [0x3ff6_4000_0000_0000, 0xbfd5_1aad_8730_0000, 0x3da0_3e97_b1b6_14fa],
    [0x3ff6_3000_0000_0000, 0xbfd4_ec97_3260_0000, 0xbd23_4d7a_af04_d104],
    [0x3ff6_2000_0000_0000, 0xbfd4_be5f_9578_0000, 0x3d80_ebe4_966c_d6c1],
    [0x3ff6_1000_0000_0000, 0xbfd4_9006_8040_0000, 0xbd43_a198_00f2_f83a],
    [0x3ff6_0000_0000_0000, 0xbfd4_618b_c220_0000, 0x3dad_09ec_17a4_2642],
    [0x3ff5_f000_0000_0000, 0xbfd4_32ef_2a08_0000, 0x3da8_bf62_5326_2e2b],
    [0x3ff5_e000_0000_0000, 0xbfd4_0430_8688_0000, 0x3d95_81c4_2f3e_d821],
    [0x3ff5_e000_0000_0000, 0xbfd4_0430_8688_0000, 0x3d95_81c4_2f3e_d821],
    [0x3ff5_d000_0000_0000, 0xbfd3_d54f_a5c0_0000, 0xbd9f_70f8_7366_8e58],
    [0x3ff5_c000_0000_0000, 0xbfd3_a64c_5568_0000, 0xbd94_5e9c_72f3_5cd7],
    [0x3ff5_b000_0000_0000, 0xbfd3_7726_62c0_0000, 0x3d63_d286_d58a_7604],
    [0x3ff5_a000_0000_0000, 0xbfd3_47dd_9a98_0000, 0xbd7f_5535_9159_d3fb],
    [0x3ff5_9000_0000_0000, 0xbfd3_1871_c958_0000, 0x3dad_f3d8_2a35_9898],
    [0x3ff5_8000_0000_0000, 0xbfd2_e8e2_bae0_0000, 0xbd91_d309_c2cc_91a8],
    [0x3ff5_7000_0000_0000, 0xbfd2_b930_3ab8_0000, 0xbd83_a493_b4a5_013d],
    [0x3ff5_6000_0000_0000, 0xbfd2_895a_13e0_0000, 0x3d97_95ca_14b6_cfb0],
    [0x3ff5_5000_0000_0000, 0xbfd2_5960_10e0_0000, 0x3d81_38c2_1eed_8ae1],
    [0x3ff5_4000_0000_0000, 0xbfd2_2941_fbd0_0000, 0x3d80_d34b_b7af_584b],
    [0x3ff5_4000_0000_0000, 0xbfd2_2941_fbd0_0000, 0x3d80_d34b_b7af_584b],
    [0x3ff5_3000_0000_0000, 0xbfd1_f8ff_9e48_0000, 0xbd84_5e51_b010_32fa],
    [0x3ff5_2000_0000_0000, 0xbfd1_c898_c168_0000, 0xbd99_9faf_bc68_e754],
    [0x3ff5_1000_0000_0000, 0xbfd1_980d_2dd8_0000, 0x3dae_e484_c585_c9e3],
    [0x3ff5_0000_0000_0000, 0xbfd1_675c_aba8_0000, 0xbdad_3070_1ce6_3eac],
    [0x3ff4_f000_0000_0000, 0xbfd1_3687_0290_0000, 0xbdad_4582_f6cc_531e],
    [0x3ff4_e000_0000_0000, 0xbfd1_058b_f9b0_0000, 0x3d9b_52ae_7605_f54b],
    [0x3ff4_e000_0000_0000, 0xbfd1_058b_f9b0_0000, 0x3d9b_52ae_7605_f54b],
    [0x3ff4_d000_0000_0000, 0xbfd0_d46b_5798_0000, 0xbda5_ba59_03ec_81c4],
    [0x3ff4_c000_0000_0000, 0xbfd0_a324_e270_0000, 0xbdac_871a_fb9f_bd01],
    [0x3ff4_b000_0000_0000, 0xbfd0_71b8_5fd0_0000, 0x3da5_3797_1747_c034],
    [0x3ff4_a000_0000_0000, 0xbfd0_4025_94b8_0000, 0x3da9_7df9_28ec_217a],
    [0x3ff4_9000_0000_0000, 0xbfd0_0e6c_45b0_0000, 0x3da5_7f19_cb95_68ff],
    [0x3ff4_9000_0000_0000, 0xbfd0_0e6c_45b0_0000, 0x3da5_7f19_cb95_68ff],
    [0x3ff4_8000_0000_0000, 0xbfcf_b918_6d60_0000, 0x3d8c_1d57_2aab_993d],
    [0x3ff4_7000_0000_0000, 0xbfcf_550a_5650_0000, 0x3da2_1323_e3a0_9203],
    [0x3ff4_6000_0000_0000, 0xbfce_f0ad_cbe0_0000, 0x3d9d_364d_6f39_0d5e],
    [0x3ff4_5000_0000_0000, 0xbfce_8c02_52b0_0000, 0x3da6_9680_5b80_e8e7],
    [0x3ff4_4000_0000_0000, 0xbfce_2707_6e30_0000, 0x3da4_3468_5855_e000],
    [0x3ff4_4000_0000_0000, 0xbfce_2707_6e30_0000, 0x3da4_3468_5855_e000],
    [0x3ff4_3000_0000_0000, 0xbfcd_c1bc_a0b0_0000, 0x3da0_4e0a_7cb3_ae66],
    [0x3ff4_2000_0000_0000, 0xbfcd_5c21_6b50_0000, 0x3d61_1ba9_1bbc_a682],
    [0x3ff4_1000_0000_0000, 0xbfcc_f635_4e10_0000, 0x3da8_e88e_dc65_f82b],
    [0x3ff4_0000_0000_0000, 0xbfcc_8ff7_c7a0_0000, 0x3da5_9779_4f68_9f84],
    [0x3ff4_0000_0000_0000, 0xbfcc_8ff7_c7a0_0000, 0x3da5_9779_4f68_9f84],
    [0x3ff3_f000_0000_0000, 0xbfcc_2968_5590_0000, 0x3d9f_39fa_e7bd_c714],
    [0x3ff3_e000_0000_0000, 0xbfcb_c286_7430_0000, 0x3d93_994e_b031_8bb8],
    [0x3ff3_d000_0000_0000, 0xbfcb_5b51_9e90_0000, 0x3d62_96e4_5d80_23e6],
    [0x3ff3_d000_0000_0000, 0xbfcb_5b51_9e90_0000, 0x3d62_96e4_5d80_23e6],
    [0x3ff3_c000_0000_0000, 0xbfca_f3c9_4e80_0000, 0xbd77_fe5b_19cc_0327],
    [0x3ff3_b000_0000_0000, 0xbfca_8bec_fc90_0000, 0x3daf_439c_f461_bc8c],
    [0x3ff3_a000_0000_0000, 0xbfca_23bc_1fe0_0000, 0xbd95_ab18_c9b8_8d84],
    [0x3ff3_a000_0000_0000, 0xbfca_23bc_1fe0_0000, 0xbd95_ab18_c9b8_8d84],
    [0x3ff3_9000_0000_0000, 0xbfc9_bb36_2e80_0000, 0x3d90_23e5_5143_9c20],
    [0x3ff3_8000_0000_0000, 0xbfc9_525a_9cf0_0000, 0xbda1_5ad1_d904_c1d5],
    [0x3ff3_7000_0000_0000, 0xbfc8_e928_de90_0000, 0x3dae_4afd_569d_851a],
    [0x3ff3_7000_0000_0000, 0xbfc8_e928_de90_0000, 0x3dae_4afd_569d_851a],
    [0x3ff3_6000_0000_0000, 0xbfc8_7fa0_6520_0000, 0xbd79_2212_0401_2030],
    [0x3ff3_5000_0000_0000, 0xbfc8_15c0_a140_0000, 0xbd9a_bf56_b41b_7f8c],
    [0x3ff3_4000_0000_0000, 0xbfc7_ab89_0210_0000, 0xbd7b_2123_7c6d_65ad],
    [0x3ff3_4000_0000_0000, 0xbfc7_ab89_0210_0000, 0xbd7b_2123_7c6d_65ad],
    [0x3ff3_3000_0000_0000, 0xbfc7_40f8_f540_0000, 0xbd5b_d264_d9bf_9d58],
    [0x3ff3_2000_0000_0000, 0xbfc6_d60f_e720_0000, 0x3da8_b78d_caae_268f],
    [0x3ff3_2000_0000_0000, 0xbfc6_d60f_e720_0000, 0x3da8_b78d_caae_268f],
    [0x3ff3_1000_0000_0000, 0xbfc6_6acd_4270_0000, 0xbd95_6a86_f6ff_1b1e],
    [0x3ff3_0000_0000_0000, 0xbfc5_ff30_70a0_0000, 0xbdae_4f4f_21cf_8828],
    [0x3ff2_f000_0000_0000, 0xbfc5_9338_d9a0_0000, 0x3daf_7de8_b2e9_1554],
    [0x3ff2_f000_0000_0000, 0xbfc5_9338_d9a0_0000, 0x3daf_7de8_b2e9_1554],
    [0x3ff2_e000_0000_0000, 0xbfc5_26e5_e3a0_0000, 0xbd8b_437a_2e40_1d6e],
    [0x3ff2_d000_0000_0000, 0xbfc4_ba36_f3a0_0000, 0x3da6_a86a_9767_e434],
    [0x3ff2_d000_0000_0000, 0xbfc4_ba36_f3a0_0000, 0x3da6_a86a_9767_e434],
    [0x3ff2_c000_0000_0000, 0xbfc4_4d2b_6cd0_0000, 0x3da2_0b86_60b0_9abc],
    [0x3ff2_b000_0000_0000, 0xbfc3_dfc2_b0f0_0000, 0x3d99_ceb1_ab3a_8e7e],
    [0x3ff2_b000_0000_0000, 0xbfc3_dfc2_b0f0_0000, 0x3d99_ceb1_ab3a_8e7e],
    [0x3ff2_a000_0000_0000, 0xbfc3_71fc_2020_0000, 0x3d87_08bc_4326_93aa],
    [0x3ff2_9000_0000_0000, 0xbfc3_03d7_18e0_0000, 0xbda1_ff4b_fa51_8e0a],
    [0x3ff2_9000_0000_0000, 0xbfc3_03d7_18e0_0000, 0xbda1_ff4b_fa51_8e0a],
    [0x3ff2_8000_0000_0000, 0xbfc2_9552_f820_0000, 0x3d35_b967_f447_1dfc],
    [0x3ff2_7000_0000_0000, 0xbfc2_266f_1910_0000, 0x3da6_94d2_0ab8_40e8],
    [0x3ff2_7000_0000_0000, 0xbfc2_266f_1910_0000, 0x3da6_94d2_0ab8_40e8],
    [0x3ff2_6000_0000_0000, 0xbfc1_b72a_d530_0000, 0x3d73_0bfa_df3f_72e3],
    [0x3ff2_5000_0000_0000, 0xbfc1_4785_8460_0000, 0xbdad_0ab1_a288_13e4],
    [0x3ff2_5000_0000_0000, 0xbfc1_4785_8460_0000, 0xbdad_0ab1_a288_13e4],
    [0x3ff2_4000_0000_0000, 0xbfc0_d77e_7cd0_0000, 0xbd71_cb2c_d2ee_2f48],
    [0x3ff2_3000_0000_0000, 0xbfc0_6715_12d0_0000, 0x3da6_9a47_579c_dc0a],
    [0x3ff2_3000_0000_0000, 0xbfc0_6715_12d0_0000, 0x3da6_9a47_579c_dc0a],
    [0x3ff2_2000_0000_0000, 0xbfbf_ec91_31e0_0000, 0x3d90_5515_5746_b998],
    [0x3ff2_1000_0000_0000, 0xbfbf_0a30_c020_0000, 0x3dad_3ab3_3d06_6d1d],
    [0x3ff2_1000_0000_0000, 0xbfbf_0a30_c020_0000, 0x3dad_3ab3_3d06_6d1d],
    [0x3ff2_0000_0000_0000, 0xbfbe_2707_6e20_0000, 0xbda5_e5cb_d3d5_1000],
    [0x3ff1_f000_0000_0000, 0xbfbd_4313_d660_0000, 0xbda9_66ba_bc86_eca9],
    [0x3ff1_f000_0000_0000, 0xbfbd_4313_d660_0000, 0xbda9_66ba_bc86_eca9],
    [0x3ff1_e000_0000_0000, 0xbfbc_5e54_8f60_0000, 0x3d90_e2f3_a8a7_a042],
    [0x3ff1_e000_0000_0000, 0xbfbc_5e54_8f60_0000, 0x3d90_e2f3_a8a7_a042],
    [0x3ff1_d000_0000_0000, 0xbfbb_78c8_2bc0_0000, 0x3dae_24bd_ef78_7310],
    [0x3ff1_c000_0000_0000, 0xbfba_926d_3a40_0000, 0xbda5_aac6_ca17_a455],
    [0x3ff1_c000_0000_0000, 0xbfba_926d_3a40_0000, 0xbda5_aac6_ca17_a455],
    [0x3ff1_b000_0000_0000, 0xbfb9_ab42_4620_0000, 0xbd49_d66d_f661_e3e8],
    [0x3ff1_a000_0000_0000, 0xbfb8_c345_d640_0000, 0x3dac_c9be_14a6_97ab],
    [0x3ff1_a000_0000_0000, 0xbfb8_c345_d640_0000, 0x3dac_c9be_14a6_97ab],
    [0x3ff1_9000_0000_0000, 0xbfb7_da76_6d80_0000, 0x3d93_b4cd_eeed_fcde],
    [0x3ff1_9000_0000_0000, 0xbfb7_da76_6d80_0000, 0x3d93_b4cd_eeed_fcde],
    [0x3ff1_8000_0000_0000, 0xbfb6_f0d2_8ae0_0000, 0xbd95_ad2e_6f92_66e8],
    [0x3ff1_7000_0000_0000, 0xbfb6_0658_a940_0000, 0x3da1_5e78_9c42_2c76],
    [0x3ff1_7000_0000_0000, 0xbfb6_0658_a940_0000, 0x3da1_5e78_9c42_2c76],
    [0x3ff1_6000_0000_0000, 0xbfb5_1b07_3f00_0000, 0xbd98_60fd_a49e_39a2],
    [0x3ff1_6000_0000_0000, 0xbfb5_1b07_3f00_0000, 0xbd98_60fd_a49e_39a2],
    [0x3ff1_5000_0000_0000, 0xbfb4_2edc_bea0_0000, 0xbd99_1bc0_eeea_7c9b],
    [0x3ff1_5000_0000_0000, 0xbfb4_2edc_bea0_0000, 0xbd99_1bc0_eeea_7c9b],
    [0x3ff1_4000_0000_0000, 0xbfb3_41d7_9620_0000, 0x3d90_b8bd_b599_f228],
    [0x3ff1_3000_0000_0000, 0xbfb2_53f6_2f00_0000, 0xbda4_282d_f1f6_d34e],
    [0x3ff1_3000_0000_0000, 0xbfb2_53f6_2f00_0000, 0xbda4_282d_f1f6_d34e],
    [0x3ff1_2000_0000_0000, 0xbfb1_6536_eea0_0000, 0xbd8b_d707_4312_e0ba],
    [0x3ff1_2000_0000_0000, 0xbfb1_6536_eea0_0000, 0xbd8b_d707_4312_e0ba],
    [0x3ff1_1000_0000_0000, 0xbfb0_7598_35a0_0000, 0x3d9c_6e3b_3f92_d666],
    [0x3ff1_0000_0000_0000, 0xbfaf_0a30_c000_0000, 0xbda1_62a6_617c_c971],
    [0x3ff1_0000_0000_0000, 0xbfaf_0a30_c000_0000, 0xbda1_62a6_617c_c971],
    [0x3ff0_f000_0000_0000, 0xbfad_276b_8ac0_0000, 0xbdab_0b52_11e3_c532],
    [0x3ff0_f000_0000_0000, 0xbfad_276b_8ac0_0000, 0xbdab_0b52_11e3_c532],
    [0x3ff0_e000_0000_0000, 0xbfab_42dd_7100_0000, 0xbda9_71be_c28d_14c8],
    [0x3ff0_e000_0000_0000, 0xbfab_42dd_7100_0000, 0xbda9_71be_c28d_14c8],
    [0x3ff0_d000_0000_0000, 0xbfa9_5c83_0ec0_0000, 0xbd91_c7d6_fad0_7403],
    [0x3ff0_d000_0000_0000, 0xbfa9_5c83_0ec0_0000, 0xbd91_c7d6_fad0_7403],
    [0x3ff0_c000_0000_0000, 0xbfa7_7458_f640_0000, 0x3d9a_4607_7396_1abc],
    [0x3ff0_b000_0000_0000, 0xbfa5_8a5b_afc0_0000, 0xbd91_c9a9_18d5_1ea6],
    [0x3ff0_b000_0000_0000, 0xbfa5_8a5b_afc0_0000, 0xbd91_c9a9_18d5_1ea6],
    [0x3ff0_a000_0000_0000, 0xbfa3_9e87_ba00_0000, 0x3d64_2a05_6fea_4dfd],
    [0x3ff0_a000_0000_0000, 0xbfa3_9e87_ba00_0000, 0x3d64_2a05_6fea_4dfd],
    [0x3ff0_9000_0000_0000, 0xbfa1_b0d9_8940_0000, 0x3dac_2680_3d35_d113],
    [0x3ff0_9000_0000_0000, 0xbfa1_b0d9_8940_0000, 0x3dac_2680_3d35_d113],
    [0x3ff0_8000_0000_0000, 0xbf9f_829b_0e80_0000, 0x3d7f_33fe_cc1c_0fb1],
    [0x3ff0_8000_0000_0000, 0xbf9f_829b_0e80_0000, 0x3d7f_33fe_cc1c_0fb1],
    [0x3ff0_7000_0000_0000, 0xbf9b_9fc0_2780_0000, 0xbda7_c8cb_fdea_32dc],
    [0x3ff0_7000_0000_0000, 0xbf9b_9fc0_2780_0000, 0xbda7_c8cb_fdea_32dc],
    [0x3ff0_6000_0000_0000, 0xbf97_b91b_0800_0000, 0x3da5_2772_ab6c_055a],
    [0x3ff0_6000_0000_0000, 0xbf97_b91b_0800_0000, 0x3da5_2772_ab6c_055a],
    [0x3ff0_5000_0000_0000, 0xbf93_cea4_4380_0000, 0x3dac_ad45_8865_ad48],
    [0x3ff0_5000_0000_0000, 0xbf93_cea4_4380_0000, 0x3dac_ad45_8865_ad48],
    [0x3ff0_4000_0000_0000, 0xbf8f_c0a8_b100_0000, 0x3d5f_e0e1_8309_2c59],
    [0x3ff0_4000_0000_0000, 0xbf8f_c0a8_b100_0000, 0x3d5f_e0e1_8309_2c59],
    [0x3ff0_3000_0000_0000, 0xbf87_dc47_6000_0000, 0x3daf_bd62_48b6_bb44],
    [0x3ff0_3000_0000_0000, 0xbf87_dc47_6000_0000, 0x3daf_bd62_48b6_bb44],
    [0x3ff0_2000_0000_0000, 0xbf7f_e02a_6c00_0000, 0x3dad_f30e_e079_12e0],
    [0x3ff0_2000_0000_0000, 0xbf7f_e02a_6c00_0000, 0x3dad_f30e_e079_12e0],
    [0x3ff0_1000_0000_0000, 0xbf6f_f00a_a400_0000, 0x3da4_ef43_fb5f_794b],
    [0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000],
    [0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000],
    [0x3fef_e000_0000_0000, 0x3f70_0805_5a00_0000, 0xbda4_ee99_5034_ce3a],
    [0x3fef_c000_0000_0000, 0x3f80_1015_7600_0000, 0xbdad_c863_b5cc_ce96],
    [0x3fef_a000_0000_0000, 0x3f88_2448_a400_0000, 0xbdad_d755_77da_74f6],
    [0x3fef_8000_0000_0000, 0x3f90_2056_5880_0000, 0x3d93_5847_49f2_3a10],
    [0x3fef_6000_0000_0000, 0x3f94_32a9_2580_0000, 0x3d98_0cc0_9cc9_431c],
    [0x3fef_4000_0000_0000, 0x3f98_4925_2900_0000, 0xbdab_9aa0_ba32_5a0c],
    [0x3fef_2000_0000_0000, 0x3f9c_63d2_ec00_0000, 0x3d94_aaf1_8c7f_3d66],
    [0x3fef_0000_0000_0000, 0x3fa0_415d_8a00_0000, 0xbda8_bbbb_8fe8_c38a],
    [0x3fee_f000_0000_0000, 0x3fa1_49e3_e400_0000, 0x3d46_a33a_b2df_4b82],
    [0x3fee_d000_0000_0000, 0x3fa3_5c8b_fac0_0000, 0xbdae_cf95_0542_0c2a],
    [0x3fee_b000_0000_0000, 0x3fa5_715c_4c00_0000, 0x3d7e_7777_2203_b89d],
    [0x3fee_9000_0000_0000, 0x3fa7_8859_5a40_0000, 0xbd95_108b_0d08_3b3a],
    [0x3fee_7000_0000_0000, 0x3fa9_a187_b580_0000, 0xbd98_4308_b93b_1364],
    [0x3fee_5000_0000_0000, 0x3fab_bceb_fc80_0000, 0xbda7_0bdf_c346_18be],
    [0x3fee_4000_0000_0000, 0x3fac_cb73_cdc0_0000, 0x3dad_b2cb_86dc_13ec],
    [0x3fee_2000_0000_0000, 0x3fae_ea31_c000_0000, 0x3d8a_e1ee_c1b0_36c5],
    [0x3fee_0000_0000_0000, 0x3fb0_8598_b5a0_0000, 0xbd7c_5f97_75c0_2641],
    [0x3fed_e000_0000_0000, 0x3fb1_973b_d140_0000, 0x3d99_559b_4553_e4c3],
    [0x3fed_d000_0000_0000, 0x3fb2_207b_5c80_0000, 0xbd9e_ad89_cc0f_bce1],
    [0x3fed_b000_0000_0000, 0x3fb3_33d7_f820_0000, 0xbd9f_02d2_56d5_0371],
    [0x3fed_9000_0000_0000, 0x3fb4_485e_03e0_0000, 0xbd90_814a_e45c_b655],
    [0x3fed_7000_0000_0000, 0x3fb5_5e10_0500_0000, 0x3dac_0707_5d03_14f2],
    [0x3fed_6000_0000_0000, 0x3fb5_e95a_4da0_0000, 0xbda0_dc69_063c_5d1d],
    [0x3fed_4000_0000_0000, 0x3fb7_00d3_0ae0_0000, 0x3da5_81c1_e8da_99df],
    [0x3fed_2000_0000_0000, 0x3fb8_197e_2f40_0000, 0x3d6c_7e03_73e5_bff8],
    [0x3fed_1000_0000_0000, 0x3fb8_a647_7aa0_0000, 0xbdac_47ae_7ea0_c852],
    [0x3fec_f000_0000_0000, 0x3fb9_c0c3_2d40_0000, 0x3daa_4a90_7ec2_f8f3],
    [0x3fec_e000_0000_0000, 0x3fba_4e76_40c0_0000, 0xbdac_8790_adae_5102],
    [0x3fec_c000_0000_0000, 0x3fbb_6ac8_8da0_0000, 0x3daa_b637_bfea_044c],
    [0x3fec_a000_0000_0000, 0x3fbc_8858_01c0_0000, 0xbd8d_a6e4_b8e6_954d],
    [0x3fec_9000_0000_0000, 0x3fbd_1797_8820_0000, 0x3d79_3643_3b5e_fbef],
    [0x3fec_7000_0000_0000, 0x3fbe_3707_ee40_0000, 0xbdaf_6f09_7b19_8995],
    [0x3fec_6000_0000_0000, 0x3fbe_c739_8300_0000, 0x3da4_223f_9750_19ba],
    [0x3fec_4000_0000_0000, 0x3fbf_e891_39e0_0000, 0xbd90_aa69_ac9f_4216],
    [0x3fec_2000_0000_0000, 0x3fc0_8598_b5a0_0000, 0xbd8c_5f97_75c0_2641],
    [0x3fec_1000_0000_0000, 0x3fc0_ce7e_cdd0_0000, 0xbd99_eb9a_d254_00ac],
    [0x3feb_f000_0000_0000, 0x3fc1_60c8_0250_0000, 0xbda3_613d_2d56_ff62],
    [0x3feb_e000_0000_0000, 0x3fc1_aa2b_7e20_0000, 0x3d9f_b94f_1c88_7132],
    [0x3feb_c000_0000_0000, 0x3fc2_3d71_2a50_0000, 0xbda8_f7f9_6e38_1610],
    [0x3feb_b000_0000_0000, 0x3fc2_8753_bc10_0000, 0x3d8a_ba4a_71ac_9817],
    [0x3feb_9000_0000_0000, 0x3fc3_1b99_4d40_0000, 0xbda6_c1ec_e238_b5f0],
    [0x3feb_8000_0000_0000, 0x3fc3_65fc_b010_0000, 0x3da6_4058_bea0_8d2e],
    [0x3feb_6000_0000_0000, 0x3fc3_fb45_a5a0_0000, 0xbdab_5cd1_d87e_6a35],
    [0x3feb_5000_0000_0000, 0x3fc4_462b_9dd0_0000, 0xbda9_3090_b14e_2361],
    [0x3feb_3000_0000_0000, 0x3fc4_dc7b_8980_0000, 0xbda0_f8e1_927d_4780],
    [0x3feb_2000_0000_0000, 0x3fc5_27e5_e4a0_0000, 0x3d8b_58cf_a395_a5f7],
    [0x3feb_1000_0000_0000, 0x3fc5_737c_c900_0000, 0x3d88_cdd5_3d35_c440],
    [0x3fea_f000_0000_0000, 0x3fc6_0b31_00b0_0000, 0x3d72_8eba_9367_707f],
    [0x3fea_e000_0000_0000, 0x3fc6_574e_be90_0000, 0xbd9f_6630_74d3_c3d2],
    [0x3fea_c000_0000_0000, 0x3fc6_f012_8b70_0000, 0x3da5_aaee_721a_63de],
    [0x3fea_b000_0000_0000, 0x3fc7_3cb9_0750_0000, 0xbd57_59aa_4340_016a],
    [0x3fea_a000_0000_0000, 0x3fc7_898d_8540_0000, 0x3da1_31cc_f7c7_b75e],
    [0x3fea_8000_0000_0000, 0x3fc8_23c1_6550_0000, 0x3d8a_3c1b_b734_c63d],
    [0x3fea_7000_0000_0000, 0x3fc8_7121_3750_0000, 0x3d7d_328e_b42f_9af7],
    [0x3fea_5000_0000_0000, 0x3fc9_0c6d_ba00_0000, 0xbd9a_1935_f577_18d8],
    [0x3fea_4000_0000_0000, 0x3fc9_5a5a_dcf0_0000, 0x3dac_05fc_8a16_2840],
    [0x3fea_3000_0000_0000, 0x3fc9_a877_8df0_0000, 0xbda1_571e_0b82_0279],
    [0x3fea_1000_0000_0000, 0x3fca_4540_82e0_0000, 0x3daa_ac14_ef90_3ee3],
    [0x3fea_0000_0000_0000, 0x3fca_93ed_3c90_0000, 0xbda4_9872_4350_5621],
    [0x3fe9_f000_0000_0000, 0x3fca_e2ca_6f60_0000, 0x3dac_af51_ab5c_a9eb],
    [0x3fe9_d000_0000_0000, 0x3fcb_8117_30c0_0000, 0xbdaf_70b7_cbe2_3194],
    [0x3fe9_c000_0000_0000, 0x3fcb_d087_3840_0000, 0xbda0_9d4b_c459_5413],
    [0x3fe9_b000_0000_0000, 0x3fcc_2028_ab10_0000, 0x3daf_e6d1_f11a_a385],
    [0x3fe9_a000_0000_0000, 0x3fcc_6ffb_c6f0_0000, 0x3d3e_e138_d3a6_9d43],
    [0x3fe9_8000_0000_0000, 0x3fcd_1037_f260_0000, 0x3da5_79ed_6062_9242],
    [0x3fe9_7000_0000_0000, 0x3fcd_60a1_7f90_0000, 0x3d5a_8a47_e40f_7cb2],
    [0x3fe9_6000_0000_0000, 0x3fcd_b13d_b0d0_0000, 0x3da2_2500_d508_ea50],
    [0x3fe9_5000_0000_0000, 0x3fce_020c_c620_0000, 0x3d9a_d5a9_fea4_8dd8],
    [0x3fe9_3000_0000_0000, 0x3fce_a444_9f00_0000, 0x3da2_abd2_2cc6_e654],
    [0x3fe9_2000_0000_0000, 0x3fce_f5ad_e4e0_0000, 0xbd98_00d1_08ab_2ddc],
    [0x3fe9_1000_0000_0000, 0x3fcf_474b_1350_0000, 0xbd90_6eb9_27c7_7ded],
    [0x3fe9_0000_0000_0000, 0x3fcf_991c_6cb0_0000, 0x3d9d_9bcb_ecca_0cdf],
    [0x3fe8_e000_0000_0000, 0x3fd0_1eae_5628_0000, 0xbd93_96f0_8c14_85e9],
    [0x3fe8_d000_0000_0000, 0x3fd0_47e6_0ce0_0000, 0xbd97_c484_1de5_8d02],
    [0x3fe8_c000_0000_0000, 0x3fd0_7138_6050_0000, 0xbda5_3cec_649d_2256],
    [0x3fe8_b000_0000_0000, 0x3fd0_9aa5_72e8_0000, 0xbd93_92bd_787a_32f3],
    [0x3fe8_a000_0000_0000, 0x3fd0_c42d_6760_0000, 0x3d96_2e31_162c_79d6],
    [0x3fe8_8000_0000_0000, 0x3fd1_178e_8228_0000, 0xbd5b_8421_cc74_be04],
    [0x3fe8_7000_0000_0000, 0x3fd1_4167_ef38_0000, 0xbd98_87cf_e1f6_c954],
    [0x3fe8_6000_0000_0000, 0x3fd1_6b5c_cbb0_0000, 0xbda8_2465_3208_5ae9],
    [0x3fe8_5000_0000_0000, 0x3fd1_956d_3b98_0000, 0x3dae_17d2_f73a_d1aa],
    [0x3fe8_4000_0000_0000, 0x3fd1_bf99_6358_0000, 0x3da3_5ca6_ed51_47be],
    [0x3fe8_3000_0000_0000, 0x3fd1_e9e1_6788_0000, 0x3d83_3e8a_8961_ba4d],
    [0x3fe8_2000_0000_0000, 0x3fd2_1445_6d10_0000, 0xbd94_72bc_10a2_dca3],
    [0x3fe8_0000_0000_0000, 0x3fd2_6962_1138_0000, 0xbda9_236c_3e20_a44c],
    [0x3fe7_f000_0000_0000, 0x3fd2_941a_fb18_0000, 0x3d7a_def3_d48c_f1d7],
    [0x3fe7_e000_0000_0000, 0x3fd2_bef0_7ce0_0000, 0xbdab_6563_05b5_affb],
    [0x3fe7_d000_0000_0000, 0x3fd2_e9e2_bce0_0000, 0x3d92_2860_1825_1a3c],
    [0x3fe7_c000_0000_0000, 0x3fd3_14f1_e1d0_0000, 0x3daa_e71d_852c_dec3],
    [0x3fe7_b000_0000_0000, 0x3fd3_401e_12b0_0000, 0xbd93_45f1_cd55_b8a4],
    [0x3fe7_a000_0000_0000, 0x3fd3_6b67_76c0_0000, 0xbd9e_ee91_324f_0e88],
    [0x3fe7_9000_0000_0000, 0x3fd3_96ce_3598_0000, 0x3dad_fa9f_18ea_6726],
    [0x3fe7_8000_0000_0000, 0x3fd3_c252_7730_0000, 0x3da9_8c1d_aa5b_035f],
    [0x3fe7_6000_0000_0000, 0x3fd4_19b4_23d8_0000, 0xbda0_b9c6_f244_dbc8],
    [0x3fe7_5000_0000_0000, 0x3fd4_4591_e050_0000, 0x3dac_fa45_a9db_5b71],
    [0x3fe7_4000_0000_0000, 0x3fd4_718d_c270_0000, 0x3d9c_41b0_63ed_3053],
    [0x3fe7_3000_0000_0000, 0x3fd4_9da7_f3c0_0000, 0xbda9_df09_9964_a169],
    [0x3fe7_2000_0000_0000, 0x3fd4_c9e0_9e18_0000, 0xbd8a_7882_246c_2b63],
    [0x3fe7_1000_0000_0000, 0x3fd4_f637_ebb8_0000, 0x3da4_c07d_4e69_9db7],
    [0x3fe7_0000_0000_0000, 0x3fd5_22ae_0738_0000, 0x3d84_7af9_c205_931d],
    [0x3fe6_f000_0000_0000, 0x3fd5_4f43_1b78_0000, 0x3daf_0d44_aa60_4885],
    [0x3fe6_e000_0000_0000, 0x3fd5_7bf7_53c8_0000, 0x3d8a_3f5b_dbdc_ba82],
    [0x3fe6_d000_0000_0000, 0x3fd5_a8ca_dbc0_0000, 0xbd92_05f1_e6c2_bdfb],
    [0x3fe6_c000_0000_0000, 0x3fd5_d5bd_df58_0000, 0x3d95_f2fa_6afb_adcd],
];

/// `ln x`, to within 0.508 ulp, the same bits on every platform. See the [module](self).
///
/// ```
/// use pantometry_core::math::ln;
/// assert_eq!(ln(1.0).to_bits(), 0.0f64.to_bits());
/// assert!((ln(2.0) - std::f64::consts::LN_2).abs() <= f64::EPSILON * std::f64::consts::LN_2);
/// assert_eq!(ln(0.0), f64::NEG_INFINITY);
/// assert!(ln(-1.0).is_nan());
/// ```
#[inline]
#[must_use]
pub fn ln(x: f64) -> f64 {
    let mut bits = x.to_bits();
    let mut k: i64 = 0;
    // One comparison takes every argument that is not a positive normal double.
    const MIN_NORMAL: u64 = 0x0010_0000_0000_0000;
    const INF: u64 = 0x7FF0_0000_0000_0000;
    if bits.wrapping_sub(MIN_NORMAL) >= INF - MIN_NORMAL {
        if x.is_nan() || bits == INF {
            return x;
        }
        if bits << 1 == 0 {
            return f64::NEG_INFINITY;
        }
        if bits >> 63 != 0 {
            return f64::NAN;
        }
        // Subnormal: 2^52 x is normal and exact.
        bits = (x * TWO_52).to_bits();
        k = -52;
    }
    // x = 2^e m with m in [0.703125, 1.40625): e is how many binades x's bits are above
    // 0.703125's, rounded down, which the arithmetic shift of the difference gives.
    let e = (bits.wrapping_sub(LN_F_FIRST) as i64) >> 52;
    k += e;
    let m_bits = bits.wrapping_sub((e as u64) << 52);
    // F is m rounded to eight fraction bits, half up: spaced 2^-9 below 1 and 2^-8 above, so
    // |m − F| ≤ 2^-10 and 2^-9.
    let f_bits = (m_bits + (1 << 43)) & !LN_CELL;
    let [c, g_hi, g_lo] = LN_TABLE[((f_bits - LN_F_FIRST) >> 44) as usize];
    let (c, g_hi, g_lo) = (
        f64::from_bits(c),
        f64::from_bits(g_hi),
        f64::from_bits(g_lo),
    );
    let m = f64::from_bits(m_bits);
    let big_f = f64::from_bits(f_bits);
    // r = m c − 1, exactly, as (F c − 1) + (m − F) c. F and c have nine significant bits each,
    // so F c is exact and so is subtracting 1 from it; m − F is exact by Sterbenz and has 43 bits
    // at most, so its product with c is exact; and the two are multiples of 2^-61 whose sum is
    // under 2^-8 — at most 1.71 · 2^-9 — so it has 53 bits at most and is exact as well.
    let r = (big_f * c - 1.0) + (m - big_f) * c;
    // k ln2_hi is exact (11 bits times 35), and so is adding (−ln c)_hi: both are multiples of
    // 2^-35 under 2^11.
    let kf = k as f64;
    let h = kf * f64::from_bits(LN2_HI) + g_hi;
    // Dekker: s = h + r and its error, exactly, since h = 0 or h's exponent is at least r's
    // (the table's test checks this cell by cell).
    let s = h + r;
    let s_err = r - (s - h);
    // ln(1 + r) − r to r⁷, by Estrin's scheme in r and r²: the truncation is under |r|⁸/8 and
    // the evaluation's error 3 · 2^-53 of p.
    let w = r * r;
    let p = w * ((L2 + r * L3) + w * ((L4 + r * L5) + w * (L6 + r * L7)));
    // The small parts: k ln2_lo and (−ln c)_lo.
    let lo = kf * f64::from_bits(LN2_LO) + g_lo;
    s + (s_err + (lo + p))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Double-double arithmetic, Dekker (1971), for checking the tables against their closed
    // forms to about 2^-104 relatively. Test-only: nothing above uses it.

    /// `a + b` and its rounding error, exactly, for any `a` and `b` (Knuth's two-sum).
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    }

    /// The same, given `|a| ≥ |b|`.
    fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        (s, b - (s - a))
    }

    /// `a` as two halves of 26 bits or fewer each (Veltkamp's split).
    fn split(a: f64) -> (f64, f64) {
        let c = 134_217_729.0 * a;
        let hi = c - (c - a);
        (hi, a - hi)
    }

    /// `a · b` and its rounding error, exactly, without a fused multiply–add.
    fn two_prod(a: f64, b: f64) -> (f64, f64) {
        let p = a * b;
        let ((ah, al), (bh, bl)) = (split(a), split(b));
        (p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
    }

    /// A double-double: the value is `hi + lo`, `|lo| ≤ ½ ulp(hi)`.
    #[derive(Clone, Copy, Debug)]
    struct Dd(f64, f64);

    impl Dd {
        fn of(a: f64) -> Dd {
            Dd(a, 0.0)
        }
        fn add(self, b: Dd) -> Dd {
            let (s, e) = two_sum(self.0, b.0);
            let (t, f) = two_sum(self.1, b.1);
            let (s, e) = fast_two_sum(s, e + t);
            let (s, e) = fast_two_sum(s, e + f);
            Dd(s, e)
        }
        fn neg(self) -> Dd {
            Dd(-self.0, -self.1)
        }
        fn sub(self, b: Dd) -> Dd {
            self.add(b.neg())
        }
        fn mul(self, b: Dd) -> Dd {
            let (p, e) = two_prod(self.0, b.0);
            let (p, e) = fast_two_sum(p, e + (self.0 * b.1 + self.1 * b.0));
            Dd(p, e)
        }
        fn div(self, b: Dd) -> Dd {
            let q1 = self.0 / b.0;
            let r = self.sub(b.mul(Dd::of(q1)));
            let q2 = r.0 / b.0;
            let r = r.sub(b.mul(Dd::of(q2)));
            let q3 = r.0 / b.0;
            let (q1, q2) = fast_two_sum(q1, q2);
            Dd(q1, q2).add(Dd::of(q3))
        }
        /// Times a power of two, exactly.
        fn scale(self, p: f64) -> Dd {
            Dd(self.0 * p, self.1 * p)
        }
        /// The value, rounded once (to within the double-double's own error).
        fn value(self) -> f64 {
            self.0 + self.1
        }
    }

    /// `|a − b|` as a double, for two double-doubles that are close.
    fn gap(a: Dd, b: Dd) -> f64 {
        a.sub(b).value().abs()
    }

    /// The spacing of doubles at the normal double `v`.
    fn ulp(v: f64) -> f64 {
        let e = ((v.to_bits() >> 52) & 0x7FF) as i64;
        f64::from_bits(((e - 52).max(1) as u64) << 52)
    }

    /// `ln 2 = Σ 1/(n 2ⁿ)`, in double-double, to about 2^-104: a series this code does not use.
    fn ln2_dd() -> Dd {
        let mut sum = Dd::of(0.0);
        let mut p = 1.0;
        for n in 1..=120 {
            p *= 0.5;
            sum = sum.add(Dd::of(1.0).div(Dd::of(n as f64)).scale(p));
        }
        sum
    }

    /// `ln F = 2 atanh((F − 1)/(F + 1))`, in double-double, for `F` with few bits (`F ± 1` exact):
    /// the series `2 Σ s²ⁱ⁺¹/(2i + 1)`, `|s| < 0.18`, to 40 terms.
    fn ln_dd(big_f: f64) -> Dd {
        let s = Dd::of(big_f - 1.0).div(Dd::of(big_f + 1.0));
        let s2 = s.mul(s);
        let (mut term, mut sum) = (s, s);
        for i in 1..40 {
            term = term.mul(s2);
            sum = sum.add(term.div(Dd::of((2 * i + 1) as f64)));
        }
        sum.scale(2.0)
    }

    /// `eˣ` by its Taylor series in double-double, for `|x| < 1`, to 40 terms.
    fn exp_dd(x: Dd) -> Dd {
        let (mut term, mut sum) = (Dd::of(1.0), Dd::of(1.0));
        for n in 1..40 {
            term = term.mul(x).div(Dd::of(n as f64));
            sum = sum.add(term);
        }
        sum
    }

    fn exp_entry(j: usize) -> Dd {
        Dd(
            f64::from_bits(EXP_TABLE[2 * j]),
            f64::from_bits(EXP_TABLE[2 * j + 1]),
        )
    }

    fn ln_f(i: usize) -> f64 {
        f64::from_bits(LN_F_FIRST + ((i as u64) << 44))
    }

    /// **`ln 2`'s split, `128/ln 2` and the shift are the doubles they say.** `LN2_HI` is a
    /// multiple of `2⁻³⁵` within `2⁻³⁶` of `ln 2`; `LN2_HI + LN2_LO` is `ln 2` to within
    /// `LN2_LO`'s own rounding; `INV_STEP` is the double nearest `128/ln 2`. The reference is
    /// `Σ 1/(n 2ⁿ)` summed here in double-double, not anything the script wrote.
    #[test]
    fn the_constants_are_their_closed_forms() {
        let ln2 = ln2_dd();
        assert_eq!(ln2.value(), std::f64::consts::LN_2, "the series is ln 2");
        let (hi, lo) = (f64::from_bits(LN2_HI), f64::from_bits(LN2_LO));
        let scaled = hi * 2f64.powi(35);
        assert_eq!(scaled, scaled.trunc(), "LN2_HI is a multiple of 2^-35");
        assert!(
            gap(Dd::of(hi), ln2) <= 2f64.powi(-36),
            "LN2_HI is the nearest"
        );
        let off = gap(Dd(hi, lo), ln2);
        assert!(
            off <= 0.5 * ulp(lo.abs()) + 2f64.powi(-100),
            "LN2_HI + LN2_LO is {off:e} from ln 2"
        );
        let inv = Dd::of(128.0).div(ln2);
        let inv_step = f64::from_bits(INV_STEP);
        assert!(
            gap(Dd::of(inv_step), inv) <= 0.5 * ulp(inv_step),
            "INV_STEP is the nearest double to 128/ln 2"
        );
        assert_eq!(SHIFT, 1.5 * 2f64.powi(52));
        assert_eq!(TWO_52, 2f64.powi(52));
        assert_eq!(TWO_M1022, 2f64.powi(-1022));
        assert_eq!(f64::from_bits(EXP_MAIN_LIMIT), 708.0);
        assert_eq!(f64::from_bits(LN_F_FIRST), 0.703125);
    }

    /// **The `2^(j/128)` table is its closed form.** Each entry against `e^(j ln 2/128)` summed in
    /// double-double here — `T_hi` the nearest double, `T_hi + T_lo` within `T_lo`'s own rounding
    /// — and, as products that are exact in double-double, `T_j · T_(128−j) = 2` and
    /// `T_j² = T_2j` (or `2 T_(2j−128)`).
    #[test]
    fn the_exp_table_is_its_closed_form() {
        let ln2 = ln2_dd();
        assert_eq!((exp_entry(0).0, exp_entry(0).1), (1.0, 0.0));
        for j in 0..128 {
            let t = exp_entry(j);
            let exact = exp_dd(ln2.mul(Dd::of(j as f64)).scale(1.0 / 128.0));
            assert!(
                gap(Dd::of(t.0), exact) <= 0.5 * ulp(t.0),
                "T_hi[{j}] is not the nearest double"
            );
            let off = gap(t, exact);
            assert!(
                off <= 0.5 * ulp(t.1.abs().max(f64::MIN_POSITIVE)) + 2f64.powi(-102),
                "T[{j}] is {off:e} from 2^({j}/128)"
            );
            let two = if j == 0 {
                Dd::of(2.0)
            } else {
                t.mul(exp_entry(128 - j))
            };
            assert!(
                gap(two, Dd::of(2.0)) <= 2f64.powi(-102),
                "T[{j}] T[128 − {j}]"
            );
            let square = t.mul(t);
            let target = if 2 * j < 128 {
                exp_entry(2 * j)
            } else {
                exp_entry(2 * j - 128).scale(2.0)
            };
            assert!(gap(square, target) <= 2f64.powi(-101), "T[{j}]²");
        }
    }

    /// **The `ln` table is its closed form, and holds what `ln` assumes of it.** Per `F`:
    ///
    /// - `c` has nine significant bits — a multiple of `2⁻⁹` below 1, of `2⁻⁸` from 1 up — and is
    ///   the nearest such to `1/F`, checked exactly since `c F` is; or is 1, in the two cells
    ///   about 1 and only there;
    /// - over the cell, `|m c − 1| < 2⁻⁸`, computed exactly at both edges: what makes `r` exact;
    /// - where `c ≠ 1`, `(−ln c)_hi`'s exponent is at least that bound's: what makes Dekker's sum
    ///   exact for `k = 0`;
    /// - `(−ln c)_hi` is the multiple of `2⁻³⁵` nearest `−ln c`, and `(−ln c)_hi + (−ln c)_lo`
    ///   is within `(−ln c)_lo`'s own rounding of `−2 atanh((c−1)/(c+1))`, summed here in
    ///   double-double.
    ///
    /// And across the table, every pair with `c_a = 2 c_b` gives `ln 2`.
    #[test]
    fn the_ln_table_is_its_closed_form() {
        assert_eq!(ln_f(0), 0.703125);
        assert_eq!(ln_f(151), 511.0 / 512.0);
        assert_eq!(ln_f(152), 1.0);
        assert_eq!(ln_f(256), 1.40625);
        let exponent = |v: f64| ((v.to_bits() >> 52) & 0x7FF) as i64;
        let mut widest = 0.0f64;
        for (i, &[c, g_hi, g_lo]) in LN_TABLE.iter().enumerate() {
            let big_f = ln_f(i);
            let (c, g_hi, g_lo) = (
                f64::from_bits(c),
                f64::from_bits(g_hi),
                f64::from_bits(g_lo),
            );
            let grid = if c < 1.0 {
                2f64.powi(-9)
            } else {
                2f64.powi(-8)
            };
            assert_eq!(
                (c / grid).fract(),
                0.0,
                "c = {c} for F = {big_f} is off its grid"
            );
            if i == 151 || i == 152 {
                assert_eq!(c, 1.0, "the cells about 1");
                assert_eq!((g_hi, g_lo), (0.0, 0.0));
            } else {
                // |c F − 1| ≤ ½ grid · F, the nearest: exact, since c F has 18 bits.
                let miss = (c * big_f - 1.0).abs();
                assert!(
                    miss <= 0.5 * grid * big_f,
                    "c is not the nearest to 1/{big_f}"
                );
                assert_ne!(c, 1.0, "c = 1 only about 1");
            }
            // The cell: m from F less half the spacing below it to F plus half the spacing above,
            // within [0.703125, 1.40625).
            let below = if big_f <= 1.0 {
                2f64.powi(-9)
            } else {
                2f64.powi(-8)
            };
            let above = if big_f < 1.0 {
                2f64.powi(-9)
            } else {
                2f64.powi(-8)
            };
            let lo = if i == 0 { big_f } else { big_f - below / 2.0 };
            let hi = if i == 256 { big_f } else { big_f + above / 2.0 };
            let r_max = (lo * c - 1.0).abs().max((hi * c - 1.0).abs());
            widest = widest.max(r_max);
            if c != 1.0 {
                assert!(
                    r_max < 2f64.powi(-8),
                    "r is not exact in F = {big_f}'s cell"
                );
                assert!(
                    exponent(g_hi) >= exponent(r_max),
                    "Dekker's sum is not exact in F = {big_f}'s cell"
                );
            } else {
                assert!(r_max <= 3.0 * 2f64.powi(-10), "r = m − 1 about 1");
            }
            let exact = ln_dd(c).neg();
            assert_eq!((g_hi * 2f64.powi(35)).fract(), 0.0, "(−ln {c})_hi on 2^-35");
            assert!(
                gap(Dd::of(g_hi), exact) <= 2f64.powi(-36),
                "(−ln {c})_hi is not the nearest"
            );
            let off = gap(Dd(g_hi, g_lo), exact);
            let allowed = if g_lo == 0.0 {
                2f64.powi(-104)
            } else {
                0.5 * ulp(g_lo.abs()) + 2f64.powi(-104)
            };
            assert!(off <= allowed, "−ln {c} is {off:e} out");
        }
        assert!(widest <= 1.71 * 2f64.powi(-9), "|r| reaches {widest:e}");
        let entry = |i: usize| {
            let [c, hi, lo] = LN_TABLE[i];
            (
                f64::from_bits(c),
                Dd(f64::from_bits(hi), f64::from_bits(lo)),
            )
        };
        let mut pairs = 0;
        for a in 0..257 {
            for b in 0..257 {
                let ((ca, ga), (cb, gb)) = (entry(a), entry(b));
                if ca == 2.0 * cb {
                    pairs += 1;
                    assert!(
                        gap(gb.sub(ga), ln2_dd()) <= 2f64.powi(-88),
                        "c = {ca} and {cb}"
                    );
                }
            }
        }
        assert!(pairs > 0, "no pair checked");
    }

    /// One double up, and one down.
    fn up(x: f64) -> f64 {
        if x == 0.0 {
            f64::from_bits(1)
        } else if x > 0.0 {
            f64::from_bits(x.to_bits() + 1)
        } else {
            f64::from_bits(x.to_bits() - 1)
        }
    }

    fn down(x: f64) -> f64 {
        -up(-x)
    }

    /// **Every special value of `exp` is the standard's**, and the two thresholds are where
    /// they are, to the double: `709.782712893384` (`0x40862e42fefa39ef`) is the last finite
    /// result's argument, and `−745.1332191019411` (`0xc0874910d52d3051`) the last that gives the
    /// smallest subnormal rather than zero — both found with mpmath, and both what rounding the
    /// exact value would give.
    #[test]
    fn exp_of_every_special_value_is_the_standards() {
        assert!(exp(f64::NAN).is_nan());
        assert!(exp(-f64::NAN).is_nan());
        assert_eq!(exp(f64::INFINITY), f64::INFINITY);
        assert_eq!(exp(f64::NEG_INFINITY).to_bits(), 0, "+0");
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(exp(-0.0), 1.0);
        for b in [
            1u64,
            0x000F_FFFF_FFFF_FFFF,
            0x0010_0000_0000_0000,
            0x3C80_0000_0000_0000,
        ] {
            let x = f64::from_bits(b);
            assert_eq!(exp(x), 1.0, "exp({x:e})");
            assert_eq!(exp(-x), 1.0, "exp(-{x:e})");
        }
        let last_finite = f64::from_bits(0x4086_2E42_FEFA_39EF);
        assert!(exp(last_finite).is_finite(), "the largest result");
        assert_eq!(exp(up(last_finite)), f64::INFINITY);
        for x in [709.8, 709.800_000_1, 710.0, 1e300, f64::MAX] {
            assert_eq!(exp(x), f64::INFINITY, "exp({x})");
        }
        let last_nonzero = f64::from_bits(0xC087_4910_D52D_3051);
        assert_eq!(exp(last_nonzero).to_bits(), 1, "the smallest subnormal");
        assert_eq!(exp(down(last_nonzero)).to_bits(), 0, "+0 past it");
        for x in [-745.2, -745.200_000_1, -746.0, -1e300, f64::MIN] {
            assert_eq!(exp(x).to_bits(), 0, "exp({x})");
        }
        // Subnormal results are subnormal, and the first normal binade is reached.
        let s = exp(-740.0);
        assert!(s > 0.0 && s < f64::MIN_POSITIVE);
        assert!(exp(-708.0) >= f64::MIN_POSITIVE);
    }

    /// **Every special value of `ln` is the standard's**: NaN, `+∞`, `−∞`, both zeros, negative
    /// numbers down to the smallest subnormal, and the extremes of the finite range — the
    /// largest double and the smallest subnormal, `ln 2⁻¹⁰⁷⁴ = −744.44007…`.
    #[test]
    fn ln_of_every_special_value_is_the_standards() {
        assert!(ln(f64::NAN).is_nan());
        assert!(ln(-f64::NAN).is_nan());
        assert_eq!(ln(f64::INFINITY), f64::INFINITY);
        assert!(ln(f64::NEG_INFINITY).is_nan());
        assert_eq!(ln(0.0), f64::NEG_INFINITY);
        assert_eq!(ln(-0.0), f64::NEG_INFINITY);
        for x in [-1.0, -f64::MIN_POSITIVE, -f64::from_bits(1), f64::MIN, -0.5] {
            assert!(ln(x).is_nan(), "ln({x:e})");
        }
        // -1074 ln 2 and 1024 ln 2 less a bit: each within its ulp, from the split of ln 2.
        let ln2 = ln2_dd();
        let smallest = ln(f64::from_bits(1));
        let exact = ln2.mul(Dd::of(-1074.0));
        assert!(
            gap(Dd::of(smallest), exact) <= 0.501 * ulp(744.0),
            "ln 2^-1074"
        );
        let largest = ln(f64::MAX);
        assert!(
            (largest - 709.782_712_893_384).abs() <= ulp(709.0),
            "ln MAX"
        );
        for e in -1074..1024 {
            let x = 2f64.powi(e);
            let exact = ln2.mul(Dd::of(f64::from(e)));
            let y = ln(x);
            if e == 0 {
                assert_eq!(y.to_bits(), 0);
            } else {
                assert!(
                    gap(Dd::of(y), exact) <= 0.501 * ulp(y.abs()),
                    "ln 2^{e} = {y:e}"
                );
            }
        }
    }

    /// **`exp(±0) = 1` and `ln 1 = +0`, exactly**, the sign of the zero included.
    #[test]
    fn exp_of_zero_and_ln_of_one_are_exact() {
        assert_eq!(exp(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(exp(-0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(ln(1.0).to_bits(), 0.0f64.to_bits());
    }

    fn mix(seed: u64, i: u64) -> u64 {
        let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A bit pattern uniform in `[lo, hi]`.
    fn between(seed: u64, i: u64, lo: u64, hi: u64) -> u64 {
        lo + mix(seed, i) % (hi - lo + 1)
    }

    /// The bounds the identities are derived from: [`exp`]'s and [`ln`]'s largest, as fractions
    /// of `2⁻⁵²` relatively, since an ulp is at most `2⁻⁵²` of the value.
    const B_EXP: f64 = 0.520;
    const B_LN: f64 = 0.508;
    const EPS: f64 = f64::EPSILON;

    /// **`exp(ln x) = x`**, for `x` log-uniform over `[2⁻¹⁰⁰⁰, MAX]`: `ln x` is off by
    /// `B_ln ε |ln x|`, which `exp` turns into that relative error of the result, and `exp`
    /// adds `B_exp ε`; so `(B_ln |ln x| + B_exp) ε` relatively, and the second-order terms, under
    /// `ε²|ln x|²`, are covered by one more `ε`.
    #[test]
    fn exp_undoes_ln() {
        let lo = 2f64.powi(-1000).to_bits();
        for i in 0..20_000 {
            let x = f64::from_bits(between(21, i, lo, f64::MAX.to_bits()));
            let l = ln(x);
            let y = exp(l);
            let allowed = (B_LN * l.abs() + B_EXP + EPS) * EPS * x;
            assert!((y - x).abs() <= allowed, "exp(ln {x:e}) = {y:e}");
        }
    }

    /// **`ln(exp x) = x`**, for `x` uniform in `[−708, 709]`: `exp`'s relative error `B_exp ε`
    /// is an absolute error `B_exp ε` after `ln`, and `ln` adds `B_ln ε |x|`.
    #[test]
    fn ln_undoes_exp() {
        for i in 0..20_000 {
            let x = -708.0 + 1417.0 * ((mix(22, i) >> 11) as f64 / (1u64 << 53) as f64);
            let y = ln(exp(x));
            let allowed = (B_EXP + B_LN * x.abs()) * EPS * (1.0 + EPS) + EPS * EPS;
            assert!((y - x).abs() <= allowed, "ln(exp {x:e}) = {y:e}");
        }
    }

    /// **`ln(xy) = ln x + ln y`**, for `x` and `y` of 26 significant bits, so that `xy` is
    /// exact, log-uniform over `[2⁻⁵⁰⁰, 2⁵⁰⁰]`: each of the three logarithms is off by at most
    /// `B_ln ε` of itself, and the sum rounds once more.
    #[test]
    fn ln_turns_products_into_sums() {
        let (lo, hi) = (2f64.powi(-500).to_bits(), 2f64.powi(500).to_bits());
        let mask = !((1u64 << 27) - 1);
        for i in 0..20_000 {
            let x = f64::from_bits(between(23, i, lo, hi) & mask);
            let y = f64::from_bits(between(24, i, lo, hi) & mask);
            let (lx, ly, lxy) = (ln(x), ln(y), ln(x * y));
            let sum = lx + ly;
            let allowed = B_LN * EPS * (lx.abs() + ly.abs() + lxy.abs()) + 0.5 * ulp(sum.abs());
            assert!((lxy - sum).abs() <= allowed, "ln({x:e} · {y:e})");
        }
    }

    /// Whether `f` is nondecreasing over `xs`, which are increasing.
    fn rises(f: fn(f64) -> f64, xs: &[f64]) -> Result<(), String> {
        for w in xs.windows(2) {
            let (a, b) = (f(w[0]), f(w[1]));
            if b < a {
                return Err(format!("f({:e}) = {a:e} > f({:e}) = {b:e}", w[0], w[1]));
            }
        }
        Ok(())
    }

    /// `reach` doubles either side of `at`, in order.
    fn about(at: f64, reach: usize) -> Vec<f64> {
        let mut v = Vec::with_capacity(2 * reach + 1);
        let mut x = at;
        for _ in 0..reach {
            x = down(x);
            v.push(x);
        }
        v.reverse();
        v.push(at);
        x = at;
        for _ in 0..reach {
            x = up(x);
            v.push(x);
        }
        v
    }

    /// **Both are monotone where their methods have seams**: `exp` across every `j` boundary
    /// `(n + ½) ln 2/128` for `|n| ≤ 4096`, across every `k` boundary, at the edges of its
    /// paths (±708, the subnormal and first normal binades, both thresholds); `ln` on `2¹⁷`
    /// doubles about 1, across every cell edge for `k = −1, 0, 1`, at `0.703125 · 2ᵉ` in every
    /// binade, and where subnormals meet normals; and both on sorted grids over their domains.
    #[test]
    fn both_are_monotone() {
        let step = std::f64::consts::LN_2 / 128.0;
        for n in -4096i64..=4096 {
            rises(exp, &about((n as f64 + 0.5) * step, 6)).unwrap();
        }
        for k in -1075i64..=1024 {
            rises(exp, &about((k as f64 * 128.0 + 0.5) * step, 6)).unwrap();
        }
        for at in [
            -708.0,
            708.0,
            -1022.0 * std::f64::consts::LN_2,
            -1021.0 * std::f64::consts::LN_2,
            709.782_712_893_384,
            -745.133_219_101_941_1,
        ] {
            rises(exp, &about(at, 64)).unwrap();
        }
        let mut grid: Vec<f64> = (0..100_000)
            .map(|i| -745.2 + 1455.0 * ((mix(31, i) >> 11) as f64 / (1u64 << 53) as f64))
            .collect();
        grid.sort_by(f64::total_cmp);
        rises(exp, &grid).unwrap();

        rises(ln, &about(1.0, 65_536)).unwrap();
        for k in [-1i32, 0, 1] {
            for i in 0..256u64 {
                let edge = f64::from_bits(LN_F_FIRST + (i << 44) + (1 << 43)) * 2f64.powi(k);
                rises(ln, &about(edge, 8)).unwrap();
            }
        }
        for e in -1022..1023 {
            rises(ln, &about(0.703125 * 2f64.powi(e), 4)).unwrap();
        }
        rises(ln, &about(f64::MIN_POSITIVE, 64)).unwrap();
        let mut bits: Vec<u64> = (0..100_000)
            .map(|i| between(32, i, 1, f64::MAX.to_bits()))
            .collect();
        bits.sort_unstable();
        let grid: Vec<f64> = bits.into_iter().map(f64::from_bits).collect();
        rises(ln, &grid).unwrap();
    }

    /// FNV-1a over words.
    fn fnv1a(words: impl Iterator<Item = u64>) -> u64 {
        let mut h: u64 = 0xCBF2_9CE4_8422_2325;
        for w in words {
            for byte in w.to_le_bytes() {
                h ^= u64::from(byte);
                h = h.wrapping_mul(0x0000_0100_0000_01B3);
            }
        }
        h
    }

    /// **The bits are pinned**: an FNV-1a digest of `exp` at 8192 arguments over its whole range,
    /// both signs, and of `ln` at 8192 over every positive double, subnormals included, and both
    /// at their special values. A NaN, whose payload Rust does not fix, is hashed as one word.
    /// CI runs this on Linux, macOS on arm64, Windows, release with LTO and `wasm32-wasip1`, and
    /// the digest is the claim that every one of them computes the same function.
    ///
    /// **This is a new function's digest.** The implementation it pins replaced the first one,
    /// whose own digest — `0x150a0cfa04b62b07`, over other points — pinned another function: the
    /// two differ by an ulp at some arguments, and a digest is not a tolerance.
    #[test]
    fn the_bits_are_pinned() {
        let word = |y: f64| {
            if y.is_nan() {
                0x7FF8_0000_0000_0000
            } else {
                y.to_bits()
            }
        };
        let exp_args = (0..8192u64).map(|i| {
            let b = between(41, i, 2f64.powi(-60).to_bits(), 746.0f64.to_bits());
            f64::from_bits(b | ((i & 1) << 63))
        });
        let ln_args =
            (0..8192u64).map(|i| f64::from_bits(between(42, i, 1, 0x7FEF_FFFF_FFFF_FFFF)));
        let specials = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            f64::MAX,
            709.782_712_893_384,
            -745.133_219_101_941_1,
            -708.5,
        ];
        let digest = fnv1a(
            exp_args
                .chain(specials)
                .map(|x| word(exp(x)))
                .chain(ln_args.chain(specials).map(|x| word(ln(x)))),
        );
        assert_eq!(digest, 0x8453_2fb8_80d9_8a57, "digest {digest:#018x}");
    }
}
