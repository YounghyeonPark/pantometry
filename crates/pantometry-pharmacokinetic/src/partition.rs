//! Tissue partition coefficients from a tissue's composition: Rodgers and Rowland's equation for
//! a neutral molecule.
//!
//! [`CompartmentModel::perfuse`](crate::CompartmentModel::perfuse) takes a `K_p` and says nothing
//! about where it came from. This module is one answer, the mechanistic one a physiologically
//! based model is usually built with: a tissue is so much water, so much neutral lipid, so much
//! neutral phospholipid and so much binding protein, and a drug that is not ionised at
//! physiological pH dissolves in the water, partitions into the lipids as it would into an organic
//! solvent, and binds the protein as it binds the same protein in plasma.
//!
//! ```text
//! Kpu = f_W + P·f_NL + (0.3P + 0.7)·f_NP
//!     + (1/f_u − 1 − P·f_NL,P − (0.3P + 0.7)·f_NP,P) · [PR]_T/[PR]_P
//! K_p = f_u · Kpu
//! ```
//!
//! That is Eq. 14 of Rodgers and Rowland 2006 with both of its ionisation factors `X` and `Y` at
//! one, which is what they are for a neutral molecule (p. 1240). Read term by term:
//!
//! - **`f_W`**: unbound drug in the tissue's water, at the unbound plasma concentration. The
//!   equation has it as `f_EW + f_IW`, extracellular and intracellular apart, and for a neutral
//!   molecule both are at the unbound plasma concentration, so only their sum enters. Every
//!   fraction is of the tissue's volume.
//! - **`P·f_NL`**: neutral lipid, at `P` times the concentration of the water beside it. `P` is
//!   the n-octanol:water partition coefficient — the antilog, not log P — in every tissue but
//!   adipose, where Rodgers and Rowland take a vegetable oil:water coefficient instead (2005,
//!   p. 1263). It is one argument here, so the caller decides which coefficient a tissue gets.
//! - **`(0.3P + 0.7)·f_NP`**: neutral phospholipid, taken to behave as 30% neutral lipid and 70%
//!   water (2005, Eq. 15).
//! - **The last term**: binding to the tissue's share of the plasma's binding protein, which for a
//!   neutral molecule is lipoprotein (2006, p. 1241). The bracket is the plasma's bound-to-unbound
//!   ratio, `1/f_u − 1`, less the part of it the plasma's own lipids already account for; what is
//!   left is put down to protein, and the tissue has `[PR]_T/[PR]_P` of the plasma's protein.
//!
//! `K_p` is then tissue to **plasma**, total over total, because `Kpu` is over the *unbound*
//! plasma concentration. A perfusion's `K_p` in this crate is tissue to *blood*, which is
//! `K_p / (B:P)` with [`blood_to_plasma`] — the blood pool in a model here holds blood, and the
//! well-stirred tissue equilibrates with what leaves it.
//!
//! # What it is checked against
//!
//! The equation is a sum of terms each linear in `P` or constant, so it has exact limits and they
//! are what the tests hold it to: as `P → 0` only the water and the phospholipid's watery 70%
//! remain; `Kpu` is a straight line in `P` whose slope is the tissue's lipid less its protein's
//! share of the plasma's; at `f_u = 1` with no plasma lipid the protein term is zero; and `B:P` is
//! one at a haematocrit of zero and `K_rbc:p` at a haematocrit of one. None of that says the
//! equation is *right* — it is a model, and its authors report 84% of 370 predicted rat `Kpu`
//! values within a factor of three of the measured ones, across every class of drug they tried
//! (2006, Table 4) — only that this is the equation they wrote.
//!
//! # The compositions, and where each number is printed
//!
//! The constants below are **rat** tissues. Rodgers and Rowland had no human table, and a human
//! model built on them is using rat composition as a proxy, which is the usual practice and is a
//! choice the caller is making.
//!
//! - **Neutral lipid, neutral phospholipid, water, lipoprotein ratio**: Rodgers and Rowland 2006,
//!   Table 1 (p. 1242). The three fractions are the columns adjusted for residual blood.
//! - **Plasma neutral lipid and neutral phospholipid, 0.0023 and 0.0013**: Rodgers and Rowland
//!   2006, p. 1241, in the text and not in any table.
//! - **Blood cells**: Rodgers, Leahy and Rowland 2005, Table 1.
//!
//! Coutinho et al. 2024, Table 3, open access under CC BY 4.0, has the same lipids and lipoprotein
//! ratios for the four tissues here, and is the copy a reader without a subscription can check
//! them against.
//!
//! **The water is the 2006 table's, and that was measured rather than assumed.** The 2006 table
//! has one column of total tissue water, adjusted for residual blood; the 2005 table has the
//! extracellular and intracellular water apart and not adjusted, and Coutinho et al. take theirs
//! from the 2005 table. The two differ — liver 0.642 against 0.161 + 0.573 = 0.734 — and the
//! authors' own predictions decide it: their Table 5 lists `Kpu` predicted by this equation for
//! three neutral drugs, and with the 2006 water this module reproduces ethoxybenzamide's to within
//! their printed digits, where the 2005 water puts ftorafur's liver 12% high. That reproduction is
//! a test below.
//!
//! **Adipose is printed transposed, and is not transposed here.** The 2006 table gives adipose a
//! neutral lipid of 0.0016 and a neutral phospholipid of 0.853, and Coutinho et al. reproduce it.
//! The 2005 table has the same tissue the other way round, 0.853 and 0.0016, and so does physics:
//! adipose tissue is mostly triglyceride, which is neutral lipid. [`RAT_ADIPOSE`] carries 0.853 as
//! neutral lipid. The choice is not cosmetic — swapped, the lipid term at `P = 100` falls from 85
//! to 26 — and it is not only a reading either: the authors' own Table 5 adipose predictions come
//! out of the equation with 0.853 as neutral lipid, and the other way round ethoxybenzamide's is
//! 0.878 against the printed 0.48.
//!
//! **Not in either paper:** a haematocrit. Both cite one and print none, so [`blood_to_plasma`]
//! takes it as an argument rather than this module guessing.
//!
//! # References
//!
//! - Rodgers T, Leahy D, Rowland M. Physiologically based pharmacokinetic modeling 1: predicting
//!   the tissue distribution of moderate-to-strong bases. *J Pharm Sci* 2005;94(6):1259–1276.
//!   doi:10.1002/jps.20322.
//! - Rodgers T, Rowland M. Physiologically based pharmacokinetic modelling 2: predicting the tissue
//!   distribution of acids, very weak bases, neutrals and zwitterions. *J Pharm Sci*
//!   2006;95(6):1238–1257. doi:10.1002/jps.20502.
//! - Coutinho AL, Cristofoletti R, Wu F, Al Shoyaib A, Dressman J, Polli JE. Relative performance
//!   of volume of distribution prediction methods for lipophilic drugs with uncertainty in LogP
//!   value. *Pharm Res* 2024;41(6):1121–1138. doi:10.1007/s11095-024-03703-4. CC BY 4.0.

/// The part of a tissue's volume that is each thing a neutral molecule can be in.
///
/// Every field is a fraction of the tissue's volume except [`lipoprotein_ratio`], which is the
/// tissue's concentration of lipoprotein over the plasma's.
///
/// [`lipoprotein_ratio`]: TissueComposition::lipoprotein_ratio
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TissueComposition {
    /// Neutral lipid — triglyceride, cholesterol and its esters. `f_NL`.
    pub neutral_lipid: f64,
    /// Phospholipid with no net charge, mostly phosphatidylcholine and -ethanolamine. `f_NP`.
    pub neutral_phospholipid: f64,
    /// All of the tissue's water, inside the cells and out. `f_EW + f_IW`.
    pub water: f64,
    /// Lipoprotein in the tissue over lipoprotein in plasma. `[PR]_T/[PR]_P` for a neutral
    /// molecule, whose binding protein Rodgers and Rowland take to be lipoprotein.
    pub lipoprotein_ratio: f64,
}

/// The lipid in plasma, as fractions of the plasma's volume. `f_NL,P` and `f_NP,P`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasmaLipids {
    /// Neutral lipid in plasma. `f_NL,P`.
    pub neutral_lipid: f64,
    /// Neutral phospholipid in plasma. `f_NP,P`.
    pub neutral_phospholipid: f64,
}

/// What a blood cell is, as fractions of its volume. A cell has no extracellular water, and
/// Rodgers and Rowland take it to hold no albumin or lipoprotein (2006, Appendix B), so these
/// three are all there is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BloodCellComposition {
    /// Water inside the cell. `f_IW,BC`.
    pub intracellular_water: f64,
    /// Neutral lipid. `f_NL,BC`.
    pub neutral_lipid: f64,
    /// Neutral phospholipid. `f_NP,BC`.
    pub neutral_phospholipid: f64,
}

/// Rat liver, from Rodgers and Rowland 2006, Table 1.
pub const RAT_LIVER: TissueComposition = TissueComposition {
    neutral_lipid: 0.0135,
    neutral_phospholipid: 0.0238,
    water: 0.642,
    lipoprotein_ratio: 0.161,
};

/// Rat kidney. Sources as [`RAT_LIVER`].
pub const RAT_KIDNEY: TissueComposition = TissueComposition {
    neutral_lipid: 0.0121,
    neutral_phospholipid: 0.0240,
    water: 0.672,
    lipoprotein_ratio: 0.137,
};

/// Rat skeletal muscle. Sources as [`RAT_LIVER`].
pub const RAT_MUSCLE: TissueComposition = TissueComposition {
    neutral_lipid: 0.0100,
    neutral_phospholipid: 0.0072,
    water: 0.726,
    lipoprotein_ratio: 0.059,
};

/// Rat adipose, from Rodgers and Rowland 2006, Table 1, with its two lipid columns put back the
/// way the 2005 table prints them — see the [module documentation](self).
///
/// Rodgers and Rowland use a vegetable oil:water partition coefficient for this tissue in place
/// of octanol:water; [`unbound_tissue_partition`] takes whichever it is given.
pub const RAT_ADIPOSE: TissueComposition = TissueComposition {
    neutral_lipid: 0.853,
    neutral_phospholipid: 0.0016,
    water: 0.144,
    lipoprotein_ratio: 0.068,
};

/// Rat plasma's lipid: 0.0023 neutral lipid and 0.0013 neutral phospholipid, the values Rodgers
/// and Rowland 2006 state in their text (p. 1241).
pub const RAT_PLASMA_LIPIDS: PlasmaLipids = PlasmaLipids {
    neutral_lipid: 0.0023,
    neutral_phospholipid: 0.0013,
};

/// Rat blood cells, from Rodgers, Leahy and Rowland 2005, Table 1.
pub const RAT_BLOOD_CELLS: BloodCellComposition = BloodCellComposition {
    intracellular_water: 0.603,
    neutral_lipid: 0.0017,
    neutral_phospholipid: 0.0029,
};

/// What a volume with these lipid fractions holds over the water beside it, at partition
/// coefficient `p`: `P·f_NL + (0.3P + 0.7)·f_NP`.
fn lipid_uptake(neutral_lipid: f64, neutral_phospholipid: f64, p: f64) -> f64 {
    p * neutral_lipid + (0.3 * p + 0.7) * neutral_phospholipid
}

/// `Kpu`: the concentration in a tissue over the **unbound** concentration in plasma, for a
/// neutral molecule.
///
/// `p` is the partition coefficient into the tissue's neutral lipid — n-octanol:water for most
/// tissues, a vegetable oil:water value for adipose — as the antilog, and `fu` is the fraction
/// unbound in plasma. The same `p` is used for the plasma's lipid, as in Rodgers and Rowland's
/// equation, where one `P` stands in both places.
///
/// # Its domain, which it does not check
///
/// It is the equation for `p ≥ 0` and `fu` in `(0, 1]`, and validates neither, being a formula
/// with nothing to refuse in. One case is worth knowing because it is reachable with sensible
/// inputs: when the plasma's lipid alone would bind more than `1/fu − 1` says the plasma binds —
/// a very lipophilic molecule with a large `fu` — the protein term goes negative. Rodgers and
/// Rowland do not clamp it and neither does this; it is the equation saying its inputs disagree.
pub fn unbound_tissue_partition(
    tissue: &TissueComposition,
    plasma: &PlasmaLipids,
    p: f64,
    fu: f64,
) -> f64 {
    let dissolved =
        tissue.water + lipid_uptake(tissue.neutral_lipid, tissue.neutral_phospholipid, p);
    let bound_to_protein_in_plasma =
        1.0 / fu - 1.0 - lipid_uptake(plasma.neutral_lipid, plasma.neutral_phospholipid, p);
    dissolved + bound_to_protein_in_plasma * tissue.lipoprotein_ratio
}

/// `K_p`: the concentration in a tissue over the **total** concentration in plasma.
/// [`unbound_tissue_partition`] times `fu`.
///
/// Divide by [`blood_to_plasma`] for the tissue-to-blood coefficient a
/// [perfusion](crate::CompartmentModel::perfuse) in this crate takes.
pub fn tissue_partition(tissue: &TissueComposition, plasma: &PlasmaLipids, p: f64, fu: f64) -> f64 {
    fu * unbound_tissue_partition(tissue, plasma, p, fu)
}

/// `Kpu,BC`: the concentration in blood cells over the unbound concentration in plasma, for a
/// neutral molecule. `f_IW,BC + P·f_NL,BC + (0.3P + 0.7)·f_NP,BC` — the tissue equation with no
/// extracellular water and no binding protein (Rodgers and Rowland 2006, Eq. B2 with `X = Y = 1`).
pub fn unbound_blood_cell_partition(cells: &BloodCellComposition, p: f64) -> f64 {
    cells.intracellular_water + lipid_uptake(cells.neutral_lipid, cells.neutral_phospholipid, p)
}

/// `B:P`: the concentration in whole blood over the concentration in plasma.
///
/// ```text
/// B:P = (1 − H) + H · fu · Kpu,BC = 1 + H (K_rbc:p − 1)
/// ```
///
/// with `H` the haematocrit — Rodgers and Rowland 2006, Eq. B3, divided through by the plasma
/// concentration. Neither paper prints the haematocrit it used, so it is an argument.
pub fn blood_to_plasma(cells: &BloodCellComposition, p: f64, fu: f64, haematocrit: f64) -> f64 {
    let cells_to_plasma = fu * unbound_blood_cell_partition(cells, p);
    1.0 + haematocrit * (cells_to_plasma - 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tissue with every fraction different from every other, so a field read in place of
    /// another changes the answer. Not a real tissue.
    const T: TissueComposition = TissueComposition {
        neutral_lipid: 0.05,
        neutral_phospholipid: 0.02,
        water: 0.7,
        lipoprotein_ratio: 0.3,
    };
    const PLASMA: PlasmaLipids = PlasmaLipids {
        neutral_lipid: 0.004,
        neutral_phospholipid: 0.002,
    };
    const NO_PLASMA_LIPID: PlasmaLipids = PlasmaLipids {
        neutral_lipid: 0.0,
        neutral_phospholipid: 0.0,
    };

    fn close(got: f64, want: f64, rel: f64) -> bool {
        (got - want).abs() <= rel * want.abs()
    }

    /// `10^x`, through the kernel's `exp` rather than `powf`, which is not bit-identical across
    /// platforms.
    fn antilog(x: f64) -> f64 {
        pantometry_core::math::exp(x * std::f64::consts::LN_10)
    }

    /// **`P → 0`, unbound in plasma, no plasma lipid:** nothing dissolves in lipid, nothing binds,
    /// and a tissue holds its water and the 70% of its phospholipid that behaves as water.
    #[test]
    fn with_no_lipophilicity_and_no_binding_only_the_water_is_left() {
        let got = unbound_tissue_partition(&T, &NO_PLASMA_LIPID, 0.0, 1.0);
        assert!(close(got, 0.7 + 0.7 * 0.02, 1e-15), "{got}");
    }

    /// **`Kpu` is a straight line in `P`**, and its slope is the tissue's lipid less its share of
    /// the plasma's: `f_NL + 0.3 f_NP − LR (f_NL,P + 0.3 f_NP,P)`. Read from the function across
    /// four decades of `P`, so a term with the wrong power of `P` bends the line.
    #[test]
    fn kpu_is_a_straight_line_in_p_with_the_lipid_slope() {
        let fu = 0.4;
        let at = |p: f64| unbound_tissue_partition(&T, &PLASMA, p, fu);
        let slope = 0.05 + 0.3 * 0.02 - 0.3 * (0.004 + 0.3 * 0.002);
        for (p0, p1) in [(0.0, 10.0), (10.0, 1e3), (1e3, 1e4)] {
            let measured = (at(p1) - at(p0)) / (p1 - p0);
            assert!(
                close(measured, slope, 1e-12),
                "{p0}..{p1}: {measured} against {slope}"
            );
        }
    }

    /// **Plasma partitions into itself at one.** A "tissue" that is plasma — all water, as the
    /// equation treats plasma's, its lipid the plasma's lipid, its lipoprotein the plasma's — has
    /// the total concentration plasma has, whatever the drug. This is what holds the plasma's
    /// lipid to being *subtracted* from the binding: added, or left out, it is counted twice or
    /// not at all, and the answer is not one.
    #[test]
    fn plasma_partitions_into_itself_at_one() {
        let plasma_as_tissue = TissueComposition {
            neutral_lipid: PLASMA.neutral_lipid,
            neutral_phospholipid: PLASMA.neutral_phospholipid,
            water: 1.0,
            lipoprotein_ratio: 1.0,
        };
        for p in [0.0, 1.0, 135.0, 1e4] {
            for fu in [1.0, 0.4, 0.05] {
                let got = tissue_partition(&plasma_as_tissue, &PLASMA, p, fu);
                // Terms up to `P·f_NL = 40` cancel to leave a number near `1/fu`.
                assert!(close(got, 1.0, 1e-13), "P {p}, fu {fu}: {got}");
            }
        }
    }

    /// **The protein term is what plasma binding leaves over**: `Kpu` rises with the tissue's
    /// lipoprotein by `1/fu − 1 − (P f_NL,P + (0.3P + 0.7) f_NP,P)` per unit of ratio, and at
    /// `fu = 1` with no plasma lipid nothing is bound and the ratio does not matter.
    #[test]
    fn the_protein_term_is_what_plasma_binding_leaves_over() {
        let (p, fu) = (50.0, 0.25);
        let none = TissueComposition {
            lipoprotein_ratio: 0.0,
            ..T
        };
        let per_ratio = (unbound_tissue_partition(&T, &PLASMA, p, fu)
            - unbound_tissue_partition(&none, &PLASMA, p, fu))
            / T.lipoprotein_ratio;
        let want = 1.0 / fu - 1.0 - (p * 0.004 + (0.3 * p + 0.7) * 0.002);
        assert!(close(per_ratio, want, 1e-13), "{per_ratio} against {want}");

        let unbound = unbound_tissue_partition(&T, &NO_PLASMA_LIPID, p, 1.0);
        let unbound_none = unbound_tissue_partition(&none, &NO_PLASMA_LIPID, p, 1.0);
        assert_eq!(unbound, unbound_none);
    }

    /// **`B:P` runs from plasma to cells**: one at a haematocrit of zero, the cells' own
    /// `fu·Kpu,BC` at a haematocrit of one, a straight line between. And a cell that is all water,
    /// in a plasma that binds nothing, is plasma at any haematocrit.
    #[test]
    fn blood_to_plasma_runs_from_plasma_to_cells() {
        let cells = BloodCellComposition {
            intracellular_water: 0.6,
            neutral_lipid: 0.002,
            neutral_phospholipid: 0.003,
        };
        let (p, fu) = (135.0, 0.4);
        let all_cells = fu * (0.6 + 135.0 * 0.002 + (0.3 * 135.0 + 0.7) * 0.003);
        assert_eq!(blood_to_plasma(&cells, p, fu, 0.0), 1.0);
        assert!(close(blood_to_plasma(&cells, p, fu, 1.0), all_cells, 1e-15));
        let half = blood_to_plasma(&cells, p, fu, 0.5);
        assert!(close(half, 0.5 * (1.0 + all_cells), 1e-15), "{half}");

        let water = BloodCellComposition {
            intracellular_water: 1.0,
            neutral_lipid: 0.0,
            neutral_phospholipid: 0.0,
        };
        for h in [0.0, 0.45, 1.0] {
            assert_eq!(blood_to_plasma(&water, p, 1.0, h), 1.0);
        }
    }

    /// **The compositions are volumes**, so no tissue's water and lipid exceed the whole of it,
    /// and adipose is mostly neutral lipid. The second is the transposition the 2006 table prints;
    /// the first is what mixing the 2006 lipids with the 2005 water gave adipose, 1.0066 of itself.
    #[test]
    fn every_tissue_fits_in_itself_and_adipose_is_fat() {
        for (name, t) in [
            ("liver", RAT_LIVER),
            ("kidney", RAT_KIDNEY),
            ("muscle", RAT_MUSCLE),
            ("adipose", RAT_ADIPOSE),
        ] {
            let filled = t.water + t.neutral_lipid + t.neutral_phospholipid;
            assert!(filled <= 1.0, "{name} is {filled} of itself");
            if name == "adipose" {
                assert!(t.neutral_lipid > 100.0 * t.neutral_phospholipid, "{t:?}");
            }
        }
    }

    /// **The authors' own predictions, reproduced.** Rodgers and Rowland 2006, Table 5, prints the
    /// `Kpu` this equation predicts with experimental log P (footnote b) for the four neutral
    /// drugs of their Table 2. These are all fourteen of those values in the four tissues this
    /// module carries.
    ///
    /// Adipose takes a vegetable oil:water coefficient, `log P_vo:w = 1.115 log P_o:w − 1.35`
    /// (Poulin and Theil 2000, the 2005 paper's Reference 1). That paper's Table 2 prints both
    /// coefficients for its compounds, and this gives seven of them to the printed digit.
    ///
    /// **The tolerance is the printed precision, worst case.** Each published value is good to
    /// half its last digit, and so is each input it was computed from — log P to 0.05, `fu` to
    /// half its last digit, each composition fraction to half its. Each value's bound is its output
    /// half-digit plus the change each input's half-digit makes, added. It was first written with
    /// log P, `fu` and the output only, and ftorafur's adipose and muscle then missed it by 0.0001
    /// and 0.0002. The composition's digits are inputs too, and with them the worst value uses 94%
    /// of its bound. That is a tolerance widened after a failure, and this is the record of it;
    /// the width is the published digits, not the failure.
    ///
    /// **Low in thirteen of the fourteen.** Only ethoxybenzamide's adipose comes out above the
    /// printed value, by 2.0%; cyclosporine is 3.4–4.7% under, ftorafur 0.4–4.3%, digoxin
    /// 1.1–2.5%, inside the digits but on one side of them. Something in the authors' arithmetic is not in their paper — a
    /// log P or `fu` carried to more digits than printed is the obvious candidate — and this
    /// test cannot say which.
    ///
    /// What it separates is a structural mistake from rounding: the 2005 water puts ftorafur's
    /// liver at 0.806 against 0.72, dropping the protein term at 0.669, and a transposed adipose
    /// puts ethoxybenzamide's at 0.878 against 0.48.
    #[test]
    fn reproduces_rodgers_and_rowlands_own_predictions() {
        let vegetable_oil = |log_p: f64| 1.115 * log_p - 1.35;
        // (drug, log P o:w, fu, half of fu's last digit, tissue, Kpu as printed, half its last
        // digit). Digoxin has no adipose value and ethoxybenzamide no liver.
        let cases: [(&str, f64, f64, f64, &str, f64, f64); 14] = [
            ("cyclosporine", 2.9, 0.083, 0.0005, "adipose", 69.4, 0.05),
            ("cyclosporine", 2.9, 0.083, 0.0005, "kidney", 17.9, 0.05),
            ("cyclosporine", 2.9, 0.083, 0.0005, "liver", 19.2, 0.05),
            ("cyclosporine", 2.9, 0.083, 0.0005, "muscle", 11.3, 0.05),
            ("digoxin", 1.2, 0.61, 0.005, "kidney", 1.10, 0.005),
            ("digoxin", 1.2, 0.61, 0.005, "liver", 1.11, 0.005),
            ("digoxin", 1.2, 0.61, 0.005, "muscle", 0.97, 0.005),
            ("ethoxybenzamide", 0.8, 0.59, 0.005, "adipose", 0.48, 0.005),
            ("ethoxybenzamide", 0.8, 0.59, 0.005, "kidney", 0.91, 0.005),
            ("ethoxybenzamide", 0.8, 0.59, 0.005, "muscle", 0.85, 0.005),
            ("ftorafur", -0.3, 0.78, 0.005, "adipose", 0.19, 0.005),
            ("ftorafur", -0.3, 0.78, 0.005, "kidney", 0.74, 0.005),
            ("ftorafur", -0.3, 0.78, 0.005, "liver", 0.72, 0.005),
            ("ftorafur", -0.3, 0.78, 0.005, "muscle", 0.76, 0.005),
        ];
        for (drug, log_p, fu, half_fu, tissue, printed, half_printed) in cases {
            let (t, half_lipid) = match tissue {
                "liver" => (RAT_LIVER, 0.00005),
                "kidney" => (RAT_KIDNEY, 0.00005),
                "muscle" => (RAT_MUSCLE, 0.00005),
                // 0.853 is printed to three decimals; every other lipid fraction to four.
                "adipose" => (RAT_ADIPOSE, 0.0005),
                _ => unreachable!(),
            };
            let p_of = |log_p: f64| {
                antilog(if tissue == "adipose" {
                    vegetable_oil(log_p)
                } else {
                    log_p
                })
            };
            let kpu = |t: &TissueComposition, log_p: f64, fu: f64| {
                unbound_tissue_partition(t, &RAT_PLASMA_LIPIDS, p_of(log_p), fu)
            };
            let got = kpu(&t, log_p, fu);
            let moved = |other: TissueComposition| (kpu(&other, log_p, fu) - got).abs();
            let bound = half_printed
                + (kpu(&t, log_p + 0.05, fu) - got).abs()
                + (kpu(&t, log_p, fu - half_fu) - got).abs()
                + moved(TissueComposition {
                    neutral_lipid: t.neutral_lipid + half_lipid,
                    ..t
                })
                + moved(TissueComposition {
                    neutral_phospholipid: t.neutral_phospholipid + 0.00005,
                    ..t
                })
                + moved(TissueComposition {
                    water: t.water + 0.0005,
                    ..t
                })
                + moved(TissueComposition {
                    lipoprotein_ratio: t.lipoprotein_ratio + 0.0005,
                    ..t
                });
            assert!(
                (got - printed).abs() <= bound,
                "{drug} in {tissue}: {got:.4} against the printed {printed}, bound {bound:.4}"
            );
        }
    }
}
