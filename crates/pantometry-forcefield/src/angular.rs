//! The Universal Force Field's angle bend, torsion and inversion, with their analytic forces.
//!
//! From Rappé et al., *J. Am. Chem. Soc.* **114**, 10024 (1992); page and equation numbers are that
//! paper's, read off the scanned pages on 2026-10-01. [`crate::energy::ForceField`] builds one
//! [`Bend`] per pair of bonds sharing an atom, one [`Torsion`] per `I–J–K–L` path about every bond
//! the rules give a barrier, and three [`Inversion`]s per centre the rules give one.
//!
//! # Angle bend (eqs 8–13, pp. 10027–28)
//!
//! **General nonlinear centre** (eq 11, p. 10027): `E = K [C₀ + C₁ cos θ + C₂ cos 2θ]` with
//! (eq 12) `C₂ = 1 / (4 sin² θ₀)`, `C₁ = −4 C₂ cos θ₀`, `C₀ = C₂ (2 cos² θ₀ + 1)`. These make
//! `E(θ₀) = 0`, `E′(θ₀) = 0` and `E″(θ₀) = 4 K C₂ sin² θ₀ = K` — the paper's "second derivative at
//! θ₀ equal to the force constant" — which the tests check rather than assume.
//!
//! **Linear and trigonal-planar centres** (eq 10, p. 10027): `E = K/n² [1 − cos nθ]` with n = 3
//! for trigonal planar. Square planar and octahedral (n = 4) are not reached by any type here: no
//! type in [`UffType::ALL`] has θ₀ = 90°.
//!
//! **A correction: the linear form is `K (1 + cos θ)`, not eq 10's `K (1 − cos θ)`.** Eq 10 with
//! n = 1 is `K (1 − cos θ)`, which is zero at θ = 0 and largest, `2K`, at 180° — the opposite of
//! a linear centre, whose natural angle (Table I) is 180°. No reading of the printed form puts its
//! minimum at 180°, so it is taken as a misprint and the sign is the one the physics requires;
//! RDKit and Open Babel are reported to use the same `K (1 + cos θ)` (not checked for this
//! crate). The n = 3 form needs no
//! correction: `1 − cos 3θ` is zero at 120°. Both special forms, like the general one, have
//! curvature exactly `K` at their minimum: `d²/dθ² [K (1 + cos θ)] = K` at 180° and
//! `d²/dθ² [K/9 (1 − cos 3θ)] = K` at 120°.
//!
//! **Which centre takes which form is decided by the centre's θ₀ in Table I** ([`bend_form`]):
//! 180° → linear (`C_1`, `N_1`, `O_1`, and the terminal `H_` and halogens, which are never a
//! centre), exactly 120° → trigonal (`C_R`, `C_2`, `N_R`, `O_2`, `S_2`), anything else → general.
//! The paper says the special forms are for "linear, trigonal-planar, square-planar, and
//! octahedral coordination environments" (p. 10027), which is a statement about geometry, and
//! the table's θ₀ is that geometry. The reading matters for `N_2` (111.2°), `O_R` (110.0°) and
//! `S_R` (92.2°): p. 10028 says `O_R`'s and `N_2`'s angles were *fitted* — to methyl vinyl ether
//! and dimethyldiazene — and a fitted θ₀ only does anything in a form that reads it, so they take
//! the general form even though they are sp² types.
//!
//! **Force constant** (eq 13, p. 10028):
//! `K_IJK = β (Z_I* Z_K* / r_IK⁵) r_IJ r_JK [3 r_IJ r_JK (1 − cos² θ₀) − r_IK² cos θ₀]` with
//! `β = 664.12 / (r_IJ r_JK)` and `r_IK² = r_IJ² + r_JK² − 2 r_IJ r_JK cos θ₀`, in
//! kcal mol⁻¹ rad⁻² with lengths in Å and Z* in e. `r_IJ` and `r_JK` are the two bonds' natural
//! lengths (eq 2, with each bond's own order), θ₀ the centre's, Z* the two ends'.
//!
//! # Torsion (eqs 14–17, pp. 10028–29)
//!
//! `E = ½ V [1 − cos(n φ₀) cos(n φ)]` (eq 15), φ the I–J–K–L dihedral with 0 eclipsed. "All
//! torsions about this bond are considered, with each torsional barrier being divided by the
//! number of torsions present about this J–K bond" (p. 10028): every path I–J–K–L with I ≠ K,
//! L ≠ J and **I ≠ L** — the last excludes the non-torsion a three-membered ring would make — gets
//! `V / m`, m the number of such paths. The cases, decided by [`torsion_parameters`] from the two
//! central types' [`Hybridisation`]:
//!
//! | case | n | φ₀ | V | where |
//! | --- | --- | --- | --- | --- |
//! | sp³–sp³ | 3 | 180° | `√(V_J V_K)`, Table III (eq 16) | p. 10028 (a) |
//! | sp³–sp³, both group 6 (O, S) | 2 | 90° | `√(V_J V_K)`, V = 2 for O and 6.8 for S | p. 10028 |
//! | sp²–sp² | 2 | 180° | `5 √(U_J U_K) (1 + 4.18 ln BO)` (eq 17) | p. 10028 (c) |
//! | sp³ group 6 – sp² not group 6 | 2 | 90° | eq 17 | p. 10028 |
//! | sp³ – sp², the sp² bonded to another sp² | 3 | 180° | 2 kcal/mol | p. 10029, "propene" |
//! | sp³–sp², otherwise | 6 | 0° | 1 kcal/mol | p. 10028 (b) |
//! | either centre sp, or H or a halogen | — | — | 0 | p. 10029 |
//!
//! The rows are tried in that order, and the order is a decision: the paper describes the group-6
//! exception first and then calls propene "the remaining exception", so an sp³ oxygen on an sp²
//! carbon that is itself bonded to another sp² atom — an ester's ether oxygen on its carbonyl
//! carbon — takes the group-6 row. "Torsional potentials for central bonds involving
//! non-main-group elements were assigned a value of zero" (p. 10028) never applies: every element
//! here is main group. **"Another sp² atom" includes an `O_2`**, by the same reading of sp² as
//! [`UffType::hybridisation`] makes, so every carbonyl carbon satisfies it, and with this crate's
//! types the general row (b) is reached only by an sp² atom whose other neighbours are all sp,
//! sp³ or hydrogen — the `C_2` of a cumulene's end, the `N_2` of an isocyanate.
//!
//! **A consequence of the propene row, recorded rather than changed:** with n = 3 and a planar
//! sp² centre, the two torsions from each sp³ substituent to the sp² atom's two other neighbours
//! differ in φ by 180°, so their `cos 3φ` cancel and the bond's six terms sum to `V/2` whatever
//! the rotation. Propene's and acetaldehyde's methyl barriers therefore come from van der Waals
//! alone in UFF as written. Aspirin's acetyl C8–C9 reads exactly 1.0000 kcal/mol for that reason.
//!
//! **The one row to distrust** is "sp³ group 6 – sp²". The paper is explicit — "eq 17 is used
//! directly. For this case the periodicity (n) is 2 and the equilibrium angle φ₀ is 90°" — and
//! that is what is implemented; RDKit is believed to do the same (not checked for this crate). But with eq 15 it puts the
//! minimum of the C–O torsion of an ester, an anisole or a carboxylic acid at 90°, with a 10
//! kcal/mol barrier at planarity, and esters are planar. Measured on aspirin at the dictionary's
//! (planar) geometry: its three such bonds — the acid's O1–C7 and the ester's C2–O3 and O3–C8 —
//! carry 29.7 of its 30.7 kcal/mol of torsion energy, each within 0.25 kcal/mol of its maximum. The paper also says its periodicities
//! and minima "are the same as those described in the recently published DREIDING force field";
//! if DREIDING's rule for this case is planar, the 90° here is a misprint. DREIDING has not been
//! read for this crate, and minimising anisole against Table II is what will settle it.
//!
//! **Near a collinear centre the torsion is switched off, continuously.** "When angles about the
//! central atoms approach 180°, the potential energy and derivative terms are set to zero as is
//! conventionally done" (p. 10029). The dihedral is undefined when I–J–K or J–K–L is straight, so
//! something must be done; setting the term to zero below a threshold, as written, would put a
//! jump in the energy. Instead each torsion is multiplied by `s(θ_IJK) s(θ_JKL)`, where `s` is 1
//! for angles between 10° and 170° and falls to 0 at 180° (and, for the same reason, though the
//! paper does not mention it, at 0°) as the smoothstep `t² (3 − 2t)` in `t = (1 − |cos θ|) /
//! (1 − cos 10°)`. Both `s` and its derivative are continuous, and `s` vanishes like `sin⁴ θ`,
//! faster than the dihedral's gradient grows (`1 / sin θ`), so the force goes to zero smoothly
//! too. The 170° is this crate's choice — the paper gives no number — and it changes nothing at
//! any angle below it.
//!
//! # Inversion (eqs 18–19, p. 10029)
//!
//! For an atom I bonded to exactly three others, `E = K [C₀ + C₁ cos ω + C₂ cos 2ω]` (eq 18), ω
//! the angle between the I–L bond and the I–J–K plane, taken about each of the three bonds in turn
//! and **each divided by three**. The paper gives:
//!
//! - **`C_2` and `C_R`**: C₀ = 1, C₁ = −1, C₂ = 0, so `E = K (1 − cos ω)`, zero when planar; K = 6
//!   kcal/mol, or **50 when the carbon is bonded to an `O_2`** (fitted to formaldehyde's b₂ wag).
//! - **Nitrogen: zero** — "inversion terms corresponding to E_barrier = 0 for nitrogen". The
//!   sentence is about the group-5 hydrides' pyramidal inversion, but it names the element, not a
//!   type, and the paper gives sp² nitrogen nothing else; so no nitrogen — `N_3`, `N_2` or `N_R`,
//!   an amide's included — has an inversion term. An amide nitrogen is held planar here by its
//!   angles and its 1.41-order torsion, not by an inversion.
//! - **Phosphorus: a barrier of 22 kcal/mol**, with "the ω₀ … obtained from standard reference
//!   structures of the hydrides and the Cₙ fit to a minimum with E = 0 at ω₀ and that E for the
//!   maximum at ω = 0° be equal to E_barrier". Applied to `P_3+3` (a phosphine), whose θ₀ is PH₃'s
//!   93.8°. ω₀ is the angle between one bond and the plane of the other two in a pyramid with
//!   that bond angle: `sin ω₀ = √(1 − 3c² + 2c³) / sin θ₀`, c = cos θ₀, which is 84.4339°.
//!   With C₂ = 1, the three conditions give C₁ = −4 cos ω₀, C₀ = 2 cos² ω₀ + 1, and
//!   `K = 22 / (2 (1 − cos ω₀)²)` kcal/mol, so that the three terms at planarity total 22.
//!   A `P_3+5` with three neighbours gets none: the barrier is PH₃'s, which is trivalent.
//! - **Every other three-coordinate atom: none.** Oxygen (an oxonium) is not given a barrier by
//!   the paper; a three-coordinate sulfur is refused before this point.
//!
//! The printed relation between eq 18's ω and eq 19's γ, `ω = γ − π`, is inconsistent with eq 19
//! itself, whose `sin γ` stands where eq 18 has `cos ω`; `ω = π/2 − γ` is the relation that makes
//! the two agree. Nothing here uses eq 19.
//!
//! ω is signed — the side of the plane the bond is on — but the energy reads only `cos ω =
//! √(1 − sin² ω)` and `cos 2ω`. At ω = ±90°, a bond perpendicular to the plane of the other two,
//! `cos ω` has a kink and the force on the `C₁` term is unbounded; the derivative is evaluated
//! with `cos ω` no smaller than 10⁻⁸ there. No sp² carbon reaches it.
//!
//! # Units
//!
//! SI, as in [`crate::energy`]: angles in radians, energies in joules per molecule, `K_IJK` in
//! J rad⁻², forces in newtons.

use crate::ccd::{Element, ANGSTROM};
use crate::uff::{Hybridisation, UffType, KCAL_PER_MOL};

/// The numerator of eq 13's β = 664.12 / (r_IJ r_JK), in kcal mol⁻¹ Å e⁻² (p. 10028).
pub const ANGLE_FORCE_PREFACTOR: f64 = 664.12;

/// The 5 of eq 17 (p. 10028), kcal/mol per kcal/mol of `√(U_J U_K)`.
pub const SP2_BARRIER_SCALE: f64 = 5.0;

/// The 4.18 of eq 17 (p. 10028).
pub const SP2_BARRIER_BOND_ORDER: f64 = 4.18;

/// The sp²–sp³ barrier of case (b), p. 10028: 1 kcal/mol, in joules per molecule.
pub const SP2_SP3_BARRIER: f64 = 1.0 * KCAL_PER_MOL;

/// The propene exception's barrier, p. 10029: 2 kcal/mol, in joules per molecule.
pub const PROPENE_BARRIER: f64 = 2.0 * KCAL_PER_MOL;

/// `V_J` for an sp³ oxygen in the group-6 exception, p. 10028: 2 kcal/mol, in joules per molecule.
pub const GROUP6_OXYGEN_BARRIER: f64 = 2.0 * KCAL_PER_MOL;

/// `V_J` for the other sp³ group-6 elements (here, sulfur) in the group-6 exception, p. 10028:
/// 6.8 kcal/mol, in joules per molecule.
pub const GROUP6_OTHER_BARRIER: f64 = 6.8 * KCAL_PER_MOL;

/// The inversion force constant of a `C_2` or `C_R` carbon, p. 10029: 6 kcal/mol, in joules per
/// molecule.
pub const CARBON_INVERSION: f64 = 6.0 * KCAL_PER_MOL;

/// The inversion force constant of a `C_2` or `C_R` carbon bonded to an `O_2`, p. 10029: 50
/// kcal/mol, in joules per molecule.
pub const CARBONYL_INVERSION: f64 = 50.0 * KCAL_PER_MOL;

/// Phosphorus's inversion barrier, p. 10029: 22 kcal/mol, in joules per molecule.
pub const PHOSPHORUS_INVERSION_BARRIER: f64 = 22.0 * KCAL_PER_MOL;

/// The angle at a torsion's central atom past which the torsion is switched off, 170°, in
/// radians. This crate's choice; see the module documentation.
pub const TORSION_SWITCH_ANGLE: f64 = 170.0 * std::f64::consts::PI / 180.0;

/// The smallest `cos ω` an inversion's derivative is evaluated at. See the module documentation.
const MIN_COS_OMEGA: f64 = 1e-8;

// ---------------------------------------------------------------- angle bend

/// Which of the paper's angle functions a centre takes. See [`bend_form`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BendForm {
    /// Eq 11, `K [C₀ + C₁ cos θ + C₂ cos 2θ]`, with the coefficients of eq 12.
    General {
        /// C₀ = C₂ (2 cos² θ₀ + 1).
        c0: f64,
        /// C₁ = −4 C₂ cos θ₀.
        c1: f64,
        /// C₂ = 1 / (4 sin² θ₀).
        c2: f64,
    },
    /// `K (1 + cos θ)`: eq 10 with n = 1 and its sign corrected so the minimum is at 180°.
    Linear,
    /// Eq 10 with n = 3, `K/9 (1 − cos 3θ)`.
    Trigonal,
}

/// The angle function for a centre of type `j`: linear when its θ₀ is 180°, trigonal when it is
/// exactly 120°, the general form of eq 11 otherwise. See the module documentation for why.
pub fn bend_form(j: UffType) -> BendForm {
    let theta0 = j.table_i().theta0;
    if theta0 == 180.0 {
        BendForm::Linear
    } else if theta0 == 120.0 {
        BendForm::Trigonal
    } else {
        let t = j.parameters().angle;
        let (s, c) = t.sin_cos();
        let c2 = 1.0 / (4.0 * s * s);
        BendForm::General {
            c0: c2 * (2.0 * c * c + 1.0),
            c1: -4.0 * c2 * c,
            c2,
        }
    }
}

/// Eq 13, the force constant of an angle I–J–K with types `i`, `j`, `k` and bond natural lengths
/// `r_ij` and `r_jk` (metres), in J rad⁻². θ₀ is `j`'s; Z* are `i`'s and `k`'s.
pub fn bend_force_constant(i: UffType, j: UffType, k: UffType, r_ij: f64, r_jk: f64) -> f64 {
    let c = j.parameters().angle.cos();
    let (a, b) = (r_ij / ANGSTROM, r_jk / ANGSTROM);
    let r_ik2 = a * a + b * b - 2.0 * a * b * c;
    let r_ik = r_ik2.sqrt();
    let r_ik5 = r_ik2 * r_ik2 * r_ik;
    let beta = ANGLE_FORCE_PREFACTOR / (a * b);
    let kcal = beta * i.table_i().z1 * k.table_i().z1 / r_ik5
        * a
        * b
        * (3.0 * a * b * (1.0 - c * c) - r_ik2 * c);
    kcal * KCAL_PER_MOL
}

/// One angle-bend term, I–J–K with J the centre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bend {
    /// I, J, K as indices into the component; J is the centre.
    pub atoms: [usize; 3],
    /// The centre's θ₀, radians.
    pub natural_angle: f64,
    /// `K_IJK` of eq 13, J rad⁻².
    pub force_constant: f64,
    /// Which function of θ.
    pub form: BendForm,
}

impl Bend {
    /// The term for types `[i, j, k]` whose two bonds have natural lengths `r_ij` and `r_jk`
    /// (metres).
    pub fn new(atoms: [usize; 3], types: [UffType; 3], r_ij: f64, r_jk: f64) -> Bend {
        let [i, j, k] = types;
        Bend {
            atoms,
            natural_angle: j.parameters().angle,
            force_constant: bend_force_constant(i, j, k, r_ij, r_jk),
            form: bend_form(j),
        }
    }

    /// The energy at angle `theta` (radians), joules per molecule.
    pub fn energy(&self, theta: f64) -> f64 {
        self.at_cos(theta.cos()).0
    }

    /// The energy at positions `at` (metres), joules per molecule.
    pub fn energy_at(&self, at: &[[f64; 3]]) -> f64 {
        let [i, j, k] = self.atoms;
        let (u, v) = (sub(at[i], at[j]), sub(at[k], at[j]));
        self.at_cos(dot(u, v) / (norm(u) * norm(v))).0
    }

    /// The energy and `dE/d(cos θ)`: the one place each form is written. Every form is a
    /// polynomial in `cos θ`, so the derivative has no `1 / sin θ` in it.
    fn at_cos(&self, c: f64) -> (f64, f64) {
        let k = self.force_constant;
        match self.form {
            BendForm::General { c0, c1, c2 } => (
                k * (c0 + c1 * c + c2 * (2.0 * c * c - 1.0)),
                k * (c1 + 4.0 * c2 * c),
            ),
            BendForm::Linear => (k * (1.0 + c), k),
            // cos 3θ = 4c³ − 3c.
            BendForm::Trigonal => (
                k / 9.0 * (1.0 - (4.0 * c * c * c - 3.0 * c)),
                k / 9.0 * (3.0 - 12.0 * c * c),
            ),
        }
    }

    /// Adds this term's forces to `forces` and returns its energy.
    pub(crate) fn accumulate(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let [i, j, k] = self.atoms;
        let (u, v) = (sub(at[i], at[j]), sub(at[k], at[j]));
        let (lu, lv) = (norm(u), norm(v));
        let c = dot(u, v) / (lu * lv);
        let (e, de_dc) = self.at_cos(c);
        // ∂c/∂u = v/(|u||v|) − c u/|u|², and the same with u and v exchanged.
        let fi = scale(
            sub(scale(v, 1.0 / (lu * lv)), scale(u, c / (lu * lu))),
            -de_dc,
        );
        let fk = scale(
            sub(scale(u, 1.0 / (lu * lv)), scale(v, c / (lv * lv))),
            -de_dc,
        );
        add_to(&mut forces[i], fi);
        add_to(&mut forces[k], fk);
        add_to(&mut forces[j], scale(add(fi, fk), -1.0));
        e
    }
}

// ---------------------------------------------------------------- torsion

/// Which of the paper's torsion rules a central bond falls under. See [`torsion_parameters`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TorsionCase {
    /// Two sp³ centres: n = 3, φ₀ = 180°, eq 16 with Table III.
    Sp3Sp3,
    /// Two sp³ group-6 centres (O, S): n = 2, φ₀ = 90°, eq 16 with V = 2 (O) or 6.8.
    Group6Sp3Sp3,
    /// Two sp² centres: n = 2, φ₀ = 180°, eq 17.
    Sp2Sp2,
    /// An sp³ group-6 centre and an sp² centre of another group: n = 2, φ₀ = 90°, eq 17.
    Group6Sp3Sp2,
    /// An sp³ centre and an sp² centre bonded to another sp² atom: n = 3, φ₀ = 180°, 2 kcal/mol.
    Propene,
    /// Any other sp³–sp² pair: n = 6, φ₀ = 0°, 1 kcal/mol.
    Sp2Sp3,
}

/// What [`torsion_parameters`] decides for one central bond.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsionParameters {
    /// Which rule.
    pub case: TorsionCase,
    /// The bond's barrier V, joules per molecule, before division among its torsions.
    pub barrier: f64,
    /// The periodicity n.
    pub periodicity: u32,
    /// The equilibrium dihedral φ₀, radians.
    pub equilibrium: f64,
}

fn group6(t: UffType) -> bool {
    matches!(t.element(), Element::O | Element::S)
}

/// The torsion rule for a bond J–K of UFF bond order `order` between types `j` and `k`, where
/// `j_has_sp2_neighbour` says whether J has an sp² neighbour other than K, and likewise for K —
/// the propene exception's test. `None` when the paper gives the bond no torsion: either centre
/// sp, or a type with no hybridisation (hydrogen, a halogen). See the module documentation for
/// the table and the order the rules are tried in.
pub fn torsion_parameters(
    j: UffType,
    k: UffType,
    order: f64,
    j_has_sp2_neighbour: bool,
    k_has_sp2_neighbour: bool,
) -> Option<TorsionParameters> {
    use Hybridisation::{Sp2, Sp3};
    const DEG: f64 = std::f64::consts::PI / 180.0;
    let eq17 = || {
        let (uj, uk) = (
            j.sp2_torsion_barrier().expect("an sp2 type has a period"),
            k.sp2_torsion_barrier().expect("an sp2 type has a period"),
        );
        SP2_BARRIER_SCALE * (uj * uk).sqrt() * (1.0 + SP2_BARRIER_BOND_ORDER * order.ln())
    };
    let p = |case, barrier, periodicity, equilibrium_degrees: f64| {
        Some(TorsionParameters {
            case,
            barrier,
            periodicity,
            equilibrium: equilibrium_degrees * DEG,
        })
    };
    match (j.hybridisation()?, k.hybridisation()?) {
        (Sp3, Sp3) => {
            if group6(j) && group6(k) {
                let v = |t: UffType| {
                    if t.element() == Element::O {
                        GROUP6_OXYGEN_BARRIER
                    } else {
                        GROUP6_OTHER_BARRIER
                    }
                };
                p(TorsionCase::Group6Sp3Sp3, (v(j) * v(k)).sqrt(), 2, 90.0)
            } else {
                let (vj, vk) = (
                    j.sp3_torsion_barrier()
                        .expect("every sp3 type has a Table III value"),
                    k.sp3_torsion_barrier()
                        .expect("every sp3 type has a Table III value"),
                );
                p(TorsionCase::Sp3Sp3, (vj * vk).sqrt(), 3, 180.0)
            }
        }
        (Sp2, Sp2) => p(TorsionCase::Sp2Sp2, eq17(), 2, 180.0),
        (Sp3, Sp2) | (Sp2, Sp3) => {
            let (sp3, sp2, sp2_has_sp2) = if j.hybridisation() == Some(Sp3) {
                (j, k, k_has_sp2_neighbour)
            } else {
                (k, j, j_has_sp2_neighbour)
            };
            if group6(sp3) && !group6(sp2) {
                p(TorsionCase::Group6Sp3Sp2, eq17(), 2, 90.0)
            } else if sp2_has_sp2 {
                p(TorsionCase::Propene, PROPENE_BARRIER, 3, 180.0)
            } else {
                p(TorsionCase::Sp2Sp3, SP2_SP3_BARRIER, 6, 0.0)
            }
        }
        _ => None,
    }
}

/// One torsion term, I–J–K–L about the bond J–K.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torsion {
    /// I, J, K, L as indices into the component.
    pub atoms: [usize; 4],
    /// The central bond's rule, barrier, periodicity and φ₀.
    pub parameters: TorsionParameters,
    /// How many torsions share the central bond, and so divide its barrier.
    pub torsions_about_bond: usize,
}

impl Torsion {
    /// `V / m`, the barrier this one term carries, joules per molecule.
    pub fn share(&self) -> f64 {
        self.parameters.barrier / self.torsions_about_bond as f64
    }

    /// Eq 15 at dihedral `phi` (radians), divided among the bond's torsions, joules per molecule.
    /// **Without** the collinear switch, which depends on the two bond angles and not on φ; see
    /// [`Torsion::energy_at`] for the term as the force field evaluates it.
    pub fn energy(&self, phi: f64) -> f64 {
        self.at_cos(phi.cos()).0
    }

    /// The term at positions `at` (metres), switch included, joules per molecule.
    pub fn energy_at(&self, at: &[[f64; 3]]) -> f64 {
        let mut scratch = [[0.0; 3]; 4];
        let local: Vec<[f64; 3]> = self.atoms.iter().map(|&a| at[a]).collect();
        Torsion {
            atoms: [0, 1, 2, 3],
            ..*self
        }
        .accumulate(&local, &mut scratch)
    }

    /// Eq 15 and its derivative in `cos φ`, through `cos nφ = T_n(cos φ)`, the Chebyshev
    /// polynomial: the one place the formula is written.
    fn at_cos(&self, c: f64) -> (f64, f64) {
        let n = self.parameters.periodicity;
        let phase = (f64::from(n) * self.parameters.equilibrium).cos();
        let (t, dt) = chebyshev(n, c);
        let half = 0.5 * self.share();
        (half * (1.0 - phase * t), -half * phase * dt)
    }

    /// Adds this term's forces to `forces` and returns its energy, switch included.
    pub(crate) fn accumulate(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let [i, j, k, l] = self.atoms;
        let b1 = sub(at[j], at[i]);
        let b2 = sub(at[k], at[j]);
        let b3 = sub(at[l], at[k]);
        let a = cross(b1, b2);
        let b = cross(b2, b3);
        let (la, lb) = (norm(a), norm(b));
        if la == 0.0 || lb == 0.0 {
            return 0.0;
        }
        let (l1, l2, l3) = (norm(b1), norm(b2), norm(b3));
        let x1 = dot(b1, b2) / (l1 * l2);
        let x2 = dot(b2, b3) / (l2 * l3);
        let (s1, ds1) = switch(x1);
        let (s2, ds2) = switch(x2);
        if s1 == 0.0 || s2 == 0.0 {
            return 0.0;
        }
        let c = dot(a, b) / (la * lb);
        let (e, de_dc) = self.at_cos(c);
        let s = s1 * s2;
        // ∂E/∂a and ∂E/∂b of the unswitched term, then through a = b1 × b2 and b = b2 × b3.
        let ga = scale(
            sub(scale(b, 1.0 / (la * lb)), scale(a, c / (la * la))),
            s * de_dc,
        );
        let gb = scale(
            sub(scale(a, 1.0 / (la * lb)), scale(b, c / (lb * lb))),
            s * de_dc,
        );
        let mut g1 = cross(b2, ga);
        let mut g2 = add(cross(ga, b1), cross(b3, gb));
        let mut g3 = cross(gb, b2);
        // The switch: ∂x1/∂b1 = b2/(l1 l2) − x1 b1/l1², and likewise.
        let w1 = e * s2 * ds1;
        if w1 != 0.0 {
            g1 = add(
                g1,
                scale(
                    sub(scale(b2, 1.0 / (l1 * l2)), scale(b1, x1 / (l1 * l1))),
                    w1,
                ),
            );
            g2 = add(
                g2,
                scale(
                    sub(scale(b1, 1.0 / (l1 * l2)), scale(b2, x1 / (l2 * l2))),
                    w1,
                ),
            );
        }
        let w2 = e * s1 * ds2;
        if w2 != 0.0 {
            g2 = add(
                g2,
                scale(
                    sub(scale(b3, 1.0 / (l2 * l3)), scale(b2, x2 / (l2 * l2))),
                    w2,
                ),
            );
            g3 = add(
                g3,
                scale(
                    sub(scale(b2, 1.0 / (l2 * l3)), scale(b3, x2 / (l3 * l3))),
                    w2,
                ),
            );
        }
        // b1 = r_J − r_I, b2 = r_K − r_J, b3 = r_L − r_K; the force is minus the gradient.
        add_to(&mut forces[i], g1);
        add_to(&mut forces[j], sub(g2, g1));
        add_to(&mut forces[k], sub(g3, g2));
        add_to(&mut forces[l], scale(g3, -1.0));
        s * e
    }
}

/// The collinear switch `s` and `ds/dx` at `x = cos` of the angle between two consecutive bond
/// vectors (`|x| = 1` when the bend angle is 180° or 0°). See the module documentation.
pub fn switch(x: f64) -> (f64, f64) {
    let on = -TORSION_SWITCH_ANGLE.cos();
    let m = x.abs();
    if m <= on {
        return (1.0, 0.0);
    }
    if m >= 1.0 {
        return (0.0, 0.0);
    }
    let t = (1.0 - m) / (1.0 - on);
    let s = t * t * (3.0 - 2.0 * t);
    let ds_dm = -6.0 * t * (1.0 - t) / (1.0 - on);
    (s, ds_dm * x.signum())
}

/// `T_n(c)` and `T_n′(c) = n U_{n−1}(c)`, by the three-term recurrences.
fn chebyshev(n: u32, c: f64) -> (f64, f64) {
    let (mut t0, mut t1) = (1.0, c);
    let (mut u0, mut u1) = (0.0, 1.0); // U_{−1}, U_0
    if n == 0 {
        return (1.0, 0.0);
    }
    for _ in 1..n {
        (t0, t1) = (t1, 2.0 * c * t1 - t0);
        (u0, u1) = (u1, 2.0 * c * u1 - u0);
    }
    (t1, f64::from(n) * u1)
}

// ---------------------------------------------------------------- inversion

/// The inversion `K` (joules per molecule, before the division by three) and `[C₀, C₁, C₂]` for a
/// centre of type `centre` with exactly three neighbours, one of which is an `O_2` when
/// `bonded_to_o2`. `None` where the paper gives no barrier. See the module documentation.
pub fn inversion_parameters(centre: UffType, bonded_to_o2: bool) -> Option<(f64, [f64; 3])> {
    match centre {
        UffType::C2 | UffType::CR => Some((
            if bonded_to_o2 {
                CARBONYL_INVERSION
            } else {
                CARBON_INVERSION
            },
            [1.0, -1.0, 0.0],
        )),
        UffType::P3Trivalent => {
            let w0 = pyramid_inversion_angle(centre.parameters().angle);
            let c = w0.cos();
            let k = PHOSPHORUS_INVERSION_BARRIER / (2.0 * (1.0 - c) * (1.0 - c));
            Some((k, [2.0 * c * c + 1.0, -4.0 * c, 1.0]))
        }
        _ => None,
    }
}

/// ω₀ for a three-fold pyramid whose bonds meet at `theta` (radians): the angle between one bond
/// and the plane of the other two, `asin(√(1 − 3c² + 2c³) / sin θ)` with `c = cos θ`. The square
/// root is the volume of the parallelepiped of three unit bonds (a Gram determinant), and
/// `sin θ` the area of the face two of them span.
pub fn pyramid_inversion_angle(theta: f64) -> f64 {
    let (s, c) = theta.sin_cos();
    ((1.0 - 3.0 * c * c + 2.0 * c * c * c).sqrt() / s).asin()
}

/// One inversion term at centre I, about the axis I–L, against the plane I–J–K.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inversion {
    /// I (the centre), J, K, L (the axis) as indices into the component.
    pub atoms: [usize; 4],
    /// K of eq 18, joules per molecule, **before** the division by three.
    pub force_constant: f64,
    /// `[C₀, C₁, C₂]` of eq 18.
    pub coefficients: [f64; 3],
}

impl Inversion {
    /// Eq 18 at angle `omega` (radians) between the axis and the plane, divided by three.
    pub fn energy(&self, omega: f64) -> f64 {
        self.at_sin(omega.sin()).0
    }

    /// The term at positions `at` (metres), joules per molecule.
    pub fn energy_at(&self, at: &[[f64; 3]]) -> f64 {
        let mut scratch = [[0.0; 3]; 4];
        let local: Vec<[f64; 3]> = self.atoms.iter().map(|&a| at[a]).collect();
        Inversion {
            atoms: [0, 1, 2, 3],
            ..*self
        }
        .accumulate(&local, &mut scratch)
    }

    /// Eq 18 over three and its derivative in `sin ω`: the one place the formula is written.
    /// `cos ω = √(1 − sin² ω)` and `cos 2ω = 1 − 2 sin² ω`.
    fn at_sin(&self, s: f64) -> (f64, f64) {
        let [c0, c1, c2] = self.coefficients;
        let k = self.force_constant / 3.0;
        let cw = (1.0 - s * s).max(0.0).sqrt();
        let e = k * (c0 + c1 * cw + c2 * (1.0 - 2.0 * s * s));
        let de_ds = k * (-c1 * s / cw.max(MIN_COS_OMEGA) - 4.0 * c2 * s);
        (e, de_ds)
    }

    /// Adds this term's forces to `forces` and returns its energy. Zero, with no force, if J, I
    /// and K are collinear, when the plane is undefined.
    pub(crate) fn accumulate(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let [i, j, k, l] = self.atoms;
        let u = sub(at[j], at[i]);
        let v = sub(at[k], at[i]);
        let w = sub(at[l], at[i]);
        let n = cross(u, v);
        let (ln, lw) = (norm(n), norm(w));
        if ln == 0.0 {
            return 0.0;
        }
        let s = dot(n, w) / (ln * lw);
        let (e, de_ds) = self.at_sin(s);
        let uv = dot(u, v);
        // ∂s/∂u = (v × w)/(|n||w|) − s (|v|² u − (u·v) v)/|n|², and the same for v; and
        // ∂s/∂w = n/(|n||w|) − s w/|w|².
        let ds_du = sub(
            scale(cross(v, w), 1.0 / (ln * lw)),
            scale(sub(scale(u, dot(v, v)), scale(v, uv)), s / (ln * ln)),
        );
        let ds_dv = sub(
            scale(cross(w, u), 1.0 / (ln * lw)),
            scale(sub(scale(v, dot(u, u)), scale(u, uv)), s / (ln * ln)),
        );
        let ds_dw = sub(scale(n, 1.0 / (ln * lw)), scale(w, s / (lw * lw)));
        let (fj, fk, fl) = (
            scale(ds_du, -de_ds),
            scale(ds_dv, -de_ds),
            scale(ds_dw, -de_ds),
        );
        add_to(&mut forces[j], fj);
        add_to(&mut forces[k], fk);
        add_to(&mut forces[l], fl);
        add_to(&mut forces[i], scale(add(add(fj, fk), fl), -1.0));
        e
    }
}

// ---------------------------------------------------------------- vectors

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn add_to(f: &mut [f64; 3], d: [f64; 3]) {
    for k in 0..3 {
        f[k] += d[k];
    }
}
