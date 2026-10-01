//! **Each of the three energy terms against a number that does not come from this crate's code.**
//!
//! A natural length worked by hand from Table I, the paper's own worked electronegativity
//! correction, the Lennard-Jones well's depth and zero crossing, the Coulomb constant against
//! CODATA's ε₀, and the exclusion rule counted on chains and rings whose topology settles the
//! answer. Hand values were evaluated separately from the crate (the expressions are written next
//! to them) and typed here as literals.

use pantometry_forcefield::energy::{
    bond_order_correction, coulomb, electronegativity_correction, natural_length, Pair, Stretch,
    AMIDE_BOND_ORDER, AROMATIC_BOND_ORDER,
};
use pantometry_forcefield::uff::{self, gmp_electronegativity, ELEMENTARY_CHARGE, KCAL_PER_MOL};
use pantometry_forcefield::{Component, Element, ForceField, UffType};

const ANGSTROM: f64 = 1e-10;
const AIN: &str = include_str!("../components/AIN.cif");

/// A one-component CCD file: atoms as `(name, element, aromatic, [x, y, z] in Å)`, bonds as
/// `(first, second, order, aromatic)`.
fn entry(atoms: &[(&str, &str, bool, [f64; 3])], bonds: &[(&str, &str, &str, bool)]) -> String {
    let yn = |b: bool| if b { "Y" } else { "N" };
    let mut s = String::from(
        "data_TST\n_chem_comp.id TST\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (name, element, aromatic, [x, y, z]) in atoms {
        s += &format!("TST {name} {element} 0 {} {x} {y} {z}\n", yn(*aromatic));
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (a, b, order, aromatic) in bonds {
        s += &format!("{a} {b} {order} {}\n", yn(*aromatic));
    }
    s
}

fn force_field(c: &Component) -> ForceField {
    ForceField::new(c, &uff::assign(c)).expect("supported")
}

fn positions(c: &Component) -> Vec<[f64; 3]> {
    c.atoms().iter().map(|a| a.at).collect()
}

/// A chain of `n` carbons, single-bonded in order, no hydrogens: every atom `C_3`. A zigzag in
/// the xy plane, 1.5 Å bonds at 120°, so no two atoms coincide.
fn carbon_chain(n: usize) -> Component {
    let names: Vec<String> = (1..=n).map(|k| format!("C{k}")).collect();
    let atoms: Vec<(&str, &str, bool, [f64; 3])> = names
        .iter()
        .enumerate()
        .map(|(k, name)| {
            let y = if k % 2 == 0 { 0.0 } else { 0.75 };
            (name.as_str(), "C", false, [1.299 * k as f64, y, 0.0])
        })
        .collect();
    let bonds: Vec<(&str, &str, &str, bool)> = (1..n)
        .map(|k| (names[k - 1].as_str(), names[k].as_str(), "SING", false))
        .collect();
    Component::from_ccd(&entry(&atoms, &bonds)).expect("the chain parses")
}

/// A ring of `n` carbons, single-bonded, no hydrogens, on a circle of radius 2 Å.
fn carbon_ring(n: usize) -> Component {
    let names: Vec<String> = (1..=n).map(|k| format!("C{k}")).collect();
    let atoms: Vec<(&str, &str, bool, [f64; 3])> = names
        .iter()
        .enumerate()
        .map(|(k, name)| {
            let a = std::f64::consts::TAU * k as f64 / n as f64;
            let at = [
                (2000.0 * a.cos()).round() / 1000.0,
                (2000.0 * a.sin()).round() / 1000.0,
                0.0,
            ];
            (name.as_str(), "C", false, at)
        })
        .collect();
    let bonds: Vec<(&str, &str, &str, bool)> = (0..n)
        .map(|k| {
            (
                names[k].as_str(),
                names[(k + 1) % n].as_str(),
                "SING",
                false,
            )
        })
        .collect();
    Component::from_ccd(&entry(&atoms, &bonds)).expect("the ring parses")
}

// ---------------------------------------------------------------- bond stretch

/// **`C_3–C_3` single is 2 × 0.757 = 1.514 Å, exactly.** Equal electronegativities make `r_EN`
/// zero, and `ln 1 = 0` makes `r_BO` zero — both exactly, not to rounding. The natural length is
/// then `2 r₁` in Å times 10⁻¹⁰, against the literal 1.514e-10: two independent roundings of one
/// real number, so they agree to an ulp, 2.2e-16 relative.
#[test]
fn c3_c3_single_is_twice_the_radius() {
    assert_eq!(bond_order_correction(0.757, 0.757, 1.0), 0.0);
    assert_eq!(
        electronegativity_correction(0.757, 5.343, 0.757, 5.343),
        0.0
    );
    let r = natural_length(UffType::C3, UffType::C3, 1.0);
    assert!((r / 1.514e-10 - 1.0).abs() <= f64::EPSILON, "{r:e}");
}

/// **Between unlike atoms the electronegativity correction is subtracted** — the decision
/// `energy` documents. `C_3–O_3` single: `0.757 + 0.658 − r_EN` with
/// `r_EN = 0.757 × 0.658 × (√5.343 − √8.741)² / (5.343 × 0.757 + 8.741 × 0.658)` = 0.0211552 Å,
/// so 1.3938448 Å (with the printed `+` it would be 1.4361552). `C_3–N_3`: 1.4510710 Å (with `+`,
/// 1.4629290). These are dimethyl ether's and trimethylamine's natural lengths, the two the
/// decision rests on. A dropped `r_EN` gives 1.415 and 1.457 Å and fails here; nothing else in
/// this file would see it, because every other natural length is between equal electronegativities.
#[test]
fn a_heteroatom_bond_is_shortened_by_electronegativity() {
    let co = natural_length(UffType::C3, UffType::O3, 1.0) / ANGSTROM;
    assert!((co - 1.393_844_845_252_683_5).abs() < 1e-14, "{co}");
    let cn = natural_length(UffType::C3, UffType::N3, 1.0) / ANGSTROM;
    assert!((cn - 1.451_071_040_722_294).abs() < 1e-14, "{cn}");
}

/// **Eq 3 by hand for a double and an aromatic bond**, and the bond orders the dictionary's bonds
/// get. `C_R–C_R` aromatic: `−0.1332 × 1.458 × ln 1.5 = −0.0787435945992` Å. `C_2–O_2` double:
/// `−0.1332 × 1.366 × ln 2 = −0.1261189612795` Å. Agreement to 1e-15 Å is a few ulps of the
/// operands. Aspirin's ring bond C3–C4 is written `SING` with the aromatic flag and must get 1.5,
/// not 1; its natural length is then `1.458 − 0.0787436` Å, since `r_EN` is zero between carbons.
#[test]
fn bond_order_correction_by_hand() {
    let aromatic = bond_order_correction(0.729, 0.729, 1.5);
    assert!(
        (aromatic - -0.078_743_594_599_210_9).abs() < 1e-15,
        "{aromatic}"
    );
    let double = bond_order_correction(0.732, 0.634, 2.0);
    assert!(
        (double - -0.126_118_961_279_498_7).abs() < 1e-15,
        "{double}"
    );

    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    let name = |i: usize| c.atoms()[i].name.as_str();
    let find = |a: &str, b: &str| {
        ff.stretches()
            .iter()
            .find(|s| {
                let [i, j] = s.atoms;
                (name(i), name(j)) == (a, b) || (name(i), name(j)) == (b, a)
            })
            .unwrap_or_else(|| panic!("no {a}-{b} bond"))
    };
    let ring = find("C3", "C4");
    assert_eq!(ring.order, AROMATIC_BOND_ORDER);
    assert!((ring.natural_length / ANGSTROM - 1.379_256_405_400_789).abs() < 1e-14);
    // Kekulé DOUB inside the ring is 1.5 as well, not 2.
    assert_eq!(find("C3", "C2").order, AROMATIC_BOND_ORDER);
    assert_eq!(find("C7", "O2").order, 2.0);
    assert_eq!(find("C7", "O1").order, 1.0);
    assert_eq!(find("O3", "C8").order, 1.0);
}

/// **The paper's own worked example (p. 10027): the Si–`O_3_z` correction is 0.0533 Å.** With
/// χ_Si 4.168, χ_O 8.741, r_Si 1.117, r_O_3_z 0.528 — Si and `O_3_z` are not types this crate
/// ships, so their numbers live here — eq 4 gives 0.053255 Å. The paper prints three figures, so
/// the earned tolerance is half the last one: 0.00005 Å. And the oxygen χ used is the one the
/// crate ships.
#[test]
fn the_papers_silicon_oxygen_example() {
    let chi_o = gmp_electronegativity(Element::O);
    assert_eq!(chi_o, 8.741);
    let r_en = electronegativity_correction(1.117, 4.168, 0.528, chi_o);
    assert!((r_en - 0.0533).abs() <= 0.000_05, "{r_en}");
}

/// **`E = ½ k Δr²`, and exactly zero at the natural length.** `k` for `C_3–C_3` is eq 6 by hand:
/// `664.12 × 1.912² / 1.514³ = 699.5917987` kcal mol⁻¹ Å⁻². Through [`ForceField::evaluate`] on a
/// two-carbon molecule 1.614 Å long, the stretch is 0.1 Å and the energy `½ k (0.1)²` = 3.4980
/// kcal/mol. The stretch is a difference of two lengths of order 1.5e-10 m, each good to an ulp,
/// so it is good to about 3e-26 m in 1e-11 m — 3e-15 relative, doubled by squaring; 1e-13 is that
/// with room for the conversions.
#[test]
fn a_bond_is_harmonic_about_its_natural_length() {
    let s = Stretch::new([0, 1], UffType::C3, UffType::C3, 1.0);
    assert_eq!(s.energy(s.natural_length), 0.0);
    let k_kcal = s.force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
    assert!(
        (k_kcal / 699.591_798_712_679 - 1.0).abs() < 1e-13,
        "{k_kcal}"
    );

    let c = Component::from_ccd(&entry(
        &[
            ("C1", "C", false, [0.0, 0.0, 0.0]),
            ("C2", "C", false, [1.614, 0.0, 0.0]),
        ],
        &[("C1", "C2", "SING", false)],
    ))
    .expect("parses");
    let e = force_field(&c).energy(&positions(&c));
    let expected = 0.5 * 699.591_798_712_679 * 0.1 * 0.1;
    assert!(
        (e.bond / KCAL_PER_MOL / expected - 1.0).abs() < 1e-13,
        "{e:?}"
    );
    assert_eq!(e.van_der_waals, 0.0, "a bonded pair has no van der Waals");
}

/// **An amide C–N is order 1.41 (p. 10026), and an amine C–N is 1.** Acetamide and methylamine
/// differ at the nitrogen only by the neighbouring C=O.
#[test]
fn an_amide_bond_is_order_one_point_four_one() {
    let acetamide = Component::from_ccd(&entry(
        &[
            ("C1", "C", false, [0.0, 0.0, 0.0]),
            ("C2", "C", false, [1.5, 0.0, 0.0]),
            ("O", "O", false, [2.1, 1.0, 0.0]),
            ("N", "N", false, [2.2, -1.2, 0.0]),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C2", "O", "DOUB", false),
            ("C2", "N", "SING", false),
        ],
    ))
    .expect("parses");
    let ff = force_field(&acetamide);
    let orders: Vec<f64> = ff.stretches().iter().map(|s| s.order).collect();
    assert_eq!(orders, [1.0, 2.0, AMIDE_BOND_ORDER]);

    let methylamine = Component::from_ccd(&entry(
        &[
            ("C1", "C", false, [0.0, 0.0, 0.0]),
            ("N", "N", false, [1.47, 0.0, 0.0]),
        ],
        &[("C1", "N", "SING", false)],
    ))
    .expect("parses");
    assert_eq!(force_field(&methylamine).stretches()[0].order, 1.0);
}

/// **Natural length and force constant by hand for two heteronuclear bonds of non-integer order**
/// — where a wrong `Z*` pairing or a `k` computed from the wrong length would show, as it cannot in
/// `C_3–C_3`. Worked outside the crate from Table I, χ 5.343 (C) and 6.899 (N), λ 0.1332 and
/// eq 6's 664.12:
///
/// - `C_R–N_R`, n = 1.5 (pyridine's ring bond): `r_IJ = 1.428 − 0.1332·1.428·ln 1.5 − r_EN` with
///   `r_EN = 0.0058039` Å, so 1.3450728 Å; `k = 664.12 · 1.912 · 2.544 / r_IJ³` = 1327.43785
///   kcal mol⁻¹ Å⁻².
/// - `C_R–N_R`, n = 1.41 (an amide, now that the amide carbon is `C_R`): 1.3568421 Å and
///   1293.19388 kcal mol⁻¹ Å⁻². Checked directly and as acetamide's C–N bond through the typing.
///
/// 1e-13 relative is a few dozen roundings of the operands. `Z_C·Z_C` instead of `Z_C·Z_N` would
/// give 997.67 for the first; `k` from `r_I + r_J` instead of `r_IJ` would give 1109.35.
#[test]
fn heteronuclear_bonds_by_hand() {
    let close = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-13;
    let kcal_per_a2 = |k: f64| k * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
    for (n, r, k) in [
        (1.5, 1.345_072_784_173_839_2, 1_327.437_853_718_414_8),
        (1.41, 1.356_842_079_964_895, 1_293.193_884_725_723_8),
    ] {
        let s = Stretch::new([0, 1], UffType::CR, UffType::NR, n);
        assert!(close(s.natural_length / ANGSTROM, r), "n {n}: {s:?}");
        assert!(close(kcal_per_a2(s.force_constant), k), "n {n}: {s:?}");
        // Symmetric in the two types.
        let t = Stretch::new([0, 1], UffType::NR, UffType::CR, n);
        assert!(close(t.force_constant, s.force_constant), "n {n}");
    }

    let acetamide = Component::from_ccd(&entry(
        &[
            ("C1", "C", false, [0.0, 0.0, 0.0]),
            ("C2", "C", false, [1.5, 0.0, 0.0]),
            ("O", "O", false, [2.1, 1.0, 0.0]),
            ("N", "N", false, [2.2, -1.2, 0.0]),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C2", "O", "DOUB", false),
            ("C2", "N", "SING", false),
        ],
    ))
    .expect("parses");
    assert_eq!(
        uff::assign(&acetamide),
        [UffType::C3, UffType::CR, UffType::O2, UffType::NR]
    );
    let amide = force_field(&acetamide).stretches()[2];
    assert!(close(
        amide.natural_length / ANGSTROM,
        1.356_842_079_964_895
    ));
    assert!(close(
        kcal_per_a2(amide.force_constant),
        1_293.193_884_725_723_8
    ));
}

/// **An N-acyl aromatic nitrogen is not an amide** — a recorded choice (see `energy`): in
/// N-acetylimidazole the ring nitrogen's lone pair is in the aromatic sextet. Its exocyclic
/// N–C(=O) bond, written single and not aromatic, is order 1, not 1.41, and the acetyl carbon
/// stays `C_2`. The ring bonds are 1.5.
#[test]
fn an_acylated_aromatic_nitrogen_is_not_an_amide() {
    let c = Component::from_ccd(&entry(
        &[
            ("N1", "N", true, [0.0, 0.0, 0.0]),
            ("C2", "C", true, [1.1, 0.7, 0.0]),
            ("N3", "N", true, [2.1, 0.0, 0.0]),
            ("C4", "C", true, [1.7, -1.3, 0.0]),
            ("C5", "C", true, [0.4, -1.3, 0.0]),
            ("C6", "C", false, [-1.3, 0.5, 0.0]),
            ("O", "O", false, [-1.5, 1.7, 0.0]),
            ("C7", "C", false, [-2.4, -0.5, 0.0]),
        ],
        &[
            ("N1", "C2", "SING", true),
            ("C2", "N3", "DOUB", true),
            ("N3", "C4", "SING", true),
            ("C4", "C5", "DOUB", true),
            ("C5", "N1", "SING", true),
            ("N1", "C6", "SING", false),
            ("C6", "O", "DOUB", false),
            ("C6", "C7", "SING", false),
        ],
    ))
    .expect("parses");
    let types = uff::assign(&c);
    assert_eq!(types[0], UffType::NR);
    assert_eq!(types[5], UffType::C2, "the acetyl carbon is not resonant");
    let orders: Vec<f64> = force_field(&c)
        .stretches()
        .iter()
        .map(|s| s.order)
        .collect();
    assert_eq!(orders, [1.5, 1.5, 1.5, 1.5, 1.5, 1.0, 2.0, 1.0]);
}

/// **One natural length per element against `C_3`**, so that every χ the crate ships enters a
/// check. Single bonds, so `r_BO` is zero and each length is `0.757 + r_X − r_EN`, worked outside
/// the crate. **A transcription check, not a primary-source check**: the χ here are Open Babel's,
/// typed a second time, and the source paper has not been read. A digit transposed in χ_H (4.528 →
/// 4.258) moves C–H by 1.4e-3 Å, far outside 1e-13 relative.
#[test]
fn each_element_against_carbon() {
    for (t, r) in [
        (UffType::H, 1.109_400_794_877_744),
        (UffType::F, 1.381_519_561_409_680_7),
        (UffType::P3Trivalent, 1.857_944_793_344_279_6),
        (UffType::S3Divalent, 1.813_747_405_971_574_4),
        (UffType::Cl, 1.777_985_485_373_360_3),
        (UffType::Br, 1.933_432_295_730_878),
        (UffType::I, 2.131_992_569_490_815_3),
    ] {
        let got = natural_length(UffType::C3, t, 1.0) / ANGSTROM;
        assert!((got / r - 1.0).abs() < 1e-13, "C_3-{t}: {got} vs {r}");
    }
}

// ---------------------------------------------------------------- van der Waals

/// **The minimum is `−D_IJ`, exactly, at `x_IJ`.** There `x_IJ/x` is 1.0 exactly and so is every
/// power of it. And the well is a minimum: 0.1% either side is higher.
#[test]
fn the_well_is_minus_d_at_x() {
    let p = Pair::new([0, 1], UffType::C3, UffType::O2);
    assert_eq!(p.energy(p.distance), -p.well);
    assert!(p.energy(p.distance * 1.001) > -p.well);
    assert!(p.energy(p.distance * 0.999) > -p.well);
}

/// **It crosses zero at `x_IJ / 2^(1/6)`**, where `(x_IJ/x)⁶ = 2`. Not exactly in floating
/// point: `2^(1/6)`, the division and six multiplications each round, so `s⁶` is 2 to about 20 ε
/// and `E = D s⁶ (s⁶ − 2)` is about `2 D × 20 ε` — 1e-13 D is twice that.
#[test]
fn it_crosses_zero_at_x_over_the_sixth_root_of_two() {
    let p = Pair::new([0, 1], UffType::C3, UffType::C3);
    let zero = p.distance / 2f64.powf(1.0 / 6.0);
    assert!(p.energy(zero).abs() < 1e-13 * p.well, "{}", p.energy(zero));
    assert!(p.energy(zero * 0.99) > 0.0 && p.energy(zero * 1.01) < 0.0);
}

/// **Geometric combination, by hand for C–O** (eqs 21b and 22): `x = √(3.851 × 3.500) =
/// 3.671307669` Å and `D = √(0.105 × 0.060) = 0.0793725393` kcal/mol. The arithmetic means would
/// be 3.6755 Å and 0.0825 kcal/mol, 0.11% and 3.9% away — far outside 1e-14.
#[test]
fn the_combination_is_geometric() {
    let p = Pair::new([0, 1], UffType::C3, UffType::O3);
    assert!((p.distance / ANGSTROM / 3.671_307_668_937_595 - 1.0).abs() < 1e-14);
    assert!((p.well / KCAL_PER_MOL / 0.079_372_539_331_937_72 - 1.0).abs() < 1e-14);
}

// ---------------------------------------------------------------- electrostatics

/// **Two unit charges at 1 Å are 332.0637 kcal/mol**, repulsive; opposite charges are the same
/// magnitude, attractive. And the paper's constant is `e²/(4π ε₀)`: CODATA 2022's ε₀ gives
/// 332.063713 in these units, so the seven printed figures agree to half their last unit,
/// 1.5e-7 relative.
#[test]
fn two_unit_charges_at_one_angstrom() {
    let e = coulomb(1.0, 1.0, ANGSTROM) / KCAL_PER_MOL;
    assert!((e / 332.0637 - 1.0).abs() < 1e-15, "{e}");
    assert!(coulomb(1.0, -1.0, ANGSTROM) < 0.0);
    assert_eq!(coulomb(1.0, -1.0, ANGSTROM), -coulomb(1.0, 1.0, ANGSTROM));

    let epsilon_0 = 8.854_187_818_8e-12;
    let codata = ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        / (4.0 * std::f64::consts::PI * epsilon_0)
        / KCAL_PER_MOL
        / ANGSTROM;
    assert!((codata / 332.0637 - 1.0).abs() < 1.5e-7, "{codata}");
}

/// **No charges, no electrostatics — exactly zero**, the default; and one charged atom alone has
/// nothing to interact with.
#[test]
fn zero_charges_give_exactly_zero() {
    let c = carbon_chain(4);
    let ff = force_field(&c);
    assert_eq!(ff.energy(&positions(&c)).electrostatic, 0.0);
    let one = ff.with_charges(vec![1.0, 0.0, 0.0, 0.0]);
    assert_eq!(one.energy(&positions(&c)).electrostatic, 0.0);
}

// ---------------------------------------------------------------- exclusions

/// **In a three-atom chain every pair is 1-2 or 1-3, so nothing interacts non-bonded**, however
/// charged.
#[test]
fn a_three_atom_chain_has_no_nonbonded_pair() {
    let c = carbon_chain(3);
    let ff = force_field(&c).with_charges(vec![1.0, -1.0, 1.0]);
    assert!(ff.pairs().is_empty());
    let e = ff.energy(&positions(&c));
    assert_eq!((e.van_der_waals, e.electrostatic), (0.0, 0.0));
}

/// **In a four-atom chain only the 1-4 pair interacts — by count and by value.** Every atom is
/// charged, so a 1-2 or 1-3 pair let through would add its own Coulomb term and miss the hand
/// value. The 1-4 distance in the zigzag is `√(3.897² + 0.75²)` Å; the van der Waals energy is
/// `0.105 {(3.851/r)¹² − 2 (3.851/r)⁶}` and the Coulomb `332.0637 × 0.5 × −0.25 / r`, typed
/// here from the paper, not from the crate. Agreement to 1e-12 relative is the coordinates'
/// three-decimal parse and a few dozen roundings.
#[test]
fn a_four_atom_chain_has_only_its_one_four_pair() {
    let c = carbon_chain(4);
    let ff = force_field(&c).with_charges(vec![0.5, 1.0, -1.0, -0.25]);
    let pairs: Vec<[usize; 2]> = ff.pairs().iter().map(|p| p.atoms).collect();
    assert_eq!(pairs, [[0, 3]]);
    let e = ff.energy(&positions(&c));
    let r: f64 = (3.897f64 * 3.897 + 0.75 * 0.75).sqrt();
    let s6 = (3.851 / r).powi(6);
    let vdw = 0.105 * (s6 * s6 - 2.0 * s6);
    let q = 332.0637 * 0.5 * -0.25 / r;
    assert!((e.van_der_waals / KCAL_PER_MOL / vdw - 1.0).abs() < 1e-12);
    assert!((e.electrostatic / KCAL_PER_MOL / q - 1.0).abs() < 1e-12);
}

/// **Rings: a pair is excluded if either path is short.** In a five-membered ring every pair is
/// 1-2 or 1-3 one way round; in a six-membered ring only the three para pairs are 1-4 both ways.
#[test]
fn ring_pairs_are_counted_by_their_shortest_path() {
    assert!(force_field(&carbon_ring(4)).pairs().is_empty());
    assert!(force_field(&carbon_ring(5)).pairs().is_empty());
    let six: Vec<[usize; 2]> = force_field(&carbon_ring(6))
        .pairs()
        .iter()
        .map(|p| p.atoms)
        .collect();
    assert_eq!(six, [[0, 3], [1, 4], [2, 5]]);
}
