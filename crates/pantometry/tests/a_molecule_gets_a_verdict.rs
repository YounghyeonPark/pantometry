//! `pantometry::candidate`, held to what can be counted by hand and to the closed forms it is
//! built from.
//!
//! - **The descriptors** against a reading of each dictionary entry's atom list. The molecular
//!   weights are also checked against the figure a chemist would print.
//! - **Each call at its boundary**, with values chosen so the interval's ends are exact in binary.
//!   A limit is then met exactly, not approximately.
//! - **`K_d`** against `exp(ΔG°/RT)` with `R` written out separately, at a `ΔG°` chosen so that
//!   `K_d` is 1 µM. **`K_p`** at log P = 0, where `P` is exactly one and the equation can be summed
//!   by hand.
//! - **The verdict** as the worst of its two calls, in every combination that can be reached.

use pantometry::candidate::{
    assess, AssessError, BindingFreeEnergy, Call, Descriptors, Distribution, Evidence, Measured,
    Partition, Rule, Thresholds, Tissue,
};
use pantometry::forcefield::Component;
use pantometry::pharmacokinetic::partition::{
    tissue_partition, RAT_KIDNEY, RAT_LIVER, RAT_MUSCLE, RAT_PLASMA_LIPIDS,
};
use pantometry::units::Temperature;

const AIN: &str = include_str!("../../pantometry-forcefield/components/AIN.cif");
const ALA: &str = include_str!("../../pantometry-forcefield/components/ALA.cif");
const ARG: &str = include_str!("../../pantometry-forcefield/components/ARG.cif");
const ASN: &str = include_str!("../../pantometry-forcefield/components/ASN.cif");
const CYS: &str = include_str!("../../pantometry-forcefield/components/CYS.cif");
const MET: &str = include_str!("../../pantometry-forcefield/components/MET.cif");
const BNZ: &str = include_str!("../../pantometry-forcefield/components/BNZ.cif");
const MOH: &str = include_str!("../../pantometry-forcefield/components/MOH.cif");
const PXY: &str = include_str!("../../pantometry-forcefield/components/PXY.cif");

fn component(text: &str) -> Component {
    Component::from_ccd(text).unwrap()
}

/// The two usual thresholds, chosen by the caller, which is this test.
const USUAL: Thresholds = Thresholds {
    violations_for_alert: 2,
    ligand_efficiency: 0.3,
};

fn with_log_p(value: f64, sigma: f64) -> Evidence {
    Evidence {
        log_p: Some(Measured::new(value, sigma)),
        ..Evidence::default()
    }
}

fn binding(kcal: f64, sigma: f64) -> Option<BindingFreeEnergy> {
    Some(BindingFreeEnergy {
        kcal_per_mol: Measured::new(kcal, sigma),
        temperature: Temperature::kelvin(298.15),
    })
}

/// The three partitions, or a panic naming what came instead.
fn reported(d: Distribution) -> [Partition; 3] {
    match d {
        Distribution::Reported(p) => p,
        other => panic!("no partitions: {other:?}"),
    }
}

/// `(text, id, donors, acceptors, heavy atoms, [C, H, N, O, S])`.
type Case = (&'static str, &'static str, usize, usize, usize, [usize; 5]);

#[test]
fn the_descriptors_are_the_counts_of_the_atom_lists() {
    // Hand counts from each entry's `_chem_comp_atom` and `_chem_comp_bond` loops:
    // aspirin C9H8O4, one O-H (the carboxyl O1-HO1) and four oxygens;
    // alanine C3H7NO2 as the dictionary writes it, uncharged, with H and H2 on N and HXT on OXT:
    //   three donors, so an NH2 counts twice; one N and two O;
    // arginine C6H15N4O2 (+1): H, H2, HE, HH11, HH12, HH21, HH22, HXT, eight donors;
    // asparagine: H, H2, HD21, HD22, HXT, five donors, exactly Lipinski's limit;
    // cysteine C3H7NO2S: H, H2, HXT, three donors, and its S-H (SG-HG) is not one; three acceptors;
    // methionine C5H11NO2S: the same three and three, its sulfur a thioether;
    // benzene and p-xylene: nothing polar.
    let cases: [Case; 8] = [
        // (text, id, donors, acceptors, heavy, [C, H, N, O, S])
        (AIN, "AIN", 1, 4, 13, [9, 8, 0, 4, 0]),
        (ALA, "ALA", 3, 3, 6, [3, 7, 1, 2, 0]),
        (ARG, "ARG", 8, 6, 12, [6, 15, 4, 2, 0]),
        (ASN, "ASN", 5, 5, 9, [4, 8, 2, 3, 0]),
        (CYS, "CYS", 3, 3, 7, [3, 7, 1, 2, 1]),
        (MET, "MET", 3, 3, 9, [5, 11, 1, 2, 1]),
        (BNZ, "BNZ", 0, 0, 6, [6, 6, 0, 0, 0]),
        (PXY, "PXY", 0, 0, 8, [8, 10, 0, 0, 0]),
    ];
    for (text, id, donors, acceptors, heavy, [c, h, n, o, sulfur]) in cases {
        let comp = component(text);
        assert_eq!(comp.id(), id);
        let d = Descriptors::of(&comp);
        assert_eq!(d.donors, donors, "{id} donors");
        assert_eq!(d.acceptors, acceptors, "{id} acceptors");
        assert_eq!(d.heavy_atoms, heavy, "{id} heavy atoms");
        // The weight from the composition counted here, not from the atoms the function walked.
        let want = c as f64 * 12.011
            + h as f64 * 1.0080
            + n as f64 * 14.007
            + o as f64 * 15.999
            + sulfur as f64 * 32.06;
        assert!(
            (d.molecular_weight - want).abs() <= 1e-12 * want,
            "{id}: {} against {want}",
            d.molecular_weight
        );
        // And the rule of five is fed these, each to its own limit.
        let limits = assess(&comp, &Evidence::default(), &USUAL)
            .unwrap()
            .absorption
            .limits;
        let fed = limits.map(|l| l.value);
        assert_eq!(
            fed,
            [
                Some(Measured::new(donors as f64, 0.0)),
                Some(Measured::new(acceptors as f64, 0.0)),
                Some(Measured::new(d.molecular_weight, 0.0)),
                None,
            ],
            "{id}"
        );
    }
    // And against the printed figures, which round to two places.
    for (text, printed) in [(AIN, 180.16), (BNZ, 78.11), (PXY, 106.17)] {
        let mw = Descriptors::of(&component(text)).molecular_weight;
        assert!((mw - printed).abs() <= 0.005, "{mw} against {printed}");
    }
}

/// Benzene with H1 made a fluorine and H2 a chlorine. No shipped component has a halogen, and a
/// halogen is neither donor nor acceptor in Lipinski's counts.
#[test]
fn a_halogen_is_heavy_and_not_polar() {
    let text = BNZ
        .replace("BNZ H1 H1 H 0", "BNZ H1 H1 F 0")
        .replace("BNZ H2 H2 H 0", "BNZ H2 H2 CL 0")
        .replace(r#""C6 H6""#, r#""C6 H4 F CL""#);
    assert_eq!(
        text.matches(" F 0").count(),
        1,
        "the fluorine line was not written"
    );
    assert_eq!(
        text.matches(" CL 0").count(),
        1,
        "the chlorine line was not written"
    );
    let d = Descriptors::of(&component(&text));
    assert_eq!((d.donors, d.acceptors, d.heavy_atoms), (0, 0, 8));
    let want = 6.0 * 12.011 + 4.0 * 1.0080 + 18.998 + 35.45;
    assert!((d.molecular_weight - want).abs() <= 1e-12 * want);
    // Its formula was rewritten with it, so it is judged.
    let a = assess(&component(&text), &Evidence::default(), &USUAL).unwrap();
    assert_eq!(
        a.absorption.limit(Rule::Acceptors).value,
        Some(Measured::new(0.0, 0.0))
    );
}

/// Methanol's file cut down to two hydrogens bonded to each other: a molecule with no heavy atom.
fn dihydrogen() -> String {
    let mut out = String::new();
    for line in MOH.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let atom_row = w.len() == 21 && w[0] == "MOH";
        let bond_row = w.len() == 7 && w[0] == "MOH";
        if bond_row || (atom_row && !["H1", "H2"].contains(&w[1])) {
            continue;
        }
        out.push_str(line);
        out.push('\n');
        if line.trim() == "_chem_comp_bond.pdbx_ordinal" {
            out.push_str("MOH H1 H2 SING N N 1\n");
        }
    }
    assert_eq!(
        out.matches(r#""C H4 O""#).count(),
        1,
        "methanol's formula line moved"
    );
    out.replace(r#""C H4 O""#, "H2")
}

#[test]
fn a_molecule_with_no_heavy_atom_has_no_ligand_efficiency() {
    let h2 = component(&dihydrogen());
    assert_eq!(h2.atoms().len(), 2);
    assert_eq!(h2.bonds().len(), 1);
    assert_eq!(Descriptors::of(&h2).heavy_atoms, 0);
    // Judged without a binding free energy, it is judged.
    assert!(assess(&h2, &Evidence::default(), &USUAL).is_ok());
    let evidence = Evidence {
        binding: binding(-1.0, 0.1),
        ..Evidence::default()
    };
    assert_eq!(
        assess(&h2, &evidence, &USUAL).unwrap_err(),
        AssessError::NoHeavyAtoms
    );
}

#[test]
fn a_donor_is_found_whichever_way_round_the_bond_is_written() {
    let forward = Descriptors::of(&component(MOH));
    assert_eq!(forward.donors, 1);
    let reversed_text = MOH.replace("MOH O HO SING", "MOH HO O SING");
    assert_ne!(reversed_text, MOH, "the bond line to reverse was not found");
    let reversed = Descriptors::of(&component(&reversed_text));
    assert_eq!(reversed.donors, 1);
}

#[test]
fn each_limit_is_called_at_its_boundary() {
    let benzene = component(BNZ);
    let alert_at_one = Thresholds {
        violations_for_alert: 1,
        ..USUAL
    };
    let log_p_call = |value: f64, sigma: f64| {
        assess(&benzene, &with_log_p(value, sigma), &alert_at_one)
            .unwrap()
            .absorption
            .limit(Rule::LogP)
            .call
    };
    // 4.5 + 2 x 0.25 is 5 exactly: at the limit passes, as "greater than 5" fails.
    assert_eq!(log_p_call(4.5, 0.25), Call::Pass);
    assert_eq!(log_p_call(4.5, 0.25 + 1e-12), Call::Undecided);
    // 5.5 - 2 x 0.25 is 5 exactly: the interval touches the limit, so it is not wholly beyond.
    assert_eq!(log_p_call(5.5, 0.25), Call::Undecided);
    assert_eq!(log_p_call(5.5, 0.25 - 1e-12), Call::Fail);
    assert_eq!(log_p_call(5.0, 0.0), Call::Pass);
    assert_eq!(
        log_p_call(f64::from_bits(5f64.to_bits() + 1), 0.0),
        Call::Fail
    );

    // Asparagine has five donors, which is the limit, and arginine eight.
    let donors = |text: &str| {
        assess(&component(text), &Evidence::default(), &USUAL)
            .unwrap()
            .absorption
            .limit(Rule::Donors)
            .to_owned()
    };
    assert_eq!(donors(ASN).call, Call::Pass);
    assert_eq!(donors(ARG).call, Call::Fail);
    let names = assess(&benzene, &Evidence::default(), &USUAL)
        .unwrap()
        .absorption
        .limits
        .map(|l| (l.rule, l.rule.to_string(), l.rule.most()));
    assert_eq!(
        names,
        [
            (Rule::Donors, "H-bond donors".to_string(), 5.0),
            (Rule::Acceptors, "H-bond acceptors".to_string(), 10.0),
            (Rule::MolecularWeight, "molecular weight".to_string(), 500.0),
            (Rule::LogP, "log P".to_string(), 5.0),
        ]
    );
    assert_eq!(names.map(|n| n.0), Rule::ALL);
}

#[test]
fn the_alert_counts_exceeded_and_open_limits() {
    let benzene = component(BNZ);
    let arginine = component(ARG);
    let call = |c: &Component, evidence: Evidence, alert: usize| {
        let t = Thresholds {
            violations_for_alert: alert,
            ..USUAL
        };
        let a = assess(c, &evidence, &t).unwrap().absorption;
        (a.exceeded, a.open, a.call)
    };
    let none = Evidence::default;
    // Benzene, nothing exceeded. An unmeasured log P can only matter if one limit is enough.
    assert_eq!(call(&benzene, none(), 2), (0, 1, Call::Pass));
    assert_eq!(call(&benzene, none(), 1), (0, 1, Call::NotMeasured));
    assert_eq!(
        call(&benzene, with_log_p(5.0, 0.5), 1),
        (0, 1, Call::Undecided)
    );
    assert_eq!(call(&benzene, with_log_p(6.0, 0.0), 1), (1, 0, Call::Fail));
    assert_eq!(call(&benzene, with_log_p(6.0, 0.0), 2), (1, 0, Call::Pass));
    // Arginine, one limit certainly exceeded: the second decides at two.
    assert_eq!(call(&arginine, none(), 1), (1, 1, Call::Fail));
    assert_eq!(call(&arginine, none(), 2), (1, 1, Call::NotMeasured));
    assert_eq!(
        call(&arginine, with_log_p(5.0, 0.5), 2),
        (1, 1, Call::Undecided)
    );
    assert_eq!(
        call(&arginine, with_log_p(-4.0, 0.1), 2),
        (1, 0, Call::Pass)
    );
    assert_eq!(call(&arginine, with_log_p(6.0, 0.1), 2), (2, 0, Call::Fail));
    assert_eq!(call(&arginine, with_log_p(6.0, 0.1), 3), (2, 0, Call::Pass));
}

#[test]
fn the_dissociation_constant_is_the_closed_form() {
    // Two temperatures, so the one given is the one used.
    for t in [298.15, 310.15] {
        dissociation_constant_at(t);
    }
}

fn dissociation_constant_at(t: f64) {
    // R = N_A k_B, both exact since 2019, over 4184 J per kcal.
    let r = 6.022_140_76e23 * 1.380_649e-23 / 4184.0;
    let rt = r * t;
    // The ΔG° at which K_d is 1 µM.
    let dg = rt * (1e-6f64).ln();
    let evidence = Evidence {
        binding: Some(BindingFreeEnergy {
            kcal_per_mol: Measured::new(dg, 0.25),
            temperature: Temperature::kelvin(t),
        }),
        ..Evidence::default()
    };
    let a = assess(&component(BNZ), &evidence, &USUAL)
        .unwrap()
        .affinity
        .unwrap();
    assert!(
        (a.dissociation_constant - 1e-6).abs() <= 1e-13,
        "{}",
        a.dissociation_constant
    );
    let [tight, loose] = a.dissociation_interval;
    let factor = (0.5 / rt).exp();
    assert!((tight * factor / 1e-6 - 1.0).abs() <= 1e-12, "{tight}");
    assert!((loose / factor / 1e-6 - 1.0).abs() <= 1e-12, "{loose}");
    // Benzene has six heavy atoms.
    assert!((a.ligand_efficiency.value - (-dg / 6.0)).abs() <= 1e-15);
    assert_eq!(a.ligand_efficiency.sigma, 0.25 / 6.0);
}

#[test]
fn ligand_efficiency_is_called_at_its_boundary() {
    // Benzene's six heavy atoms: -3 ± 0.75 kcal/mol is 0.5 ± 0.125, an interval [0.25, 0.75],
    // and 0 ± 0.75 is [-0.25, 0.25]. Every one of those is exact in binary.
    let benzene = component(BNZ);
    let call = |dg: f64, threshold: f64| {
        let evidence = Evidence {
            binding: binding(dg, 0.75),
            ..Evidence::default()
        };
        let t = Thresholds {
            violations_for_alert: 2,
            ligand_efficiency: threshold,
        };
        assess(&benzene, &evidence, &t)
            .unwrap()
            .affinity
            .unwrap()
            .call
    };
    assert_eq!(call(-3.0, 0.25), Call::Pass);
    assert_eq!(call(-3.0, 0.25 + 1e-12), Call::Undecided);
    assert_eq!(call(-3.0, 0.75), Call::Undecided);
    assert_eq!(call(0.0, 0.25), Call::Undecided);
    assert_eq!(call(0.0, 0.25 + 1e-12), Call::Fail);
}

#[test]
fn the_partition_at_log_p_zero_is_the_sum_by_hand() {
    // P = 1 exactly, and with fu = 1 the bound term is the plasma lipid's uptake alone:
    // Kp = f_W + f_NL + f_NP - (f_NL,P + f_NP,P) x ratio.
    let evidence = Evidence {
        log_p: Some(Measured::new(0.0, 0.0)),
        fraction_unbound: Some(1.0),
        ..Evidence::default()
    };
    let a = assess(&component(BNZ), &evidence, &USUAL).unwrap();
    let partitions = reported(a.distribution);
    assert_eq!(partitions.map(|p| p.tissue), Tissue::ALL);
    assert_eq!(
        Tissue::ALL.map(|t| t.to_string()),
        ["liver", "kidney", "muscle"]
    );
    for (p, t) in partitions.iter().zip([RAT_LIVER, RAT_KIDNEY, RAT_MUSCLE]) {
        assert_eq!(p.tissue.composition(), t);
        assert_eq!(a.distribution.partition(p.tissue), Some(p));
        let by_hand = t.water + t.neutral_lipid + t.neutral_phospholipid
            - (0.0023 + 0.0013) * t.lipoprotein_ratio;
        assert!((p.kp - by_hand).abs() <= 1e-15, "{}: {}", p.tissue, p.kp);
        assert_eq!(p.kp_interval, [p.kp, p.kp]);
    }
    let liver = 0.642 + 0.0135 + 0.0238 - 0.0036 * 0.161;
    assert!((partitions[0].kp - liver).abs() <= 1e-15);
}

#[test]
fn the_partition_interval_spans_log_p_plus_and_minus_two_sigma() {
    let (fu, lp) = (0.25, Measured::new(2.0, 0.5));
    let evidence = Evidence {
        log_p: Some(lp),
        fraction_unbound: Some(fu),
        ..Evidence::default()
    };
    let a = assess(&component(PXY), &evidence, &USUAL).unwrap();
    for (p, t) in reported(a.distribution)
        .iter()
        .zip([RAT_LIVER, RAT_KIDNEY, RAT_MUSCLE])
    {
        let at = |x: f64| tissue_partition(&t, &RAT_PLASMA_LIPIDS, x, fu);
        // 10^1, 10^2, 10^3, against the kernel's exp of a multiple of ln 10.
        for (got, p) in [
            (p.kp_interval[0], 10.0),
            (p.kp, 100.0),
            (p.kp_interval[1], 1000.0),
        ] {
            assert!((got / at(p) - 1.0).abs() <= 1e-13, "P = {p}: {got}");
        }
    }
    // Not given either input, nothing is reported.
    let only_fu = Evidence {
        fraction_unbound: Some(fu),
        ..Evidence::default()
    };
    for evidence in [with_log_p(2.0, 0.5), Evidence::default(), only_fu] {
        let d = assess(&component(PXY), &evidence, &USUAL)
            .unwrap()
            .distribution;
        assert_eq!(d, Distribution::NotMeasured);
        assert_eq!(d.partition(Tissue::Liver), None);
    }
}

#[test]
fn the_verdict_is_the_worse_of_the_two_screens() {
    // Absorption is set by log P with the alert at one, binding by ΔG° on benzene's six atoms
    // against 0.25: -3 ± 0 passes, 0 ± 0 fails, -1.5 ± 0.75 is undecided.
    let benzene = component(BNZ);
    let t = Thresholds {
        violations_for_alert: 1,
        ligand_efficiency: 0.25,
    };
    let absorption = [
        (Some(Measured::new(1.0, 0.0)), Call::Pass),
        (Some(Measured::new(9.0, 0.0)), Call::Fail),
        (Some(Measured::new(5.0, 1.0)), Call::Undecided),
        (None, Call::NotMeasured),
    ];
    let affinity = [
        (binding(-3.0, 0.0), Call::Pass),
        (binding(0.0, 0.0), Call::Fail),
        (binding(-1.5, 0.75), Call::Undecided),
        (None, Call::NotMeasured),
    ];
    // The documented order, written out rather than read from `Ord`.
    let rank = |c: Call| match c {
        Call::Pass => 0,
        Call::NotMeasured => 1,
        Call::Undecided => 2,
        Call::Fail => 3,
    };
    let mut sorted = [Call::Fail, Call::Pass, Call::Undecided, Call::NotMeasured];
    sorted.sort();
    assert_eq!(sorted.map(rank), [0, 1, 2, 3]);
    let mut seen = 0;
    for (log_p, want_a) in absorption {
        for (b, want_b) in affinity {
            let evidence = Evidence {
                log_p,
                binding: b,
                ..Evidence::default()
            };
            let a = assess(&benzene, &evidence, &t).unwrap();
            assert_eq!(a.absorption.call, want_a);
            assert_eq!(a.affinity.map_or(Call::NotMeasured, |x| x.call), want_b);
            let want = if rank(want_a) >= rank(want_b) {
                want_a
            } else {
                want_b
            };
            assert_eq!(a.verdict, want, "{want_a:?} with {want_b:?}");
            seen += 1;
        }
    }
    assert_eq!(seen, 16);
}

#[test]
fn every_refusal_is_its_own_variant() {
    let benzene = component(BNZ);
    let refused = |evidence: Evidence, t: Thresholds| assess(&benzene, &evidence, &t).unwrap_err();
    assert!(matches!(
        refused(with_log_p(f64::NAN, 0.1), USUAL),
        AssessError::NotAMeasurement { what: "log P", .. }
    ));
    assert!(matches!(
        refused(with_log_p(1.0, -0.1), USUAL),
        AssessError::NotAMeasurement { what: "log P", .. }
    ));
    let bad_binding = Evidence {
        binding: binding(-5.0, f64::INFINITY),
        ..Evidence::default()
    };
    assert!(matches!(
        refused(bad_binding, USUAL),
        AssessError::NotAMeasurement {
            what: "the binding free energy",
            ..
        }
    ));
    let cold = Evidence {
        binding: Some(BindingFreeEnergy {
            kcal_per_mol: Measured::new(-5.0, 0.1),
            temperature: Temperature::kelvin(0.0),
        }),
        ..Evidence::default()
    };
    assert_eq!(refused(cold, USUAL), AssessError::Temperature(0.0));
    let hot = Evidence {
        binding: Some(BindingFreeEnergy {
            kcal_per_mol: Measured::new(-5.0, 0.1),
            temperature: Temperature::kelvin(f64::INFINITY),
        }),
        ..Evidence::default()
    };
    assert!(matches!(refused(hot, USUAL), AssessError::Temperature(_)));
    assert!(matches!(
        refused(with_log_p(f64::INFINITY, 0.1), USUAL),
        AssessError::NotAMeasurement { what: "log P", .. }
    ));
    for fu in [0.0, 1.5, f64::NAN] {
        let e = Evidence {
            fraction_unbound: Some(fu),
            ..Evidence::default()
        };
        assert!(matches!(refused(e, USUAL), AssessError::FractionUnbound(_)));
    }
    let one = Evidence {
        fraction_unbound: Some(1.0),
        ..Evidence::default()
    };
    assert!(assess(&benzene, &one, &USUAL).is_ok());
    for n in [0, 5] {
        let t = Thresholds {
            violations_for_alert: n,
            ..USUAL
        };
        assert_eq!(refused(Evidence::default(), t), AssessError::Violations(n));
    }
    for n in [1, 4] {
        let t = Thresholds {
            violations_for_alert: n,
            ..USUAL
        };
        assert!(assess(&benzene, &Evidence::default(), &t).is_ok());
    }
    for x in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let t = Thresholds {
            ligand_efficiency: x,
            ..USUAL
        };
        assert!(matches!(
            refused(Evidence::default(), t),
            AssessError::LigandEfficiency(_)
        ));
    }
}

/// Alanine with its hydrogens taken out of the atom and bond loops: lighter, no donors, and
/// refused, because its entry still says C3 H7 N O2.
#[test]
fn a_component_without_its_hydrogens_is_refused() {
    let stripped: String = ALA
        .lines()
        .filter(|line| {
            let w: Vec<&str> = line.split_whitespace().collect();
            let atom_row = w.len() == 21 && w[0] == "ALA" && w[3] == "H";
            let bond_row =
                w.len() == 7 && w[0] == "ALA" && (w[1].starts_with('H') || w[2].starts_with('H'));
            !(atom_row || bond_row)
        })
        .map(|l| format!("{l}\n"))
        .collect();
    let c = component(&stripped);
    assert_eq!(c.atoms().len(), 6);
    assert_eq!(Descriptors::of(&c).donors, 0);
    assert_eq!(
        assess(&c, &Evidence::default(), &USUAL).unwrap_err(),
        AssessError::Formula {
            stated: "C3 H7 N O2".into(),
            read: "C3 N O2".into(),
        }
    );
}

/// Arginine carries +1 as the dictionary writes it, so the neutral equation is not used on it.
#[test]
fn a_charged_molecule_gets_no_partition() {
    let evidence = Evidence {
        log_p: Some(Measured::new(-4.0, 0.1)),
        fraction_unbound: Some(0.5),
        ..Evidence::default()
    };
    let arg = component(ARG);
    assert_eq!(arg.formal_charge(), Some(1));
    let a = assess(&arg, &evidence, &USUAL).unwrap();
    assert_eq!(a.distribution, Distribution::Charged(1));
    // Neutral alanine, the same evidence: reported.
    let a = assess(&component(ALA), &evidence, &USUAL).unwrap();
    assert!(matches!(a.distribution, Distribution::Reported(_)));
}

/// Every K_d and K_p goes through the kernel's `exp`, so their bits are the same on every
/// platform, and these are those bits. A platform `exp` or `powf` in their place would move them.
#[test]
fn the_numbers_are_the_same_bits_everywhere() {
    let evidence = Evidence {
        log_p: Some(Measured::new(1.92, 0.23)),
        binding: binding(-4.61, 0.29),
        fraction_unbound: Some(1.0),
    };
    let a = assess(&component(BNZ), &evidence, &USUAL).unwrap();
    let kd = a.affinity.unwrap().dissociation_constant;
    let liver = reported(a.distribution)[0];
    let bits = [kd, liver.kp, liver.kp_interval[0], liver.kp_interval[1]].map(f64::to_bits);
    println!("{bits:#018x?}");
    assert_eq!(bits, PINNED);

    // Windows' own `exp` agrees with the kernel's for K_d at -4.61, so that pin alone would not
    // see one put in its place. At -2.06 kcal/mol it does not agree: the kernel gives
    // 3.0903631851106728e-2 and Windows' `exp` ...724e-2. Found by trying 2 000 values from -1 to
    // -11 kcal/mol through this module's RT, of which three differed: -2.06, -2.19 and -2.77.
    let kd = assess(
        &component(BNZ),
        &Evidence {
            binding: binding(-2.06, 0.0),
            ..Evidence::default()
        },
        &USUAL,
    )
    .unwrap()
    .affinity
    .unwrap()
    .dissociation_constant;
    assert_eq!(kd.to_bits(), 3.0903631851106728e-2f64.to_bits(), "{kd:e}");
}

const PINNED: [u64; 4] = [
    0x3f3b_5f8d_cbf0_eafe,
    0x4002_b6c9_4202_bd2a,
    0x3ff3_dc4f_e782_3ffa,
    0x4016_05f3_febf_5051,
];

#[test]
fn the_calls_print_as_words() {
    let words = [Call::Pass, Call::NotMeasured, Call::Undecided, Call::Fail].map(|c| c.to_string());
    assert_eq!(words, ["pass", "not measured", "undecided", "fail"]);
    assert_eq!(Measured::new(1.0, 0.0123).to_string(), "1.00 ± 0.01");
    assert_eq!(
        format!("{:.3}", Measured::new(1.0, 0.0123)),
        "1.000 ± 0.012"
    );
}

/// The same binding free energy in either unit is the same evidence. -4.61 kcal/mol is
/// -19.288 24 kJ/mol, and 0.29 is 1.213 36.
#[test]
fn a_binding_free_energy_in_kilojoules_is_the_same_one() {
    let t = Temperature::kelvin(298.15);
    let kcal = BindingFreeEnergy::from_kcal_per_mol(Measured::new(-4.61, 0.29), t);
    assert_eq!(
        kcal,
        BindingFreeEnergy {
            kcal_per_mol: Measured::new(-4.61, 0.29),
            temperature: t,
        }
    );
    let kj = BindingFreeEnergy::from_kj_per_mol(Measured::new(-19.28824, 1.21336), t);
    let close = |a: f64, b: f64| (a / b - 1.0).abs() <= 1e-15;
    assert!(close(kj.kcal_per_mol.value, -4.61), "{:?}", kj.kcal_per_mol);
    assert!(close(kj.kcal_per_mol.sigma, 0.29), "{:?}", kj.kcal_per_mol);
    assert_eq!(kj.temperature, t);
    let judge = |b| {
        assess(
            &component(BNZ),
            &Evidence {
                binding: Some(b),
                ..Evidence::default()
            },
            &USUAL,
        )
        .unwrap()
        .affinity
        .unwrap()
    };
    let (a, b) = (judge(kcal), judge(kj));
    assert!(close(a.dissociation_constant, b.dissociation_constant));
    assert!(close(a.ligand_efficiency.sigma, b.ligand_efficiency.sigma));
    assert_eq!(a.call, b.call);
}

/// Aspirin with nothing measured, written out by hand from its counts: C9H8O4 is 180.159, one
/// O-H, four oxygens, thirteen heavy atoms; log P open, so one open limit and a pass at two.
#[test]
fn an_assessment_prints_as_a_block() {
    let a = assess(&component(AIN), &Evidence::default(), &USUAL).unwrap();
    let want = concat!(
        "  molecular weight 180.16, donors 1, acceptors 4, heavy atoms 13\n",
        "    H-bond donors                  1  at most 5     pass\n",
        "    H-bond acceptors               4  at most 10    pass\n",
        "    molecular weight          180.16  at most 500   pass\n",
        "    log P               not measured  at most 5     not measured\n",
        "  absorption: pass (0 exceeded, 1 open)\n",
        "  binding: not measured\n",
        "  distribution: not measured\n",
        "  verdict: not measured\n",
    );
    assert_eq!(a.to_string(), want);

    // Each of the other branches, by the line it prints.
    let evidence = Evidence {
        log_p: Some(Measured::new(-4.0, 0.1)),
        binding: binding(-4.61, 0.29),
        fraction_unbound: Some(0.5),
    };
    let arg = assess(&component(ARG), &evidence, &USUAL)
        .unwrap()
        .to_string();
    assert!(
        arg.contains("    log P               -4.00 ± 0.10  at most 5     pass\n"),
        "{arg}"
    );
    assert!(
        arg.contains("  absorption: pass (1 exceeded, 0 open)\n"),
        "{arg}"
    );
    assert!(
        arg.contains("  distribution: net charge 1, not computed\n"),
        "{arg}"
    );
    // -4.61 ± 0.29 kcal/mol at 298.15 K is K_d 4.177e-4 M, 1.569e-4 to 1.112e-3 over ±2σ,
    // worked by hand by `physics-checker` for 360a89b, and tighter first.
    assert!(
        arg.contains("  binding: K_d 4.18e-4 M (1.57e-4 to 1.11e-3), "),
        "{arg}"
    );
    // Arginine's twelve heavy atoms: 4.61/12 = 0.38 ± 0.02, which passes 0.3.
    assert!(
        arg.contains("ligand efficiency 0.38 ± 0.02 kcal/mol per heavy atom: pass\n"),
        "{arg}"
    );
    assert!(arg.ends_with("  verdict: pass\n"), "{arg}");
    let ala = assess(&component(ALA), &evidence, &USUAL).unwrap();
    let text = ala.to_string();
    let kp = reported(ala.distribution);
    for (p, name) in kp.iter().zip(["liver ", "kidney", "muscle"]) {
        let line = format!(
            "  K_p {name} {:.3} ({:.3} to {:.3})\n",
            p.kp, p.kp_interval[0], p.kp_interval[1]
        );
        assert!(text.contains(&line), "{line:?} in {text}");
    }
    assert_eq!(text.lines().count(), 11, "{text}");
}
