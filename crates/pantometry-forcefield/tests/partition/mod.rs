//! **What benzene's partition legs share**: the machinery 4c-2 wrote for its octanol and water legs
//! (`benzene_partitioned_into_octanol.rs`), moved here unchanged so that step C-2's decane leg
//! (`benzene_partitioned_into_decane.rs`) runs through the same code and reads the water leg back
//! through the same headers. Benzene's GAFF/AM1-BCC model, the two 4c-2 boxes, the melts, the
//! digests of the start and of the model, the resumable and headed windows, the lock, the log, the
//! schedule's refinement and the report.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`) with `tests/gaff/`, whose numbers
//! it uses.

#![allow(dead_code)]
// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

use crate::gaff::{self, Model, ONE_FOUR};
use crate::hydration::{self, Benzene, Solvated};
use crate::octanol::{
    self, production, Liquid, Octanol, Solute, ANGSTROM, DENSITY, FS, KELVIN, OC9,
};
use pantometry_forcefield::free_energy::{bennett_chain, window_seed};
use pantometry_forcefield::uff::{AVOGADRO, KCAL_PER_MOL};
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Component, Coupling, Estimate, Lambda, Mixing, MolecularDynamics,
    PeriodicDecoupling, PeriodicEnergy, PeriodicForceField, Protocol,
};
use pantometry_units::BOLTZMANN;
use std::ops::Range;

pub const EPS: f64 = f64::EPSILON;

/// The small octanol boxes: one site along the molecule's axis, four by four across.
pub const SMALL: [usize; 3] = [1, 4, 4];
/// The small boxes' cutoff, under half their shortest edge.
pub const SMALL_CUTOFF: f64 = 6.5 * ANGSTROM;

/// Every lattice octanol with an atom nearer than this to a benzene atom is removed.
pub const CLEARANCE: f64 = 2.0 * ANGSTROM;

/// **Benzene's volume in the solvent, m³ a molecule**: liquid benzene's molar volume, 89 cm³/mol
/// (quoted, as `tests/hydration/` quotes it for W3's five waters), standing in for its partial
/// molar volume in octanol. Its uncertainty, a few cm³/mol, moves the octanol's density by less
/// than 0.01% in the 192-site box.
pub const BENZENE_VOLUME: f64 = 89e-6 / AVOGADRO;

/// Experiment's `log K_ow` for benzene: see the module documentation.
pub const LOG_P: f64 = 2.13;

/// Variant 4's staged water leg, GAFF/AM1-BCC benzene's hydration free energy, kcal/mol
/// (`variant_4_staged_short_measured`, run 2026-10-20): W3's −1.570 plus the staged swap.
pub const HYDRATION_STAGED: (f64, f64) = (-0.82, 0.15);

pub fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The magnitudes the energy is a sum of, for its rounding.
pub fn magnitude(e: &PeriodicEnergy) -> f64 {
    let w = &e.ewald;
    e.bond.abs()
        + e.angle.abs()
        + e.torsion.abs()
        + e.inversion.abs()
        + e.van_der_waals.abs()
        + e.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs()
        + w.background.abs()
}

pub fn largest_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// A central difference and the tolerance it earns: the truncation by the difference at `2h`, and
/// the energy's rounding, `64 ε` of its parts over `h`. W1's.
pub fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    (d1, 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h)
}

// ---------------------------------------------------------------------------------------------
// The model and the boxes

/// Benzene with FreeSolv's AM1-BCC charges, relaxed in UFF's vacuum with them and put on its
/// constraints.
pub fn benzene() -> Benzene {
    let c = Component::from_ccd(hydration::BNZ).expect("BNZ");
    Benzene::new(Some(Model::benzene(&c).charges))
}

/// Octanol with FreeSolv's AM1-BCC charges, relaxed and in its frame.
pub fn octanol() -> Octanol {
    let c = Component::from_ccd(OC9).expect("OC9");
    Octanol::new(Some(Model::octanol(&c).charges))
}

/// **Benzene first and then `octanols` octanols, every atom GAFF's**, by FreeSolv's types and
/// charges, combined by Lorentz–Berthelot and with GAFF's 1-4 scales.
pub fn octanol_rules(
    field: PeriodicForceField,
    benzene: &Component,
    octanol: &Component,
    octanols: usize,
) -> PeriodicForceField {
    let (b, o) = (Model::benzene(benzene), Model::octanol(octanol));
    let mut vdw = b.vdw.clone();
    let mut q = b.charges.clone();
    for _ in 0..octanols {
        vdw.extend_from_slice(&o.vdw);
        q.extend_from_slice(&o.charges);
    }
    let n = vdw.len();
    field
        .with_nonbonded(0..n, vdw)
        .expect("GAFF's values are valid")
        .with_charges(q)
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1])
}

/// What a leg samples: its field (GAFF's, with the production's mesh and list where it has them),
/// the decoupling, the dynamics, the positions as built, each molecule's atoms and the masses (for
/// the pressure), how many solvent molecules, and how it is melted.
pub struct System {
    pub name: &'static str,
    pub field: PeriodicForceField,
    pub decoupling: PeriodicDecoupling,
    pub template: MolecularDynamics,
    pub built: Vec<[f64; 3]>,
    pub molecules: Vec<Range<usize>>,
    pub masses: Vec<f64>,
    pub solvent: usize,
    pub melt: &'static [MeltStage],
}

/// **Benzene in octanol**: [`Liquid::around`] with benzene at the box's centre and every octanol
/// within [`CLEARANCE`] of it removed, GAFF's everywhere. With `compensate`, **the lattice built at
/// the density that leaves the octanols at [`DENSITY`] in the volume benzene leaves them**,
/// `N_sites m / (N_oct m / ρ + v_B)`, by a fixed point in the number removed, which the density
/// moves; refused if it does not settle in ten rounds. The lattice density is returned too.
pub fn octanol_box(
    b: &Benzene,
    o: &Octanol,
    counts: [usize; 3],
    cutoff: f64,
    accuracy: f64,
    seed: u64,
    compensate: bool,
) -> (Liquid, PeriodicForceField, f64) {
    let sites = counts[0] * counts[1] * counts[2];
    let build = |density: f64| {
        let (cell, _, _) = Liquid::cell_of(o, counts, density);
        let centre = cell.lengths().map(|x| 0.5 * x);
        Liquid::around(
            Solute {
                component: &b.component,
                types: &b.types,
                charges: &b.charges,
                at: b.centred_at(centre),
                clearance: CLEARANCE,
            },
            o,
            counts,
            density,
            seed,
            cutoff,
            accuracy,
        )
    };
    let mut density = DENSITY;
    let mut l = build(density);
    if compensate {
        let mut settled = false;
        for _ in 0..10 {
            let want =
                sites as f64 * o.mass() / (l.octanols as f64 * o.mass() / DENSITY + BENZENE_VOLUME);
            if want.to_bits() == density.to_bits() {
                settled = true;
                break;
            }
            density = want;
            l = build(density);
        }
        assert!(settled, "the number removed did not settle");
    }
    let field = octanol_rules(l.field.clone(), &b.component, &o.component, l.octanols);
    (l, field, density)
}

/// [`System`] for benzene in octanol: `l` and its GAFF field `field`, the decoupling of benzene's
/// twelve atoms.
pub fn octanol_system(l: &Liquid, field: PeriodicForceField, melt: &'static [MeltStage]) -> System {
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < l.solute).collect();
    System {
        name: "octanol",
        decoupling: PeriodicDecoupling::new(&field, &group).expect("benzene decouples"),
        field,
        template: l.dynamics(),
        built: l.at.clone(),
        molecules: l.molecules.clone(),
        masses: l.masses.clone(),
        solvent: l.octanols,
        melt,
    }
}

/// The production's octanol leg: [`octanol::LIQUID`]'s 192 sites, `r_c` 9 Å, compensated, PME at
/// δ = 10⁻⁵ with the neighbour list ([`production`]), 4c-1's lattice seed.
pub fn octanol_leg() -> System {
    let (b, o) = (benzene(), octanol());
    let (l, field, _) = octanol_box(&b, &o, octanol::LIQUID, 9.0 * ANGSTROM, 1e-5, 0x0C71, true);
    let field = production(field, &l);
    octanol_system(&l, field, &OCTANOL_MELT)
}

/// [`System`] for benzene in water: `s`, whose field already carries benzene's model.
pub fn water_system(b: &Benzene, s: &Solvated, melt: &'static [MeltStage]) -> System {
    let mut molecules: Vec<Range<usize>> = Vec::with_capacity(1 + s.waters);
    molecules.push(0..b.at.len());
    for w in 0..s.waters {
        let o = b.at.len() + 3 * w;
        molecules.push(o..o + 3);
    }
    System {
        name: "water",
        decoupling: s.decoupling(),
        field: s.field.clone(),
        template: s.dynamics(b),
        built: s.at.clone(),
        molecules,
        masses: s.masses.clone(),
        solvent: s.waters,
        melt,
    }
}

/// `s` with benzene's model swapped in.
pub fn with_gaff_benzene(mut s: Solvated, b: &Benzene) -> Solvated {
    s.field = gaff::benzene_in_water(s.field.clone(), &b.component);
    s
}

/// The production's water leg: W3's box ([`hydration::w3_solvated`]) with benzene's model.
pub fn water_leg() -> System {
    let b = benzene();
    let s = with_gaff_benzene(hydration::w3_solvated(&b), &b);
    water_system(&b, &s, &WATER_MELT)
}

/// The pressure at `at`, pascals, by the molecular virial: 4c-1's ([`octanol::pressure`]) on any
/// set of whole molecules.
pub fn pressure(sys: &System, at: &[[f64; 3]]) -> f64 {
    let f = &sys.field;
    let ev = f.evaluate(at);
    let mut w = ev.virial[0][0] + ev.virial[1][1] + ev.virial[2][2];
    for r in &sys.molecules {
        let m: f64 = sys.masses[r.clone()].iter().sum();
        let com = [0, 1, 2].map(|k| r.clone().map(|i| sys.masses[i] * at[i][k]).sum::<f64>() / m);
        for i in r.clone() {
            w -= dot(sub(at[i], com), ev.forces[i]);
        }
    }
    let v = f.cell().volume();
    let n = sys.molecules.len() as f64;
    (n * BOLTZMANN.to_si() * KELVIN + w / 3.0) / v - f.dispersion_correction() / v
        + f.pressure_correction()
}

// ---------------------------------------------------------------------------------------------
// The melt

/// One stage of a melt: its state, length, step, bath and temperature.
#[derive(Clone, Copy, Debug)]
pub struct MeltStage {
    pub name: &'static str,
    pub lambda: Lambda,
    pub picoseconds: f64,
    pub dt: f64,
    pub friction: f64,
    pub temperature: f64,
}

pub const fn stage(
    name: &'static str,
    lambda: Lambda,
    picoseconds: f64,
    dt: f64,
    friction: f64,
    temperature: f64,
) -> MeltStage {
    MeltStage {
        name,
        lambda,
        picoseconds,
        dt,
        friction,
        temperature,
    }
}

pub const GROW_HALF: Lambda = Lambda {
    restraint: 0.0,
    electrostatics: 0.0,
    van_der_waals: 0.5,
};
pub const GROW_WHOLE: Lambda = Lambda {
    restraint: 0.0,
    electrostatics: 0.0,
    van_der_waals: 1.0,
};

/// **The octanol leg's melt**: benzene grown in as W3 grows it — 0.5 ps at 0.5 fs in a 50 ps⁻¹ bath
/// at each of `λ_v` = 0.5 and 1 with no charges, then coupled — then 4c-1's 10 ps at 500 K to lose
/// the lattice's order, and 60 ps at 298.15 K to settle, every 0.1 ps's `U` and `P` logged in 2.5 ps
/// blocks. 4c-1's GAFF liquid had `U/N`'s halves within 0.009 kcal/mol after 24 ps and τ 0.24 ps;
/// with benzene in and octanol eight times as viscous as water, sixty.
pub const OCTANOL_MELT: [MeltStage; 5] = [
    stage("grow λ_v 0.5", GROW_HALF, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("grow λ_v 1", GROW_WHOLE, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("coupled", Lambda::COUPLED, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("hot", Lambda::COUPLED, 10.0, 2.0 * FS, 5e12, 500.0),
    stage("settle", Lambda::COUPLED, 60.0, 2.0 * FS, 5e12, KELVIN),
];

/// **The water leg's melt**: W3's ([`hydration::w3_melt`]) stage for stage — the three growing
/// stages, then 20 ps at 2 fs in a 5 ps⁻¹ bath — with this file's seeds and its log.
pub const WATER_MELT: [MeltStage; 4] = [
    stage("grow λ_v 0.5", GROW_HALF, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("grow λ_v 1", GROW_WHOLE, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("coupled", Lambda::COUPLED, 0.5, 0.5 * FS, 50e12, KELVIN),
    stage("settle", Lambda::COUPLED, 20.0, 2.0 * FS, 5e12, KELVIN),
];

/// The melts' seeds.
pub const MELT_SEED: u64 = 0x4C2_E01;

/// A float as text that parses back to the same bits.
pub fn exact(x: f64) -> String {
    format!("{x:e}")
}

/// FNV-1a, 64-bit, over each float's bits, little-endian: W4's.
pub fn digest(xs: impl IntoIterator<Item = f64>) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for x in xs {
        for b in x.to_bits().to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

pub fn positions_digest(at: &[[f64; 3]]) -> u64 {
    digest(at.iter().flatten().copied())
}

/// **The model's digest**: every atom's `[x, D]` and charge, the mixing, the 1-4 scales, α, the
/// cutoff, `k_c` and the mesh's order and grid. What a window's header and the melt's say the
/// model was.
pub fn model_digest(f: &PeriodicForceField) -> u64 {
    let p = f.ewald().parameters();
    let mesh: Vec<f64> = f.ewald().mesh().map_or(Vec::new(), |m| {
        let q = m.parameters();
        let mut v = vec![q.order as f64];
        v.extend(q.grid.map(|g| g as f64));
        v
    });
    digest(
        f.van_der_waals()
            .iter()
            .flatten()
            .copied()
            .chain(f.charges().iter().copied())
            .chain([f64::from(u8::from(f.mixing() == Mixing::LorentzBerthelot))])
            .chain(f.one_four())
            .chain([p.alpha, p.cutoff, p.k_cutoff])
            .chain(mesh),
    )
}

/// The melt's header: every stage literally, the seed, the atom count, and digests of the box as
/// built and of the model.
pub fn melt_header(sys: &System, seed: u64) -> String {
    let st: Vec<String> = sys
        .melt
        .iter()
        .map(|s| {
            format!(
                "{}/{}/{}/{}/{}/{}/{}",
                exact(s.lambda.electrostatics),
                exact(s.lambda.van_der_waals),
                exact(s.picoseconds),
                exact(s.dt),
                exact(s.friction),
                exact(s.temperature),
                s.name.replace(' ', "_")
            )
        })
        .collect();
    format!(
        "# melted {} stages {} seed {seed:#x} atoms {} built {:016x} model {:016x}",
        sys.name,
        st.join(","),
        sys.built.len(),
        positions_digest(&sys.built),
        model_digest(&sys.field)
    )
}

/// Writes `header`, the positions [`exact`]ly and a trailer `# end N digest` to `path`'s
/// `.partial`, renamed once whole: W4's.
pub fn write_positions(path: &std::path::Path, header: &str, at: &[[f64; 3]]) {
    assert!(!header.contains('\n'), "a header is one line");
    let mut text = format!("{header}\n");
    for p in at {
        text.push_str(&format!(
            "{} {} {}\n",
            exact(p[0]),
            exact(p[1]),
            exact(p[2])
        ));
    }
    text.push_str(&format!(
        "# end {} {:016x}\n",
        at.len(),
        positions_digest(at)
    ));
    let partial = path.with_extension("partial");
    std::fs::write(&partial, text).expect("positions");
    std::fs::rename(&partial, path).expect("rename");
}

/// A positions file's header and positions, or why it is not whole: W4's.
pub fn parse_positions(text: &str) -> Result<(String, Vec<[f64; 3]>), String> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| "an empty file".to_string())?
        .to_string();
    let mut at = Vec::new();
    for line in lines {
        if let Some(end) = line.strip_prefix("# end ") {
            let mut words = end.split_whitespace();
            let n: Option<usize> = words.next().and_then(|x| x.parse().ok());
            let d = words.next().and_then(|x| u64::from_str_radix(x, 16).ok());
            return match (n, d) {
                (Some(n), _) if n != at.len() => {
                    Err(format!("{} positions where the trailer says {n}", at.len()))
                }
                (Some(_), Some(d)) if d == positions_digest(&at) => Ok((header, at)),
                (Some(_), Some(_)) => Err("the positions' digest is not the trailer's".into()),
                _ => Err(format!("a trailer cut short: {line:?}")),
            };
        }
        let v: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse).collect();
        match v {
            Ok(v) if v.len() == 3 => at.push([v[0], v[1], v[2]]),
            _ => return Err(format!("not a position: {line:?}")),
        }
    }
    Err("no trailer: the file is not whole".into())
}

/// Where the measurements write: `PANTOMETRY_OCTANOL_DIR`, or `pantometry-octanol` in the system's
/// temporary directory. Never the repository.
pub fn results_dir() -> std::path::PathBuf {
    let dir = std::env::var_os("PANTOMETRY_OCTANOL_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("pantometry-octanol"));
    std::fs::create_dir_all(&dir).expect("the results directory");
    dir
}

/// A lock that makes a second run in the same directory refuse: W4's. Removed when dropped; a run
/// that was killed leaves it, and says so.
pub struct Lock(std::path::PathBuf);

impl Lock {
    pub fn take(dir: &std::path::Path, name: &str) -> Lock {
        let path = dir.join(format!("{name}.lock"));
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap_or_else(|e| {
                panic!(
                    "{} exists ({e}): another run is writing here, or one was killed — remove it \
                     by hand once no run is",
                    path.display()
                )
            });
        Lock(path)
    }
}

impl Drop for Lock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Appends a line to `name.log` in `dir` and prints it.
pub fn log(dir: &std::path::Path, name: &str, line: &str) {
    use std::io::Write;
    println!("{line}");
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join(format!("{name}.log")))
        .expect("the log");
    writeln!(f, "{line}").expect("the log");
}

/// **`sys` melted from its positions as built**, read back from `start.txt` in `dir` when a run
/// wrote it under the same header ([`melt_header`]); a file melted differently, or not whole, is
/// refused by name. Otherwise run now under a lock `melt`: each stage from the last, the bath
/// thermalised once at the start, every 0.1 ps's `U` per solvent molecule, `P` and `T` logged to
/// `melt.log` in 2.5 ps blocks, and the end written atomically.
pub fn melted(sys: &System, dir: &std::path::Path, seed: u64) -> Vec<[f64; 3]> {
    let path = dir.join("start.txt");
    let want = melt_header(sys, seed);
    if let Ok(text) = std::fs::read_to_string(&path) {
        let (head, at) = parse_positions(&text).unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it to melt again",
                path.display()
            )
        });
        assert_eq!(
            head,
            want,
            "{} was melted differently: remove it, and every window run from it, to melt again",
            path.display()
        );
        assert_eq!(at.len(), sys.built.len());
        return at;
    }
    let _lock = Lock::take(dir, "melt");
    let t = std::time::Instant::now();
    let mut at = sys.built.clone();
    let mut md = sys.template.clone();
    let n = sys.solvent as f64;
    let mut rows: Vec<(String, [f64; 4])> = Vec::new();
    let mut time = 0.0;
    for (k, st) in sys.melt.iter().enumerate() {
        let s = window_seed(seed, k);
        md = md.with_bath(Bath::Langevin {
            temperature: st.temperature,
            friction: st.friction,
            seed: s,
        });
        if k == 0 {
            md = md.thermalised(&at, st.temperature, s);
        }
        md.forget_forces();
        let potential = AtLambda {
            hamiltonian: &sys.decoupling,
            lambda: st.lambda,
        };
        let steps = (st.picoseconds * 1e-12 / st.dt).round() as usize;
        let every = ((0.1e-12 / st.dt).round() as usize).max(1);
        let t0 = std::time::Instant::now();
        for i in 0..steps {
            md.step(&potential, &mut at, st.dt);
            if (i + 1) % every == 0 {
                time += every as f64 * st.dt;
                rows.push((
                    st.name.to_string(),
                    [
                        time * 1e12,
                        kcal(md.potential_energy().expect("stepped")) / n,
                        pressure(sys, &at) / 1e5,
                        md.half_step_temperature(),
                    ],
                ));
            }
        }
        let secs = t0.elapsed().as_secs_f64();
        log(
            dir,
            "melt",
            &format!(
                "{} stage {k} ({}): {} ps at {} fs, {} K, γ {} /ps, λ {:?}: {secs:.0} s, {:.1} ms a \
                 step",
                sys.name,
                st.name,
                st.picoseconds,
                st.dt / FS,
                st.temperature,
                st.friction / 1e12,
                st.lambda,
                1e3 * secs / steps as f64
            ),
        );
    }
    assert!(rows.iter().all(|r| r.1.iter().all(|x| x.is_finite())));
    log(
        dir,
        "melt",
        &format!(
            "{}: the melt in blocks of 2.5 ps\n| stage | t (ps) | U/N (kcal/mol) | P (bar) | T (K) |\n| --- | --- | --- | --- | --- |",
            sys.name
        ),
    );
    for block in rows.chunks(25) {
        let m = |c: usize| block.iter().map(|r| r.1[c]).sum::<f64>() / block.len() as f64;
        log(
            dir,
            "melt",
            &format!(
                "| {} | {:.1}–{:.1} | {:.4} | {:.0} | {:.1} |",
                block[block.len() - 1].0,
                block[0].1[0] - 0.1,
                block[block.len() - 1].1[0],
                m(1),
                m(2),
                m(3)
            ),
        );
    }
    // The last stage's halves, as a check on drift.
    let last = sys.melt.len() - 1;
    let tail: Vec<&[f64; 4]> = rows
        .iter()
        .filter(|r| r.0 == sys.melt[last].name)
        .map(|r| &r.1)
        .collect();
    if tail.len() >= 4 {
        let h = tail.len() / 2;
        let est =
            |c: usize, r: &[&[f64; 4]]| Estimate::of(&r.iter().map(|x| x[c]).collect::<Vec<_>>());
        for (c, what) in [(1, "U/N (kcal/mol)"), (2, "P (bar)")] {
            let (a, b, all) = (est(c, &tail[..h]), est(c, &tail[h..]), est(c, &tail));
            log(
                dir,
                "melt",
                &format!(
                    "{}: {what} over the last stage {:.4} ± {:.4} (τ {:.1} samples of 0.1 ps); \
                     halves {:.4} ± {:.4} and {:.4} ± {:.4}",
                    sys.name, all.mean, all.error, all.tau, a.mean, a.error, b.mean, b.error
                ),
            );
        }
    }
    write_positions(&path, &want, &at);
    log(
        dir,
        "melt",
        &format!("{} melted in {:.0} s", sys.name, t.elapsed().as_secs_f64()),
    );
    at
}

// ---------------------------------------------------------------------------------------------
// The windows

/// The candidate states, W3's, in path order: the charges off by eighths at full van der Waals,
/// then the van der Waals off by twentieths with no charge.
pub fn grid() -> Vec<Lambda> {
    let mut g = Vec::new();
    for e in [1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0] {
        g.push(Lambda::new(0.0, e, 1.0));
    }
    for k in 1..=20 {
        g.push(Lambda::new(0.0, 0.0, f64::from(20 - k) / 20.0));
    }
    g
}

/// W3's first schedule: the charges off in two intervals and the van der Waals in ten.
pub const SCHEDULE: [usize; 13] = [0, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];

/// An interval overlapping below this gets a window at the candidate between: W3's.
pub const OVERLAP_THRESHOLD: f64 = 0.1;

/// W3's two rounds of insertion.
pub const INSERTION_ROUNDS: usize = 2;

/// W3's protocol: 2 fs, a 1 ps⁻¹ bath at 298.15 K, 20 ps discarded, 2000 samples every 100 fs.
/// Each leg gives it its own seed ([`seeded`]).
pub const LONG: Protocol = Protocol {
    time_step: 2e-15,
    temperature: KELVIN,
    friction: 1e12,
    equilibration: 10_000,
    stride: 50,
    samples: 2000,
    seed: 0,
};

/// Variant 4's short protocol: 4 ps discarded, 320 samples every 50 fs, 16 ps.
pub const SHORT: Protocol = Protocol {
    equilibration: 2_000,
    stride: 25,
    samples: 320,
    ..LONG
};

/// The pilot that measures τ: 2 ps discarded, 300 samples every 100 fs, 30 ps.
pub const PILOT: Protocol = Protocol {
    equilibration: 1_000,
    stride: 50,
    samples: 300,
    ..LONG
};

/// `p` with leg `name`'s seed: the water's and the octanol's apart from each other and from W3's.
pub fn seeded(p: Protocol, name: &str) -> Protocol {
    Protocol {
        seed: if name == "water" {
            0x4C2_3A7E
        } else {
            0x4C2_0C71
        },
        ..p
    }
}

/// What one window recorded, read back from its file: W3's and W4's.
#[derive(Clone, Debug)]
pub struct Recorded {
    pub grid: usize,
    pub lambda: Lambda,
    pub gradients: Vec<[f64; 3]>,
    /// Per sample, the coupling energy at every candidate, joules.
    pub energies: Vec<Vec<f64>>,
    /// Per sample, what [`window_observed`] was given to observe of the positions — the decane
    /// leg's nematic order — and empty for a window run without one.
    pub observed: Vec<f64>,
    pub seconds: f64,
}

impl Recorded {
    pub fn delta(&self, to: usize, kt: f64) -> Vec<f64> {
        self.energies
            .iter()
            .map(|e| (e[to] - e[self.grid]) / kt)
            .collect()
    }

    pub fn slice(&self, range: Range<usize>) -> Recorded {
        Recorded {
            gradients: self.gradients[range.clone()].to_vec(),
            energies: self.energies[range.clone()].to_vec(),
            observed: if self.observed.is_empty() {
                Vec::new()
            } else {
                self.observed[range].to_vec()
            },
            ..self.clone()
        }
    }
}

pub fn interval_bar(a: &Recorded, b: &Recorded, kt: f64) -> pantometry_forcefield::Bar {
    let reverse: Vec<f64> = b.delta(a.grid, kt).iter().map(|x| -x).collect();
    pantometry_forcefield::free_energy::bar_correlated(&a.delta(b.grid, kt), &reverse)
}

pub fn chain(records: &[Recorded], kt: f64) -> pantometry_forcefield::FreeEnergy {
    let n = records.len();
    let windows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|j| {
            let prev = if j > 0 {
                records[j].delta(records[j - 1].grid, kt)
            } else {
                Vec::new()
            };
            let next = if j + 1 < n {
                records[j].delta(records[j + 1].grid, kt)
            } else {
                Vec::new()
            };
            (prev, next)
        })
        .collect();
    bennett_chain(&windows)
}

/// TI on the trapezoid, kcal/mol, with its error and each window's term.
pub fn trapezoid(records: &[Recorded]) -> (f64, f64, Vec<Estimate>) {
    let n = records.len();
    let terms: Vec<Estimate> = (0..n)
        .map(|k| {
            let lo = records[k.saturating_sub(1)].lambda.components();
            let hi = records[(k + 1).min(n - 1)].lambda.components();
            let w = [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]));
            let ys: Vec<f64> = records[k]
                .gradients
                .iter()
                .map(|g| (0..3).map(|c| g[c] * w[c]).sum())
                .collect();
            Estimate::of(&ys)
        })
        .collect();
    (
        kcal(terms.iter().map(|t| t.mean).sum()),
        kcal(terms.iter().map(|t| t.error * t.error).sum::<f64>().sqrt()),
        terms,
    )
}

/// One leg's inputs: its decoupling, the dynamics, the melted start, the model's digest, the
/// protocol, the candidates and the log's name.
pub struct Leg<'a> {
    pub name: &'static str,
    pub decoupling: &'a PeriodicDecoupling,
    pub template: &'a MolecularDynamics,
    pub start: Vec<[f64; 3]>,
    pub model: u64,
    pub protocol: Protocol,
    pub grid: Vec<Lambda>,
}

impl<'a> Leg<'a> {
    pub fn new(
        name: &'static str,
        decoupling: &'a PeriodicDecoupling,
        template: &'a MolecularDynamics,
        start: Vec<[f64; 3]>,
        protocol: Protocol,
        grid: Vec<Lambda>,
    ) -> Leg<'a> {
        Leg {
            name,
            model: model_digest(decoupling.field()),
            decoupling,
            template,
            start,
            protocol,
            grid,
        }
    }
}

/// Window `g`'s header: the leg, the candidate and its λ, the protocol, the atom count, the number
/// of candidates, and digests of the start and of the model.
pub fn header(leg: &Leg, g: usize) -> String {
    let p = leg.protocol;
    let l = leg.grid[g];
    format!(
        "# {} candidate {g} lambda {} {} {} protocol {} {} {} {} {} {} {:#x} atoms {} states {} \
         start {:016x} model {:016x}",
        leg.name,
        exact(l.restraint),
        exact(l.electrostatics),
        exact(l.van_der_waals),
        exact(p.time_step),
        exact(p.temperature),
        exact(p.friction),
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        leg.start.len(),
        leg.grid.len(),
        positions_digest(&leg.start),
        leg.model
    )
}

pub fn parse(text: &str) -> Recorded {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().expect("a header").split_whitespace().collect();
    let field = |name: &str| head.iter().position(|&w| w == name).expect(name) + 1;
    let grid_index: usize = head[field("candidate")].parse().unwrap();
    let l = field("lambda");
    let lambda = Lambda::new(
        head[l].parse().unwrap(),
        head[l + 1].parse().unwrap(),
        head[l + 2].parse().unwrap(),
    );
    let states: usize = head[field("states")].parse().unwrap();
    let (mut gradients, mut energies, mut seconds) = (Vec::new(), Vec::new(), f64::NAN);
    let mut observed = Vec::new();
    for line in lines {
        if let Some(rest) = line.strip_prefix("# done seconds ") {
            seconds = rest.trim().parse().unwrap();
            continue;
        }
        if let Some(rest) = line.strip_prefix("# observed ") {
            let v: Vec<&str> = rest.split_whitespace().collect();
            assert_eq!(v.len(), 2, "an observation line");
            observed.push(v[1].parse().unwrap());
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(v.len(), 3 + states, "a sample line");
        gradients.push([v[0], v[1], v[2]]);
        energies.push(v[3..].to_vec());
    }
    assert!(seconds.is_finite(), "a complete window says so");
    Recorded {
        grid: grid_index,
        lambda,
        gradients,
        energies,
        observed,
        seconds,
    }
}

/// One window as `Windows` would run it at schedule position `index`: W3's and W4's.
#[allow(clippy::too_many_arguments)]
pub fn sample_window(
    h: &PeriodicDecoupling,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    p: Protocol,
    index: usize,
    lambda: Lambda,
    states: &[Lambda],
    mut on: impl FnMut(u64, &[Coupling], &[[f64; 3]]),
) -> Vec<[f64; 3]> {
    let seed = window_seed(p.seed, index);
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: p.temperature,
            friction: p.friction,
            seed,
        })
        .thermalised(start, p.temperature, seed);
    let potential = AtLambda {
        hamiltonian: h,
        lambda,
    };
    let mut at = start.to_vec();
    md.prepare(&potential, &at);
    let mut lambdas = vec![lambda];
    lambdas.extend_from_slice(states);
    let mut taken = 0;
    // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
    // builds on 1.78.
    #[allow(clippy::manual_is_multiple_of)]
    while taken < p.samples {
        md.step(&potential, &mut at, p.time_step);
        let s = md.steps();
        if s > p.equilibration && (s - p.equilibration) % p.stride == 0 {
            on(s, &h.couplings(&at, &lambdas), &at);
            taken += 1;
        }
    }
    at
}

/// Window `g`'s record: from `window_NN.txt` in `dir` when a run completed it under the same
/// header, refused by name when it was run differently; or run now, written sample by sample to
/// `window_NN.partial`, begun afresh whatever a killed run left there, and renamed once whole. W4's.
pub fn window(leg: &Leg, dir: &std::path::Path, g: usize) -> Recorded {
    window_observed(leg, dir, g, None)
}

/// What [`window_observed`] observes of each sample's positions.
pub type Observe<'a> = &'a dyn Fn(&[[f64; 3]]) -> f64;

/// [`window`], and with `observe` each sample's positions observed too: a line `# observed STEP
/// VALUE` after the sample's, read back into [`Recorded::observed`]. A file read back must then
/// hold one observation a sample, or it is refused by name. Without, the file is [`window`]'s to
/// the byte.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
pub fn window_observed(
    leg: &Leg,
    dir: &std::path::Path,
    g: usize,
    observe: Option<Observe>,
) -> Recorded {
    use std::io::Write;
    let path = dir.join(format!("window_{g:02}.txt"));
    let p = leg.protocol;
    let want = header(leg, g);
    if let Ok(text) = std::fs::read_to_string(&path) {
        assert_eq!(
            text.lines().next().unwrap_or(""),
            want,
            "{} was run differently",
            path.display()
        );
        let r = parse(&text);
        assert_eq!(r.energies.len(), p.samples);
        if observe.is_some() {
            assert_eq!(
                r.observed.len(),
                p.samples,
                "{} was run without its observation",
                path.display()
            );
        }
        return r;
    }
    let t = std::time::Instant::now();
    let partial = path.with_extension("partial");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
    writeln!(file, "{want}").unwrap();
    let total = p.steps_per_window();
    log(
        dir,
        leg.name,
        &format!(
            "{} window at candidate {g} λ = {:?} begins: {total} steps",
            leg.name, leg.grid[g]
        ),
    );
    let mut taken = 0usize;
    sample_window(
        leg.decoupling,
        &leg.start,
        leg.template,
        p,
        g,
        leg.grid[g],
        &leg.grid,
        |s, c, at| {
            let mut line = format!(
                "{s} {} {} {}",
                exact(c[0].gradient[0]),
                exact(c[0].gradient[1]),
                exact(c[0].gradient[2])
            );
            for x in &c[1..] {
                line.push(' ');
                line.push_str(&exact(x.energy));
            }
            writeln!(file, "{line}").unwrap();
            if let Some(f) = observe {
                writeln!(file, "# observed {s} {}", exact(f(at))).unwrap();
            }
            taken += 1;
            if taken % 100 == 0 {
                file.flush().unwrap();
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    leg.name,
                    &format!(
                        "{} candidate {g}: {taken} of {} samples, step {s} of {total}, {secs:.0} \
                         s, {:.1} ms a step",
                        leg.name,
                        p.samples,
                        secs * 1e3 / s as f64
                    ),
                );
            }
        },
    );
    let seconds = t.elapsed().as_secs_f64();
    writeln!(file, "# done seconds {}", exact(seconds)).unwrap();
    file.flush().unwrap();
    drop(file);
    std::fs::rename(&partial, &path).expect("rename");
    log(
        dir,
        leg.name,
        &format!(
            "{} window at candidate {g} done in {seconds:.0} s",
            leg.name
        ),
    );
    parse(&std::fs::read_to_string(&path).unwrap())
}

/// What a panic said.
#[cfg(not(target_family = "wasm"))]
pub fn refusal(f: impl FnOnce()) -> String {
    let e = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).expect_err("refused");
    e.downcast_ref::<String>()
        .cloned()
        .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_default()
}

/// Every gradient and energy a record holds, as bits.
#[cfg(not(target_family = "wasm"))]
pub fn record_bits(r: &Recorded) -> Vec<u64> {
    let mut v: Vec<u64> = r.gradients.iter().flatten().map(|x| x.to_bits()).collect();
    for e in &r.energies {
        v.extend(e.iter().map(|x| x.to_bits()));
    }
    v
}

#[cfg(not(target_family = "wasm"))]
pub fn is_empty(dir: &std::path::Path) -> bool {
    std::fs::read_dir(dir)
        .expect("a directory")
        .next()
        .is_none()
}

/// One leg's result, kcal/mol: `ΔG_decouple` by BAR along the chain and its error, TI and its
/// error, the halves, and the lowest overlap.
pub struct LegResult {
    pub bar: (f64, f64),
    pub ti: (f64, f64),
    pub halves: [f64; 2],
    pub overlap: f64,
}

/// The chain's BAR over each window's first half of samples and over its second, kcal/mol.
pub fn halves_of(records: &[Recorded], kt: f64) -> [f64; 2] {
    let half = records[0].energies.len() / 2;
    let first: Vec<Recorded> = records.iter().map(|r| r.slice(0..half)).collect();
    let second: Vec<Recorded> = records
        .iter()
        .map(|r| r.slice(half..r.energies.len()))
        .collect();
    [
        kcal(kt * chain(&first, kt).value),
        kcal(kt * chain(&second, kt).value),
    ]
}

/// The leg's table — each window's λ, TI term, τ and seconds, each interval's BAR and overlap —
/// and its totals, logged.
pub fn report(dir: &std::path::Path, name: &str, records: &[Recorded], kt: f64) -> LegResult {
    let to_kcal = |x: f64| kcal(kt * x);
    let (ti, ti_error, terms) = trapezoid(records);
    log(
        dir,
        name,
        "| window | candidate | λ_e | λ_v | TI term (kcal/mol) | τ (samples) | wall (s) | → next: BAR (kcal/mol) | overlap |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    );
    let mut overlap = f64::INFINITY;
    for k in 0..records.len() {
        let r = &records[k];
        let next = if k + 1 < records.len() {
            let b = interval_bar(r, &records[k + 1], kt);
            overlap = overlap.min(b.overlap);
            format!(
                "{:+.3} ± {:.3} | {:.3}",
                to_kcal(b.delta),
                to_kcal(b.error()),
                b.overlap
            )
        } else {
            "| ".to_string()
        };
        log(
            dir,
            name,
            &format!(
                "| {k} | {} | {} | {} | {:+.3} ± {:.3} | {:.1} | {:.0} | {next} |",
                r.grid,
                r.lambda.electrostatics,
                r.lambda.van_der_waals,
                kcal(terms[k].mean),
                kcal(terms[k].error),
                terms[k].tau,
                r.seconds
            ),
        );
    }
    let total = chain(records, kt);
    let halves = halves_of(records, kt);
    let out = LegResult {
        bar: (to_kcal(total.value), to_kcal(total.error)),
        ti: (ti, ti_error),
        halves,
        overlap,
    };
    log(
        dir,
        name,
        &format!(
            "{name}: ΔG_decouple BAR {:+.3} ± {:.3}, TI {ti:+.3} ± {ti_error:.3} kcal/mol; halves \
             {:+.3}, {:+.3}; lowest overlap {overlap:.3}; ΔG_solv = {:+.3}",
            out.bar.0, out.bar.1, halves[0], halves[1], -out.bar.0
        ),
    );
    out
}

/// The candidates whose windows a finished run left whole in `dir`, ascending: its schedule,
/// whatever the insertion added.
pub fn finished_schedule(dir: &std::path::Path) -> Vec<usize> {
    (0..grid().len())
        .filter(|g| dir.join(format!("window_{g:02}.txt")).is_file())
        .collect()
}

/// Window `g` as a finished run left it in `dir`, header-checked by [`window`]; **refused by name
/// when it is not there**, so that reading a finished run never runs one of its windows.
pub fn read_back(leg: &Leg, dir: &std::path::Path, g: usize) -> Recorded {
    let path = dir.join(format!("window_{g:02}.txt"));
    assert!(
        path.is_file(),
        "{} is not there: a finished run is read, not run",
        path.display()
    );
    window(leg, dir, g)
}

/// `sys`'s `start.txt` in `leg_dir`, header-checked by [`melted`]; **refused by name when it is
/// not there**, so that no melt is run in its place.
pub fn read_back_start(sys: &System, leg_dir: &std::path::Path) -> Vec<[f64; 3]> {
    let path = leg_dir.join("start.txt");
    assert!(
        path.is_file(),
        "{} is not there: the start is read, not melted",
        path.display()
    );
    melted(sys, leg_dir, MELT_SEED)
}

/// A segment's BAR along its chain, kcal/mol, and its error.
pub fn segment(records: &[Recorded], kt: f64) -> (f64, f64) {
    let c = chain(records, kt);
    (kcal(kt * c.value), kcal(kt * c.error))
}
