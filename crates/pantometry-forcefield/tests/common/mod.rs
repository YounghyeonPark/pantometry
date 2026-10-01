//! Helpers shared by the minimisation tests: building a component, geometry, a Hessian's smallest
//! eigenvalue, and a best-fit RMSD. Each is plain arithmetic that the tests can check by hand.

#![allow(dead_code)]
// The Jacobi rotations and the Hessian are written as the textbook's index loops, which read
// against the formulas; iterator chains over two rows at once would not.
#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::{Component, ForceField};

pub const ANGSTROM: f64 = 1e-10;
pub const DEG: f64 = std::f64::consts::PI / 180.0;

/// A one-component CCD file: atoms as `(name, element, [x, y, z] in Å)`, bonds as
/// `(first, second, order, aromatic)`. Every atom's aromatic flag is the OR of its bonds'.
pub fn entry(atoms: &[(&str, &str, [f64; 3])], bonds: &[(&str, &str, &str, bool)]) -> Component {
    let yn = |b: bool| if b { "Y" } else { "N" };
    let mut s = String::from(
        "data_TST\n_chem_comp.id TST\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (name, element, [x, y, z]) in atoms {
        let aromatic = bonds
            .iter()
            .any(|(a, b, _, ar)| *ar && (a == name || b == name));
        s += &format!("TST {name} {element} 0 {} {x} {y} {z}\n", yn(aromatic));
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (a, b, order, aromatic) in bonds {
        s += &format!("{a} {b} {order} {}\n", yn(*aromatic));
    }
    Component::from_ccd(&s).expect("the test entry parses")
}

pub fn positions(c: &Component) -> Vec<[f64; 3]> {
    c.atoms().iter().map(|a| a.at).collect()
}

/// The index of the atom called `name`.
pub fn index(c: &Component, name: &str) -> usize {
    c.atoms()
        .iter()
        .position(|a| a.name == name)
        .unwrap_or_else(|| panic!("no atom {name}"))
}

pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Distance between atoms `i` and `j`, Å.
pub fn distance(at: &[[f64; 3]], i: usize, j: usize) -> f64 {
    len(sub(at[i], at[j])) / ANGSTROM
}

/// The angle at `j`, degrees.
pub fn angle(at: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let (u, v) = (sub(at[i], at[j]), sub(at[k], at[j]));
    (dot(u, v) / (len(u) * len(v))).clamp(-1.0, 1.0).acos() / DEG
}

/// The largest per-atom force, kcal mol⁻¹ Å⁻¹.
pub fn max_force(ff: &ForceField, at: &[[f64; 3]]) -> f64 {
    let to = ANGSTROM / pantometry_forcefield::uff::KCAL_PER_MOL;
    ff.evaluate(at)
        .forces
        .iter()
        .map(|f| len(*f) * to)
        .fold(0.0, f64::max)
}

/// The eigenvalues of a symmetric matrix by cyclic Jacobi rotations, ascending, with the
/// eigenvectors as the columns of the second result. Converges to the rounding of the matrix.
pub fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        let scale: f64 = (0..n).map(|i| a[i][i] * a[i][i]).sum::<f64>() + off;
        if off <= 1e-30 * scale {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.iter().map(|&i| a[i][i]).collect();
    let vectors = (0..n)
        .map(|r| order.iter().map(|&i| v[r][i]).collect())
        .collect();
    (values, vectors)
}

/// The smallest Hessian eigenvalue of `ff` at `at` once the rigid motions are projected out,
/// kcal mol⁻¹ Å⁻², and the largest. The Hessian is the central difference of the analytic
/// forces at h = 1e-5 Å, symmetrised; translations and rotations about the centroid are
/// projected out (P H P with P the projector onto their orthogonal complement), which leaves them
/// as exact zeros, and the smallest eigenvalue above 1e-6 of the largest is returned.
pub fn hessian_extremes(ff: &ForceField, at: &[[f64; 3]]) -> (f64, f64) {
    let n = at.len();
    let to = ANGSTROM / pantometry_forcefield::uff::KCAL_PER_MOL;
    let h = 1e-5;
    let mut hess = vec![vec![0.0; 3 * n]; 3 * n];
    for a in 0..n {
        for c in 0..3 {
            let mut p = at.to_vec();
            let mut m = at.to_vec();
            p[a][c] += h * ANGSTROM;
            m[a][c] -= h * ANGSTROM;
            let (fp, fm) = (ff.evaluate(&p).forces, ff.evaluate(&m).forces);
            for b in 0..n {
                for d in 0..3 {
                    hess[3 * a + c][3 * b + d] = -(fp[b][d] - fm[b][d]) * to / (2.0 * h);
                }
            }
        }
    }
    for i in 0..3 * n {
        for j in i + 1..3 * n {
            let s = 0.5 * (hess[i][j] + hess[j][i]);
            hess[i][j] = s;
            hess[j][i] = s;
        }
    }
    // Rigid modes, orthonormalised.
    let centre: [f64; 3] = {
        let mut c = [0.0; 3];
        for p in at {
            for k in 0..3 {
                c[k] += p[k] / ANGSTROM / n as f64;
            }
        }
        c
    };
    let mut modes: Vec<Vec<f64>> = Vec::new();
    for k in 0..3 {
        let mut v = vec![0.0; 3 * n];
        for a in 0..n {
            v[3 * a + k] = 1.0;
        }
        modes.push(v);
    }
    for k in 0..3 {
        let axis = [
            if k == 0 { 1.0 } else { 0.0 },
            if k == 1 { 1.0 } else { 0.0 },
            if k == 2 { 1.0 } else { 0.0 },
        ];
        let mut v = vec![0.0; 3 * n];
        for a in 0..n {
            let r = [
                at[a][0] / ANGSTROM - centre[0],
                at[a][1] / ANGSTROM - centre[1],
                at[a][2] / ANGSTROM - centre[2],
            ];
            let w = cross(axis, r);
            v[3 * a..3 * a + 3].copy_from_slice(&w);
        }
        modes.push(v);
    }
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for mut v in modes {
        for b in &basis {
            let d: f64 = v.iter().zip(b).map(|(x, y)| x * y).sum();
            for (x, y) in v.iter_mut().zip(b) {
                *x -= d * y;
            }
        }
        let l = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if l > 1e-8 {
            basis.push(v.iter().map(|x| x / l).collect());
        }
    }
    let m = 3 * n;
    let mut proj = vec![vec![0.0; m]; m];
    for (i, row) in proj.iter_mut().enumerate() {
        row[i] = 1.0;
        for b in &basis {
            for (j, x) in row.iter_mut().enumerate() {
                *x -= b[i] * b[j];
            }
        }
    }
    let mul = |a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>| {
        let mut o = vec![vec![0.0; m]; m];
        for i in 0..m {
            for k in 0..m {
                if a[i][k] == 0.0 {
                    continue;
                }
                for j in 0..m {
                    o[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        o
    };
    let projected = mul(&proj, &mul(&hess, &proj));
    let (values, _) = jacobi(projected);
    let largest = values.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let smallest = values
        .iter()
        .copied()
        .filter(|v| v.abs() > 1e-6 * largest)
        .fold(f64::INFINITY, f64::min);
    (smallest, largest)
}

/// The RMSD between two sets of points (Å) after the best rigid superposition, by Horn's
/// quaternion method: the largest eigenvalue of a 4×4 symmetric matrix built from the
/// cross-covariance gives the minimum of Σ|R a + t − b|² directly.
pub fn rmsd(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let n = a.len() as f64;
    let centroid = |p: &[[f64; 3]]| {
        let mut c = [0.0; 3];
        for q in p {
            for k in 0..3 {
                c[k] += q[k] / n;
            }
        }
        c
    };
    let (ca, cb) = (centroid(a), centroid(b));
    let a: Vec<[f64; 3]> = a.iter().map(|p| sub(*p, ca)).collect();
    let b: Vec<[f64; 3]> = b.iter().map(|p| sub(*p, cb)).collect();
    let mut s = [[0.0; 3]; 3];
    let mut ga = 0.0;
    let mut gb = 0.0;
    for (p, q) in a.iter().zip(&b) {
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += p[i] * q[j];
            }
        }
        ga += dot(*p, *p);
        gb += dot(*q, *q);
    }
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    let k = vec![
        vec![sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        vec![syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        vec![szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        vec![sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let (values, _) = jacobi(k);
    let top = values[3];
    ((ga + gb - 2.0 * top).max(0.0) / n).sqrt()
}

/// Rotates everything on `k`'s side of the bond `j`–`k` (found by walking the bonds without
/// crossing it) about the `j`→`k` axis so that the dihedral `[i, j, k, l]` becomes `target`.
pub fn set_dihedral(c: &Component, at: &mut [[f64; 3]], atoms: [usize; 4], target: f64) {
    let [_, j, k, _] = atoms;
    let n = at.len();
    let mut side = vec![false; n];
    let mut stack = vec![k];
    side[k] = true;
    while let Some(a) = stack.pop() {
        for (b, _) in c.neighbours(a) {
            if (a == k && b == j) || side[b] {
                continue;
            }
            assert_ne!(b, j, "the bond {j}-{k} is in a ring");
            side[b] = true;
            stack.push(b);
        }
    }
    let now = pantometry_forcefield::minimise::dihedral(at, atoms);
    rotate(at, &side, j, k, target - now);
    let got = pantometry_forcefield::minimise::dihedral(at, atoms);
    let miss = (got - target + 3.0 * std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::PI;
    if miss.abs() > 1e-9 {
        // The rotation's sense is the other one in this convention: undo and turn the other way.
        rotate(at, &side, j, k, -2.0 * (target - now));
    }
    let got = pantometry_forcefield::minimise::dihedral(at, atoms);
    let miss = (got - target + 3.0 * std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::PI;
    assert!(miss.abs() < 1e-9, "set_dihedral missed by {miss}");
}

/// Rotates the atoms marked in `side` (except `k`) by `turn` about the line from `j` to `k`
/// (Rodrigues' formula).
fn rotate(at: &mut [[f64; 3]], side: &[bool], j: usize, k: usize, turn: f64) {
    let axis = sub(at[k], at[j]);
    let l = len(axis);
    let u = [axis[0] / l, axis[1] / l, axis[2] / l];
    let (s, co) = turn.sin_cos();
    let pivot = at[k];
    for (a, p) in at.iter_mut().enumerate() {
        if !side[a] || a == k {
            continue;
        }
        let v = sub(*p, pivot);
        let uv = cross(u, v);
        let d = dot(u, v);
        for m in 0..3 {
            p[m] = pivot[m] + v[m] * co + uv[m] * s + u[m] * d * (1.0 - co);
        }
    }
}

/// The restraint stiffness the scans hold a dihedral with: 1e4 kcal mol⁻¹ rad⁻². At a point where
/// the molecule's own torque on the dihedral is τ, the dihedral sits τ/k off its target and the
/// restraint holds τ²/(2k); a scan reports that energy beside each point, and at a barrier's top
/// or bottom τ is zero by definition.
pub const RESTRAINT: f64 = 1e4;

/// The force the scans minimise to: 1e-5 kcal mol⁻¹ Å⁻¹. Its energy error is `|F|²/(2 λ_min)`,
/// below 1e-9 kcal/mol for any λ_min above 0.05 kcal mol⁻¹ Å⁻².
pub const SCAN_TOLERANCE: f64 = 1e-5;

/// One restrained point: the dihedral `atoms` set to `phi` degrees and held there while everything
/// else relaxes from `start`. Returns the force field's energy and the restraint's (kcal/mol) and
/// the relaxed positions.
pub fn held(
    ff: &ForceField,
    c: &Component,
    start: &[[f64; 3]],
    atoms: [usize; 4],
    phi: f64,
) -> (f64, f64, Vec<[f64; 3]>) {
    use pantometry_forcefield::minimise::{DihedralRestraint, KCAL_PER_MOL_ANGSTROM};
    use pantometry_forcefield::uff::KCAL_PER_MOL;
    let mut at = start.to_vec();
    set_dihedral(c, &mut at, atoms, phi * DEG);
    let r = DihedralRestraint {
        atoms,
        target: phi * DEG,
        stiffness: RESTRAINT * KCAL_PER_MOL,
    };
    let p = ff.minimise_restrained(
        &mut at,
        &[r],
        50_000,
        SCAN_TOLERANCE * KCAL_PER_MOL_ANGSTROM,
    );
    finished(&p, phi);
    (p.energy / KCAL_PER_MOL, p.restraint / KCAL_PER_MOL, at)
}

/// A scan point is finished when it converged to [`SCAN_TOLERANCE`], or stalled — no step lowers
/// the energy in floating point — with the largest per-atom force at most [`STALL_LIMIT`]. Either
/// way its energy is within `|F|² / (2 λ_min)` of the minimum, where `|F|` is the **whole** force
/// vector, at most `√N` times the per-atom maximum: so within `N · STALL_LIMIT² / (2 λ_min)`, which
/// for the 18 atoms of the largest molecule scanned and the smallest λ_min measured at any of
/// their minima (see `the_papers_barriers`) is below 1e-5 kcal/mol.
pub fn finished(p: &pantometry_forcefield::Progress, phi: f64) {
    use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
    use pantometry_forcefield::Status;
    let f = p.max_force / KCAL_PER_MOL_ANGSTROM;
    assert!(
        p.status == Status::Converged || (p.status == Status::Stalled && f <= STALL_LIMIT),
        "held at {phi}°: {p:?}"
    );
}

/// The largest force, kcal mol⁻¹ Å⁻¹, a stalled scan point may stop at.
pub const STALL_LIMIT: f64 = 1e-4;

/// The sum of the magnitudes of every term of `ff` at `at`, kcal/mol, and how many terms there
/// are — what the rounding of one evaluation of the energy scales with.
pub fn term_magnitude(ff: &ForceField, at: &[[f64; 3]]) -> (f64, usize) {
    use pantometry_forcefield::uff::KCAL_PER_MOL;
    let mut sum = 0.0;
    let mut n = 0;
    for s in ff.stretches() {
        sum += s.energy(len(sub(at[s.atoms[0]], at[s.atoms[1]]))).abs();
        n += 1;
    }
    for b in ff.bends() {
        sum += b.energy_at(at).abs();
        n += 1;
    }
    for t in ff.torsions() {
        sum += t.energy_at(at).abs();
        n += 1;
    }
    for v in ff.inversions() {
        sum += v.energy_at(at).abs();
        n += 1;
    }
    for p in ff.pairs() {
        sum += p.energy(len(sub(at[p.atoms[0]], at[p.atoms[1]]))).abs();
        n += 1;
    }
    (sum / KCAL_PER_MOL, n)
}

/// **The force below which a correct minimiser may stall**, kcal mol⁻¹ Å⁻¹. A step along a good
/// direction from a point with force `F` lowers the energy by about `|F|² / (2 λ)` for the
/// curvature λ along it, at least `|F|² / (2 λ_max)`. The Armijo comparison sees that decrease
/// only if it exceeds the rounding of the two energies compared, `2 δE`; so a minimiser can stall
/// once `|F|² < 4 λ_max δE`, that is below `2 √(λ_max δE)`. δE is the rounding of one evaluation:
/// a few dozen roundings inside each term and one per term in the sum, each of at most ε of the
/// sum of the terms' magnitudes `S` — `(n + 64) ε S` for `n` terms. λ_max is measured.
pub fn stall_floor(ff: &ForceField, at: &[[f64; 3]]) -> f64 {
    let (s, n) = term_magnitude(ff, at);
    let delta_e = (n as f64 + 64.0) * f64::EPSILON * s;
    let (_, largest) = hessian_extremes(ff, at);
    2.0 * (largest * delta_e).sqrt()
}

/// **A minimisation finished, checked independently of the minimiser's own report**: it says
/// converged or stalled, and the largest per-atom force, evaluated here from the force field, is at
/// most the tolerance — or, for a stall, at most the [`stall_floor`], the force below which no
/// correct minimiser can be asked to go. Returns that force, kcal mol⁻¹ Å⁻¹.
pub fn assert_minimised(
    ff: &ForceField,
    at: &[[f64; 3]],
    p: &pantometry_forcefield::Progress,
    tolerance: f64,
) -> f64 {
    use pantometry_forcefield::Status;
    let f = max_force(ff, at);
    let limit = match p.status {
        Status::Converged => tolerance,
        Status::Stalled => tolerance.max(stall_floor(ff, at)),
        Status::Running => panic!("still running after its step budget: {p:?}"),
    };
    assert!(f <= limit, "force {f:e} above {limit:e}: {p:?}");
    f
}

/// The relaxed minimum from `start`, unrestrained: energy (kcal/mol) and positions.
pub fn relaxed(ff: &ForceField, start: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
    use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
    use pantometry_forcefield::uff::KCAL_PER_MOL;
    let mut at = start.to_vec();
    let p = ff.minimise(&mut at, 50_000, SCAN_TOLERANCE * KCAL_PER_MOL_ANGSTROM);
    finished(&p, f64::NAN);
    (p.energy / KCAL_PER_MOL, at)
}

/// A point of a relaxed scan: the dihedral (degrees), the energy and the restraint's energy
/// (kcal/mol).
#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub phi: f64,
    pub energy: f64,
    pub restraint: f64,
}

/// A relaxed scan of the dihedral `atoms` from `start`: every 30° round the circle, each point
/// started from the one before, then the highest and lowest grid points refined by golden-section
/// search over the 60° about them (25 rounds: to 4e-4°). Returns the grid, the maximum and the
/// minimum.
pub fn scan(
    ff: &ForceField,
    c: &Component,
    start: &[[f64; 3]],
    atoms: [usize; 4],
) -> (Vec<Point>, Point, Point) {
    let (_, relaxed_start) = relaxed(ff, start);
    let first = pantometry_forcefield::minimise::dihedral(&relaxed_start, atoms) / DEG;
    let mut grid = Vec::new();
    let mut geoms = Vec::new();
    let mut at = relaxed_start;
    for k in 0..12 {
        let phi = first + 30.0 * k as f64;
        let (e, r, next) = held(ff, c, &at, atoms, phi);
        grid.push(Point {
            phi,
            energy: e,
            restraint: r,
        });
        geoms.push(next.clone());
        at = next;
    }
    let refine = |sign: f64| {
        let best = (0..grid.len())
            .max_by(|&a, &b| (sign * grid[a].energy).total_cmp(&(sign * grid[b].energy)))
            .expect("a grid");
        let from = &geoms[best];
        let f = |phi: f64| {
            let (e, r, _) = held(ff, c, from, atoms, phi);
            Point {
                phi,
                energy: e,
                restraint: r,
            }
        };
        let g = (5f64.sqrt() - 1.0) / 2.0;
        let (mut a, mut b) = (grid[best].phi - 30.0, grid[best].phi + 30.0);
        let mut x1 = b - g * (b - a);
        let mut x2 = a + g * (b - a);
        let mut f1 = f(x1);
        let mut f2 = f(x2);
        for _ in 0..25 {
            if sign * f1.energy > sign * f2.energy {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - g * (b - a);
                f1 = f(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + g * (b - a);
                f2 = f(x2);
            }
        }
        let best_point = if sign * f1.energy > sign * f2.energy {
            f1
        } else {
            f2
        };
        if sign * grid[best].energy > sign * best_point.energy {
            grid[best]
        } else {
            best_point
        }
    };
    let max = refine(1.0);
    let min = refine(-1.0);
    (grid, max, min)
}
