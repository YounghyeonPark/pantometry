//! The elastic network: springs between everything close enough, and the Hessian they make.
//!
//! # The model
//!
//! Every pair of nodes within a cutoff gets a spring of the same stiffness, at its rest length.
//! That is the entire physical content, and it is worth being exact about how little it is: no
//! chemistry, no bond types, no charges, no solvent, and every spring identical. The claim is
//! that a protein's **collective** motions are a property of the shape of the fold, and the
//! reason to believe it is not that the model is detailed — it is that the fluctuations it
//! predicts agree with the ones a crystallographer measured, which
//! [`Structure::experimental_fluctuations`] makes checkable.
//!
//! The potential is
//!
//! ```text
//! V = (γ/2) Σ (|r_ij| − |r⁰_ij|)²
//! ```
//!
//! over pairs inside the cutoff, and its Hessian at `r = r⁰` has the `3×3` blocks
//!
//! ```text
//! H_ij = −γ ê ê^T        H_ii = +γ Σ_j ê ê^T
//! ```
//!
//! with `ê` the unit vector from `i` to `j`.
//!
//! # What the two parameters do, and what they do not
//!
//! `γ` is one number for every spring, so it multiplies the whole Hessian and divides every
//! fluctuation. **It cannot change the shape of anything** — not which residues move, not the
//! ordering of the modes, not a correlation coefficient. It sets a scale and nothing else, which
//! is why a comparison against measured B-factors is a prediction with no free parameter in it
//! rather than a fit.
//!
//! The cutoff is not like that. Too small and the network falls into pieces, each of which can
//! drift away from the others for free — so the spectrum gains six zeros per piece and the
//! "lowest collective mode" is one half of the protein leaving. [`Network::components`] counts
//! them, because that failure produces numbers rather than an error.

use crate::spectrum::Symmetric;
use crate::structure::Structure;
use pantometry_units::{Length, Qty, Stiffness};

/// A structure with springs in it.
#[derive(Clone, Debug, PartialEq)]
pub struct Network {
    at: Vec<[f64; 3]>,
    residues: usize,
    cutoff: f64,
    gamma: f64,
}

impl Network {
    /// Springs between every pair of residues closer than `cutoff`.
    ///
    /// Fifteen ångström is the cutoff the anisotropic network model was published with and the
    /// one most of the literature uses; it is a parameter here rather than a constant because
    /// the right value depends on what is being asked, and because a structure with a loosely
    /// packed loop can fall apart at a value another structure is fine at — see
    /// [`Network::components`].
    ///
    /// # Panics
    ///
    /// If `cutoff` or `gamma` is not positive.
    pub fn new(structure: &Structure, cutoff: Length, gamma: Stiffness) -> Network {
        let cutoff = cutoff.to_si();
        let gamma = gamma.to_si();
        assert!(cutoff > 0.0, "a cutoff of {cutoff} m");
        assert!(gamma > 0.0, "a stiffness of {gamma} N/m");
        Network {
            at: structure.residues().iter().map(|r| r.at).collect(),
            residues: structure.len(),
            cutoff,
            gamma,
        }
    }

    /// The same network with a ligand's atoms added as nodes.
    ///
    /// A bound ligand is modelled the way everything else here is: as nodes joined by the same
    /// springs to whatever is within the cutoff. Nothing about the chemistry of the binding is
    /// in it, and nothing about *how tightly* it binds comes out — see the crate documentation
    /// for what this model does not answer.
    ///
    /// What does come out is a theorem rather than a guess, though not the obvious one.
    ///
    /// **The spectrum does not simply rise.** A ligand adds *nodes*, so the bound Hessian is
    /// `3(n+m)` square where the bare one was `3n`, and comparing them mode by mode is not
    /// Weyl's inequality and is not true — measured, mode 5 of a 24-residue helix falls from
    /// `0.204` to `0.178` N/m, because the ligand brings soft motions of its own that interleave
    /// with the protein's.
    ///
    /// **Every residue's fluctuation falls.** That is the statement about the protein's own
    /// coordinates, and it holds. Integrating the ligand out leaves the effective Hessian
    /// `H_pp − H_pl H_ll⁺ H_lp`, which is the bound potential minimised over where the ligand
    /// can go; every coupling term is a square, so that minimum is at least the bare potential,
    /// and both vanish at the reference structure. So the effective Hessian is `⪰ H_bare`, the
    /// covariance `k_BT H⁺` goes the other way, and **no residue can become more mobile on
    /// binding**. Binding stiffens; what is worth looking at is where, and how far from the site
    /// the effect reaches.
    ///
    /// # Panics
    ///
    /// If called twice — a network holds one ligand, and the second call would silently make the
    /// first ligand part of the protein.
    pub fn with_ligand(mut self, atoms: Vec<[f64; 3]>) -> Network {
        assert_eq!(
            self.at.len(),
            self.residues,
            "this network already has a ligand"
        );
        self.at.extend(atoms);
        self
    }

    /// How many nodes: residues plus any ligand atoms.
    pub fn nodes(&self) -> usize {
        self.at.len()
    }

    /// How many of the nodes are residues rather than ligand.
    pub fn residues(&self) -> usize {
        self.residues
    }

    /// Where node `i` is, in metres.
    ///
    /// # Panics
    ///
    /// If `i` is at or past [`Network::nodes`].
    pub fn at(&self, i: usize) -> [f64; 3] {
        self.at[i]
    }

    /// The cutoff.
    pub fn cutoff(&self) -> Length {
        Qty::from_si(self.cutoff)
    }

    /// The stiffness every spring has.
    pub fn stiffness(&self) -> Stiffness {
        Qty::from_si(self.gamma)
    }

    /// How many springs there are.
    pub fn springs(&self) -> usize {
        let mut n = 0;
        for i in 0..self.at.len() {
            for j in (i + 1)..self.at.len() {
                if self.distance(i, j) <= self.cutoff {
                    n += 1;
                }
            }
        }
        n
    }

    /// How many pieces the network falls into.
    ///
    /// One, for a cutoff that holds the structure together. More than one is the failure that
    /// does not announce itself: each piece contributes its own six rigid-body modes, so the
    /// spectrum's first "collective mode" is a piece drifting off, the fluctuations it predicts
    /// are enormous and shapeless, and every number downstream is a number.
    pub fn components(&self) -> usize {
        let n = self.at.len();
        let mut seen = vec![false; n];
        let mut pieces = 0;
        for start in 0..n {
            if seen[start] {
                continue;
            }
            pieces += 1;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(i) = stack.pop() {
                for (j, visited) in seen.iter_mut().enumerate() {
                    if !*visited && self.distance(i, j) <= self.cutoff {
                        *visited = true;
                        stack.push(j);
                    }
                }
            }
        }
        pieces
    }

    /// The `3n × 3n` Hessian, in N/m.
    ///
    /// The rigid-body motions are in its null space exactly, because `V` is a function of
    /// distances alone and a rigid motion changes none of them. That is not a property this
    /// assembles carefully and hopes for — it is what
    /// `a_hessian_annihilates_every_rigid_motion` measures.
    pub fn hessian(&self) -> Symmetric {
        let n = self.at.len();
        let mut h = Symmetric::zeros(3 * n);
        for i in 0..n {
            for j in (i + 1)..n {
                let d = self.distance(i, j);
                if d > self.cutoff {
                    continue;
                }
                let e = [
                    (self.at[j][0] - self.at[i][0]) / d,
                    (self.at[j][1] - self.at[i][1]) / d,
                    (self.at[j][2] - self.at[i][2]) / d,
                ];
                // The off-diagonal block covers nine positions and its transpose nine more, and
                // `add` writes one of each per call, so all nine `(a, b)` are visited.
                for a in 0..3 {
                    for b in 0..3 {
                        h.add(3 * i + a, 3 * j + b, -self.gamma * e[a] * e[b]);
                    }
                }
                // The diagonal blocks are visited for `a <= b` only. `add` writes both triangles,
                // so running the full nine here would add every off-diagonal entry twice — and
                // leave a matrix that is still symmetric, still positive semi-definite, and no
                // longer the Hessian of anything.
                for a in 0..3 {
                    for b in a..3 {
                        let g = self.gamma * e[a] * e[b];
                        h.add(3 * i + a, 3 * i + b, g);
                        h.add(3 * j + a, 3 * j + b, g);
                    }
                }
            }
        }
        h
    }

    /// The `3r × 3r` Hessian of the **protein's** coordinates, with the ligand integrated out.
    ///
    /// Equal to [`Network::hessian`] when there is no ligand. With one, it is the Schur
    /// complement `H_pp − H_pl H_ll⁺ H_lp`, and the reason it is a different matrix rather than
    /// a submatrix is what the bound system's statistics actually are.
    ///
    /// # Why the full Hessian is the wrong one to ask
    ///
    /// A residue's fluctuation is a property of the **marginal** distribution of the protein's
    /// coordinates:
    ///
    /// ```text
    /// ∫dy exp(−V(x,y)/kT)  ∝  exp(−W(x)/kT),    W(x) = min_y V(x,y)
    /// ```
    ///
    /// and the Hessian of `W` is this. The `pp` block of the full Hessian's *pseudo*-inverse is
    /// not that, because a pseudo-inverse removes the null space it is given — and the null space
    /// of the complex's Hessian is the complex's rigid motions, not the protein's.
    ///
    /// Measured, on a 24-residue helix with a three-atom ligand: through the full pseudo-inverse
    /// residue 4 comes out **1.003** times as mobile after binding, and through this it comes out
    /// **0.828**. The first of those violates the theorem in [`Network::with_ligand`], which is
    /// how the difference was found rather than argued about.
    pub fn protein_hessian(&self) -> Symmetric {
        let full = self.hessian();
        let (p, l) = (3 * self.residues, 3 * (self.at.len() - self.residues));
        if l == 0 {
            return full;
        }
        // `H_ll` is singular whenever the ligand has a direction nothing holds — one atom on one
        // spring has two — so this is a pseudo-inverse, and the directions it drops are the ones
        // the minimisation cannot use.
        let mut ll = Symmetric::zeros(l);
        for i in 0..l {
            for j in i..l {
                ll.set(i, j, full.get(p + i, p + j));
            }
        }
        let s = ll.eigen();
        let cut = f64::EPSILON.sqrt() * s.values().last().copied().unwrap_or(0.0).abs();
        let mut inverse = vec![0.0; l * l];
        for k in 0..l {
            if s.values()[k].abs() <= cut {
                continue;
            }
            let u = s.vector(k);
            for a in 0..l {
                for b in 0..l {
                    inverse[a * l + b] += u[a] * u[b] / s.values()[k];
                }
            }
        }
        // `H_pl C⁺`, once, rather than inside the `p²` loop below.
        let mut left = vec![0.0; p * l];
        for a in 0..p {
            for y in 0..l {
                left[a * l + y] = (0..l)
                    .map(|x| full.get(a, p + x) * inverse[x * l + y])
                    .sum();
            }
        }
        let mut effective = Symmetric::zeros(p);
        for a in 0..p {
            for b in a..p {
                let correction: f64 = (0..l).map(|y| left[a * l + y] * full.get(p + y, b)).sum();
                effective.set(a, b, full.get(a, b) - correction);
            }
        }
        effective
    }

    /// Between node `i` and node `j`, in metres.
    fn distance(&self, i: usize, j: usize) -> f64 {
        (0..3)
            .map(|k| (self.at[i][k] - self.at[j][k]).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}
