//! A countable number of things at places — the other half of what a domain can be.
//!
//! [`ScalarField`](crate::ScalarField) covers the domains that are continua: a room, a bar, a
//! solid. A caller can sample one anywhere without knowing what it is, which is what let the
//! renderer stop naming `Room` and `Bar1D`.
//!
//! Then the other half stayed uncovered. An orbit, a bouncing ball and a box of atoms are *not*
//! fields — they are a finite number of bodies at positions, and rasterising them would invent a
//! continuum they do not have. `as_field` returning `None` for them is the honest answer and it
//! left every caller back at `domain_as::<NBody>`, `domain_as::<ContactSystem>`,
//! `domain_as::<Fluid>`. That is `FRICTION.md` finding 11, recorded and not fixed for months.
//!
//! Splitting the layers apart is what forced it. A scene layer that had to name three domains to
//! find out where anything *was* would need editing every time a physics arrived, and that is
//! the one thing this structure exists to prevent.
//!
//! # What belongs here and what does not
//!
//! A body's **position** is physics. So is a **real wall** it cannot cross — a periodic cell is
//! a boundary condition, and an atom leaving one face genuinely enters the opposite one.
//!
//! The **extent to draw** is not. An orbit's box is only as big as the picture wants; there is
//! nothing physical at its edge. So [`Bodies::cell`] reports a wall or `None`, and a view that
//! wants framing computes it from the positions it has — over the whole run, so that a body
//! moving is a body moving rather than the frame rescaling underneath it.

use pantometry_units::LengthVec;

/// A domain that is a finite number of bodies rather than a continuum.
///
/// Indexed rather than iterated so a caller can take one body without building a vector, and
/// because the order is meaningful: index is identity here, and body 3 in one frame is body 3 in
/// the next.
pub trait Bodies {
    /// How many bodies there are. Constant through a run for every domain that implements this.
    fn count(&self) -> usize;

    /// Where body `i` is.
    ///
    /// # Panics
    ///
    /// May panic for `i >= count()`, like any indexing. A caller iterating `0..count()` is fine.
    fn position(&self, i: usize) -> LengthVec;

    /// What body `i` is worth colouring by — a speed, a height, a charge.
    ///
    /// One scalar, chosen by the domain, because only it knows which of its per-body quantities a
    /// reader wants. A view showing this on one scale across a whole run is showing the physics;
    /// one rescaling per frame is showing the scale.
    fn value(&self, i: usize) -> f64;

    /// What [`value`](Bodies::value) is in.
    fn value_unit(&self) -> &'static str;

    /// The wall these bodies live in, if it is a real one.
    ///
    /// `Some` for a periodic cell, which is a boundary condition: an atom leaving one face
    /// enters the opposite one, and drawing that box is drawing physics. `None` for an orbit or
    /// a falling ball, whose extent is a property of the picture and not of the problem — a view
    /// that wants to frame those should measure the positions rather than be told.
    fn cell(&self) -> Option<(LengthVec, LengthVec)> {
        None
    }

    /// What body `i` **is**, for a reader matching it against something this run did not produce.
    ///
    /// **Index is identity within a run and says nothing outside it.** A file holding forty-six
    /// points and forty-six numbers is a point cloud: nothing in it can be matched against a
    /// crystallographer's temperature factors, against a second structure, or against the entry it
    /// was read from. `pantometry-protein` knows each residue's chain, number and name and had
    /// nowhere to put them, so a protein arrived in a run as arithmetic about anonymous points —
    /// which is a picture and not a measurement.
    ///
    /// `None` for a domain whose bodies have no name outside their order, which is most of them:
    /// the third planet in an orbit is the third planet.
    fn label(&self, _i: usize) -> Option<String> {
        None
    }

    /// Which bodies are bonded, as index pairs into `0..count()`.
    ///
    /// **Pairs and not runs**, because a molecule's bonds are not one chain: a ligand branches, a
    /// protein has several chains and a disulphide joins two residues that are nowhere near each
    /// other in sequence. A backbone is then 45 pairs for a 46-residue protein, which is nothing.
    ///
    /// A domain that returns these is promising they are **real** — that the two bodies are joined
    /// in the thing being modelled. A chain read out of a file breaks where its residue numbering
    /// does, and bonding across that gap draws a bond nobody measured.
    ///
    /// Empty for bodies that are not joined to each other, which is most of them.
    fn bonds(&self) -> Vec<[u32; 2]> {
        Vec::new()
    }
}
