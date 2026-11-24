//! A step that carries a cell across its melting point leaves it between the coldest and the hottest
//! of what it exchanges with, as a step within one phase does.
//!
//! **It did not.** A block's limit charged each cell's row at the capacity the cell had, the slope of
//! its enthalpy, and a step that crosses the melting point moves the cell's heat at that slope and
//! gives it back at the other phase's. Where the phase it leaves holds more — water twice what ice
//! does — the cell lands past everything it exchanged with. Measured before this change, each step at
//! the block's own limit and each answered `Ok`:
//!
//! - one cell of water at 274.15 K under films to a room at 77 K came to **35.71 K**, 41 K below
//!   anything it touched; to a room at 20 K, −80.57 K; at 0 K, −121.37 K;
//! - water at 274 K between two cells of ice at 1 K came to −108.41 K;
//! - a solid holding twice what its liquid does, at 299 K between cells of its liquid at 1000 K,
//!   came to **1675 K**.
//!
//! The ones below absolute zero had been refused, since `f3a76e1`, by a pass over every cell of a
//! two-phase block after its step; the others were not asked.
//!
//! # What is checked against what
//!
//! **The limit against its closed form.** A step may carry a cell's enthalpy, `H`, no further than
//! that of the coldest thing it exchanges with, `dt·(−F) ≤ H(T) − H(T_lo)`, `F` what the step moves
//! into it — to the hottest, for a solid that holds more. Where everything the cell meets is at
//! `T_lo` and carries heat linearly, a film that convects, air or a face, that is its row over the
//! enthalpy's secant `(H(T) − H(T_lo))/(T − T_lo)` in place of its capacity: read as the ratio of
//! the limit to the same block's with that thing warm, which leaves the row as it was, or against
//! the stencil's own conductances where the geometry gives them. **Tight, not only safe**: a cell
//! stepped at the limit lands *on* that coldest thing, so a bound that charged more than it must —
//! the solid's capacity for every cell, or a radiating pair at its tangent — fails here.
//! **Refused at the slope's**, as a Fourier number by the ratio. **The same to the bit** where
//! nothing is past the cut-off, which for water is 116.70 K, and where what is past it cannot carry
//! the cell there; and only what a cell meets binds it, not what another cell does, nor a face that
//! carries nothing. And a rod of water stepped at its limit in a room at 20 K until it has frozen
//! through never goes below the room.

use pantometry_core::substance::{FusionProps, ThermalProps};
use pantometry_core::units::{
    Area, Density, LatentHeat, Length, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance, Violation};
use pantometry_thermal::{Environment, Face, Solid3D};

const FACES: [Face; 6] = [
    Face::XMin,
    Face::XMax,
    Face::YMin,
    Face::YMax,
    Face::ZMin,
    Face::ZMax,
];

/// What water holds, J/kg·K, twice what ice does.
const WATER_C: f64 = 4182.0;
/// What water conducts, a quarter of what ice does.
const WATER_K: f64 = 0.6;

/// Ice, with water as its liquid, so a block of it is two-phase; both radiating with `emissivity`,
/// so a cell part melted, which radiates as the mean of the two, does too.
fn water(emissivity: f64) -> Substance {
    let mut base = Substance::ice();
    if let Some(solid) = base.thermal.as_mut() {
        solid.emissivity = emissivity;
    }
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(WATER_K),
        specific_heat: SpecificHeat::j_per_kg_k(WATER_C),
        expansion: ThermalExpansion::ppm_per_k(207.0),
        emissivity,
    }))
}

/// A substance's enthalpy per kilogram, from its solid at its melting point, for a cell wholly one
/// phase or the other.
struct Curve {
    point: f64,
    latent: f64,
    solid: f64,
    liquid: f64,
}

impl Curve {
    fn of(s: &Substance) -> Curve {
        let fusion = s.fusion.expect("it melts");
        Curve {
            point: fusion.melting_point.to_si(),
            latent: fusion.latent_heat.to_si(),
            solid: s.thermal.expect("it conducts").specific_heat.to_si(),
            liquid: fusion
                .liquid
                .expect("it has a liquid")
                .specific_heat
                .to_si(),
        }
    }

    fn h(&self, t: f64) -> f64 {
        if t < self.point {
            self.solid * (t - self.point)
        } else {
            self.latent + self.liquid * (t - self.point)
        }
    }

    /// The secant from `from` to `to`, J/kg·K.
    fn secant(&self, from: f64, to: f64) -> f64 {
        (self.h(from) - self.h(to)) / (from - to)
    }
}

fn film(h: f64, room: f64) -> Environment {
    Environment {
        ambient: Temperature::from_si(room),
        convection_w_per_m2_k: h,
        area: Area::from_si(1e-6),
    }
}

fn limit(b: &Solid3D) -> f64 {
    b.max_stable_dt(Time::ZERO).to_si()
}

fn step(b: &mut Solid3D, dt: f64) -> Result<(), Violation> {
    b.step(Time::ZERO, Time::from_si(dt), &mut Exchange::new())
}

/// Every cell's temperature, in index order: a row along any axis.
fn row_of(b: &Solid3D) -> Vec<f64> {
    let (nx, ny, nz) = b.counts();
    let mut out = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                out.push(b.temperature_at(i, j, k).to_si());
            }
        }
    }
    out
}

/// Where a cell stepped at its limit lands, against what it lands on: within sixteen units in the
/// last place of `scale`. A cell landing on a room lands at `T_m + e/C_S`, `e` the enthalpy the step
/// left it, and its error is a few roundings of the enthalpy the step moved, over the capacity it
/// lands at — some hundreds of kelvin's worth, of the size `scale` is given. Measured, at most one
/// unit in the last place of 273.15 K, at a room at 0 K.
fn lands_on(landed: f64, on: f64, scale: f64) -> bool {
    (landed - on).abs() <= 16.0 * f64::EPSILON * scale
}

/// A step at `slope`, the limit the block had, is refused as a Fourier number, and nothing moves.
/// The ratio past its limit it reports is the limit's own, `slope · rate`, so this holds `step` to
/// the limit `max_stable_dt` reports; the ratio against its closed form is each caller's check.
fn refused_at_the_slope(why: &str, mut b: Solid3D, slope: f64, ratio: f64) {
    let bits = |b: &Solid3D| row_of(b).iter().map(|t| t.to_bits()).collect::<Vec<_>>();
    let before = bits(&b);
    let refusal = step(&mut b, slope).expect_err(why);
    let past = refusal.after / refusal.before;
    println!("  {why}: at the slope's limit, refused {past:.6}x past it, {ratio:.6}x expected");
    assert!(
        refusal.quantity == "Fourier number" && (past / ratio - 1.0).abs() < 1e-14,
        "{why}: {refusal}, {past} past, {ratio} expected"
    );
    assert_eq!(bits(&b), before, "{why}: a refused step moved the block");
}

/// **One cell of water under films lands on its room**, at a limit shorter than its slope's by
/// `C_L` over the secant to the room: ×1.2772 to a room at 0 K, ×1.1144 at 77 K, ×1.0023 at 116 K,
/// and at 120 K, past the cut-off of 116.70 K, the slope's to the bit. Films that only convect, so
/// a row is the same whatever the room, and the room changes only the capacity it is charged at.
#[test]
fn water_freezing_under_a_film_lands_on_its_room() {
    let curve = Curve::of(&water(0.0));
    let cut = curve.point - curve.latent / (curve.liquid - curve.solid);
    println!("water freezes past its slope below {cut} K");
    assert!((cut - 116.700_656_660_412_7).abs() < 1e-9, "{cut}");
    let start = 274.15;
    let one = |room: f64| {
        let mut b = Solid3D::new(
            "water",
            water(0.0),
            (1, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(start),
        );
        for face in FACES {
            b = b.losing_from(face, film(1e3, room));
        }
        b
    };
    let slope = limit(&one(300.0));
    for room in [120.0, 200.0] {
        assert_eq!(
            limit(&one(room)).to_bits(),
            slope.to_bits(),
            "a room at {room} K, past the cut-off, changed the limit"
        );
    }
    let mut worst: f64 = 0.0;
    for room in [0.0, 20.0, 77.0, 116.0] {
        let why = format!("water at {start} K to a room at {room} K");
        let dt = limit(&one(room));
        let ratio = slope / dt;
        let expected = curve.liquid / curve.secant(start, room);
        println!("  {why}: limit {dt:.9e} s, {ratio:.6}x shorter, {expected:.6}x by the secant");
        assert!(
            (ratio / expected - 1.0).abs() < 1e-14,
            "{why}: {ratio} shorter, {expected} by the secant"
        );
        let mut b = one(room);
        step(&mut b, dt).unwrap_or_else(|e| panic!("{why}: {e}"));
        let landed = b.temperature_at(0, 0, 0).to_si();
        worst = worst.max((landed - room).abs());
        println!(
            "    stepped at it: {landed:e} K, frozen {}",
            b.melted_fraction_at(0, 0, 0)
        );
        assert!(
            lands_on(landed, room, curve.point) && b.melted_fraction_at(0, 0, 0) == 0.0,
            "{why}: landed at {landed} K"
        );
        refused_at_the_slope(&why, one(room), slope, ratio);
    }
    println!("  furthest from its room: {worst:e} K");
    let mut b = one(120.0);
    step(&mut b, slope).expect("water to a room past the cut-off");
    assert!(b.temperature_at(0, 0, 0).to_si() > 120.0);
}

/// **A cell half melted lands on its room too.** At its melting point holding half its latent
/// heat, its capacity is the mean of the two phases' and its secant to the room
/// `(L/2 + C_S(T_m − T_lo))/(T_m − T_lo)`; the cut-off is the same, as it is for every fraction.
/// At the slope's limit, to a room at 77 K, it came to 56.36 K, and at 0 K to −60.68 K.
#[test]
fn a_half_melted_cell_lands_on_its_room() {
    let curve = Curve::of(&water(0.0));
    let one = |room: f64| {
        let mut b = Solid3D::new(
            "half melted",
            water(0.0),
            (1, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(curve.point),
        );
        for face in FACES {
            b = b.losing_from(face, film(1e3, room));
        }
        b.set_melted_fraction(0, 0, 0, 0.5);
        b
    };
    let slope = limit(&one(300.0));
    let mixed = 0.5 * (curve.solid + curve.liquid);
    for room in [0.0, 77.0] {
        let why = format!("a half-melted cell to a room at {room} K");
        let dt = limit(&one(room));
        let ratio = slope / dt;
        let below = curve.point - room;
        let expected = mixed / ((0.5 * curve.latent + curve.solid * below) / below);
        println!("  {why}: {ratio:.6}x shorter, {expected:.6}x by the secant");
        assert!(
            (ratio / expected - 1.0).abs() < 1e-14,
            "{why}: {ratio} shorter, {expected} by the secant"
        );
        let mut b = one(room);
        step(&mut b, dt).unwrap_or_else(|e| panic!("{why}: {e}"));
        let landed = b.temperature_at(0, 0, 0).to_si();
        println!("    stepped at it: {landed:e} K");
        assert!(
            lands_on(landed, room, curve.point) && b.melted_fraction_at(0, 0, 0) == 0.0,
            "{why}: landed at {landed} K"
        );
        refused_at_the_slope(&why, one(room), slope, ratio);
    }
}

/// **A mode released names each cell's phase**, as `set_temperature` does, so a step starts from a
/// state the limit can bound. `release_mode` wrote temperatures and left each cell the phase it
/// had: water at 300 K released to 50 K stayed liquid at 50 K, an enthalpy below its solid's at
/// absolute zero, and a step that moved no heat at all — one cell with nothing to meet — came out
/// at −19.37 K, its books unmoved. The pass after a two-phase step refused that; the limit cannot
/// see it, because nothing moves.
#[test]
fn a_mode_released_names_each_cells_phase() {
    let mut b = Solid3D::new(
        "released",
        water(0.0),
        (1, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(300.0),
    );
    assert_eq!(b.melted_fraction_at(0, 0, 0), 1.0);
    b.release_mode((0, 0, 0), Temperature::from_si(50.0), 0.0);
    assert_eq!(
        b.melted_fraction_at(0, 0, 0),
        0.0,
        "a cell released to 50 K is still liquid"
    );
    step(&mut b, 1.0).expect("a released cell with nothing to meet");
    let t = b.temperature_at(0, 0, 0).to_si();
    assert_eq!(t, 50.0, "a cell that met nothing moved, to {t} K");

    // Three cells about 100 K, one of them at 160 K and liquid again: each named by its own
    // temperature, and what a restore brings back is that.
    let mut row = Solid3D::new(
        "released",
        water(0.0),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(300.0),
    );
    row.release_mode((1, 0, 0), Temperature::from_si(250.0), 60.0);
    let named = |b: &Solid3D| {
        (0..3)
            .map(|i| b.melted_fraction_at(i, 0, 0))
            .collect::<Vec<_>>()
    };
    let temperatures = row_of(&row);
    let expected: Vec<f64> = temperatures
        .iter()
        .map(|t| if *t > 273.15 { 1.0 } else { 0.0 })
        .collect();
    println!(
        "released about 250 K: {temperatures:?}, melted {:?}",
        named(&row)
    );
    assert_eq!(named(&row), expected, "at {temperatures:?}");
    assert!(expected.contains(&1.0) && expected.contains(&0.0));
    let dt = limit(&row);
    step(&mut row, dt).expect("a released row at its limit");
    row.restore();
    assert_eq!(named(&row), expected, "a restore brought back other phases");
    assert_eq!(row_of(&row), temperatures);
}

/// **Water between two cells of ice at 1 K lands on the ice**, at the limit its own row sets once
/// charged at the secant: `1/max(k_f/(ρdx²c_S), 2k_f/(ρdx²S))`, `k_f` the faces' series mean. It
/// was the ice's, at which the water came to −108.41 K; refused there now by ×1.2506. Along each
/// axis, since a cell meets its neighbours through six faces read six ways.
#[test]
fn water_between_ice_lands_on_the_ice() {
    let ice = Substance::ice();
    let curve = Curve::of(&water(0.0));
    let rho = ice.density.to_si();
    let k_ice = ice.thermal.unwrap().conductivity.to_si();
    let dx: f64 = 1e-3;
    let k_f = 2.0 * k_ice * WATER_K / (k_ice + WATER_K);
    let per = rho * dx * dx;
    let at_ice = k_f / (per * curve.solid);
    let slope = 1.0 / at_ice.max(2.0 * k_f / (per * curve.liquid));
    let secant = curve.secant(274.0, 1.0);
    let expected = 1.0 / at_ice.max(2.0 * k_f / (per * secant));
    for (axis, counts, middle) in [
        ("x", (3, 1, 1), (1, 0, 0)),
        ("y", (1, 3, 1), (0, 1, 0)),
        ("z", (1, 1, 3), (0, 0, 1)),
    ] {
        let row = || {
            let mut b = Solid3D::new(
                "row",
                water(0.0),
                counts,
                Length::from_si(dx),
                Temperature::from_si(1.0),
            );
            b.set_temperature(middle.0, middle.1, middle.2, Temperature::from_si(274.0));
            b
        };
        let dt = limit(&row());
        println!(
            "water between ice along {axis}: limit {dt:e} s, {expected:e} s by the secant, \
             slope's {slope:e} s"
        );
        assert!(
            (dt / expected - 1.0).abs() < 1e-14,
            "along {axis}: {dt} s, {expected} s"
        );
        let mut b = row();
        step(&mut b, dt).expect("water between ice at its limit");
        let landed = row_of(&b);
        println!("  stepped at it: {landed:?}");
        assert!(
            lands_on(landed[1], 1.0, curve.point)
                && b.melted_fraction_at(middle.0, middle.1, middle.2) == 0.0,
            "along {axis}: the water landed at {} K",
            landed[1]
        );
        refused_at_the_slope("water between ice at 1 K", row(), slope, slope / dt);
    }
    assert!(
        (slope / expected - 1.2506).abs() < 1e-4,
        "{}",
        slope / expected
    );
}

/// **A solid that holds more than its liquid, melting, lands on its hottest neighbour**: the other
/// way across the melting point, `T_m + L/(C_S − C_L)` its cut-off. No catalogue entry does, and a
/// substance stated by its caller can: 4000 J/kg·K solid, 2000 liquid, 50 kJ/kg at 300 K, its
/// cut-off 325 K. At 299 K between cells of its liquid at 1000 K, stepped at the slope's limit of
/// 2 s, it came to 1675 K; its limit is the secant's, 1.0371 s, and it lands on 1000 K.
#[test]
fn a_solid_holding_more_than_its_liquid_melts_onto_its_hottest_neighbour() {
    let props = |c: f64| ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(1.0),
        specific_heat: SpecificHeat::j_per_kg_k(c),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: 0.0,
    };
    let odd = Substance::bulk("odd", Density::from_si(1000.0))
        .with_thermal(props(4000.0))
        .with_fusion(
            FusionProps::new(Temperature::from_si(300.0), LatentHeat::from_si(50_000.0))
                .with_liquid(props(2000.0)),
        );
    let curve = Curve::of(&odd);
    let row = || {
        let mut b = Solid3D::new(
            "row",
            odd.clone(),
            (3, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(1000.0),
        );
        b.set_temperature(1, 0, 0, Temperature::from_si(299.0));
        b
    };
    // `k_f = 1`, and `ρdx² = 1e-3`: the middle row is `2/(1e-3·S)`, each end's `1/(1e-3·c_L)`.
    let secant = curve.secant(1000.0, 299.0);
    let expected = 1.0 / (2.0 / (1e-3 * secant)).max(1.0 / (1e-3 * curve.liquid));
    let slope = 2.0;
    let dt = limit(&row());
    println!("a solid melting: limit {dt:e} s, {expected:e} s by the secant");
    assert!((dt / expected - 1.0).abs() < 1e-14, "{dt} s, {expected} s");
    let mut b = row();
    step(&mut b, dt).expect("a solid melting at its limit");
    let landed = row_of(&b);
    println!("  stepped at it: {landed:?}");
    assert!(
        lands_on(landed[1], 1000.0, 1000.0) && b.melted_fraction_at(1, 0, 0) == 1.0,
        "the solid landed at {} K",
        landed[1]
    );
    refused_at_the_slope("a solid melting toward 1000 K", row(), slope, slope / dt);

    // Three tenths melted at its melting point, the same way: a step carries `0.3·L` toward the
    // liquid's enthalpy at 1000 K, `F = 1.4 W`, and the slope is its capacity's, `C_S + 0.3·Δ`.
    let mushy = || {
        let mut b = row();
        b.set_temperature(1, 0, 0, Temperature::from_si(curve.point));
        b.set_melted_fraction(1, 0, 0, 0.3);
        b
    };
    let room = curve.h(1000.0) - 0.3 * curve.latent;
    let mixed = curve.solid + 0.3 * (curve.liquid - curve.solid);
    let expected = 1.0 / (2.0 * 700.0 / (1e-3 * room)).max(1.0 / (1e-3 * curve.liquid));
    let slope = 1e-3 * mixed / 2.0;
    let dt = limit(&mushy());
    println!("a solid three tenths melted: limit {dt:e} s, {expected:e} s by what it carries");
    assert!(
        (dt / expected - 1.0).abs() < 1e-14,
        "three tenths melted: {dt} s, {expected} s"
    );
    let mut b = mushy();
    step(&mut b, dt).expect("a solid three tenths melted at its limit");
    let landed = row_of(&b);
    println!("  stepped at it: {landed:?}");
    assert!(
        lands_on(landed[1], 1000.0, 1000.0) && b.melted_fraction_at(1, 0, 0) == 1.0,
        "three tenths melted, the solid landed at {} K",
        landed[1]
    );
    refused_at_the_slope(
        "three tenths melted toward 1000 K",
        mushy(),
        slope,
        slope / dt,
    );
}

/// **The air in a void and the far end of a clearance bind a cell as a room does**, by what they
/// carry and not by what the limit charges for them. Water facing air at 50 K, by convection alone,
/// lands on the air at a limit `C_L/S` shorter than with the air at 300 K. Across a clearance a
/// pair radiates, and is charged at its tangent at its hotter end, which from 274 K to 1 K is four
/// times what it carries: water across one from a heavy solid at 1 K cannot reach it in a step its
/// row allows, keeps its limit to the bit, and lands above it. A liquid that melts with almost
/// nothing can: 1 J/kg at 300 K, holding four times its solid, at 301 K across a clearance from the
/// heavy solid at 299 K, the pair carrying 0.99 of its tangent, came out of a step at its limit at
/// 296.1 K, and lands on 299 K.
#[test]
fn air_and_a_clearance_bind_a_cell_as_a_room_does() {
    let curve = Curve::of(&water(0.0));
    let start = 274.15;
    let facing = |air: f64| {
        Solid3D::new(
            "facing air",
            water(0.0),
            (2, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(start),
        )
        .empty(|i, _, _| i == 1)
        .air_in(Temperature::from_si(air), 1e3)
    };
    let warm = limit(&facing(300.0));
    let dt = limit(&facing(50.0));
    let expected = curve.liquid / curve.secant(start, 50.0);
    println!(
        "water facing air at 50 K: {:.6}x shorter, {expected:.6}x by the secant",
        warm / dt
    );
    assert!(
        ((warm / dt) / expected - 1.0).abs() < 1e-14,
        "water facing air at 50 K: {} shorter, {expected} by the secant",
        warm / dt
    );
    let mut b = facing(50.0);
    step(&mut b, dt).expect("water facing air at its limit");
    let landed = b.temperature_at(0, 0, 0).to_si();
    println!("  stepped at it: {landed} K");
    assert!(lands_on(landed, 50.0, curve.point), "landed at {landed} K");
    refused_at_the_slope("water facing air at 50 K", facing(50.0), warm, warm / dt);

    let heavy = Substance::bulk("heavy", Density::from_si(1e4)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(1.0),
        specific_heat: SpecificHeat::j_per_kg_k(1e4),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: 0.9,
    });
    let across = |liquid: Substance, from: f64, far: f64| {
        let mut b = Solid3D::new(
            "across a clearance",
            liquid,
            (3, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(from),
        )
        .empty(|i, _, _| i == 1);
        b.fill(heavy.clone(), |i, _, _| i == 2);
        b.set_temperature(2, 0, 0, Temperature::from_si(far));
        b
    };
    let dt = limit(&across(water(0.9), start, 1.0));
    assert_eq!(
        dt.to_bits(),
        limit(&across(water(0.9), start, 120.0)).to_bits(),
        "water radiating to 1 K was given a limit its row does not"
    );
    let mut b = across(water(0.9), start, 1.0);
    step(&mut b, dt).expect("water across a clearance at its limit");
    let landed = row_of(&b);
    println!("water across a clearance from 1 K, stepped at its limit: {landed:?}");
    assert!(landed[0] > 1.0 && landed[0] < start, "{landed:?}");

    let props = |c: f64| ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(1.0),
        specific_heat: SpecificHeat::j_per_kg_k(c),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: 0.9,
    };
    let slight = Substance::bulk("slight", Density::from_si(1000.0))
        .with_thermal(props(1000.0))
        .with_fusion(
            FusionProps::new(Temperature::from_si(300.0), LatentHeat::from_si(1.0))
                .with_liquid(props(4000.0)),
        );
    let dt = limit(&across(slight.clone(), 301.0, 299.0));
    let mut b = across(slight.clone(), 301.0, 299.0);
    step(&mut b, dt).expect("a slight liquid across a clearance at its limit");
    let landed = row_of(&b);
    println!("a slight liquid across a clearance from 299 K, stepped at its limit: {landed:?}");
    assert!(
        lands_on(landed[0], 299.0, 300.0) && b.melted_fraction_at(0, 0, 0) == 0.0,
        "the liquid landed at {} K",
        landed[0]
    );

    // The same in a block built of water, whose bound of 116.70 K nothing here is past: the
    // slight liquid's, 299.9997 K, is the block's, the highest of its materials'.
    let mut two = Solid3D::new(
        "two materials",
        water(0.9),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(301.0),
    )
    .empty(|i, _, _| i == 1);
    two.fill(slight, |i, _, _| i == 0);
    two.fill(heavy.clone(), |i, _, _| i == 2);
    two.set_temperature(2, 0, 0, Temperature::from_si(299.0));
    let dt = limit(&two);
    step(&mut two, dt).expect("a slight liquid in a block of water at its limit");
    let landed = row_of(&two);
    println!("the same in a block built of water: {landed:?}");
    assert!(
        lands_on(landed[0], 299.0, 300.0),
        "in a block of water, the liquid landed at {} K",
        landed[0]
    );

    // Air without limit holds each surface of a clearance at its own temperature, so the far end
    // carries nothing to the water and the air is what it meets: water across from the heavy solid
    // at 10 K, in air at 50 K, lands on 50 K.
    let mut held = Solid3D::new(
        "held by its air",
        water(0.9),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(start),
    )
    .empty(|i, _, _| i == 1)
    .air_in(Temperature::from_si(50.0), f64::INFINITY);
    held.fill(heavy, |i, _, _| i == 2);
    held.set_temperature(2, 0, 0, Temperature::from_si(10.0));
    let dt = limit(&held);
    step(&mut held, dt).expect("water in air without limit at its limit");
    let landed = row_of(&held);
    println!("water in air at 50 K without limit, across from 10 K: {landed:?}");
    assert!(
        lands_on(landed[0], 50.0, curve.point),
        "held by its air, the water landed at {} K",
        landed[0]
    );
}

/// **Only what a cell meets binds it**, and the coldest thing it meets is how far a step may carry
/// it, so anything colder that it does not meet would let it go further. Water beside a heavy
/// solid at 50 K, whose far face has a film to a room at 0 K, lands on the solid, not the room. And
/// a face that carries nothing — no film and no radiation — is not met: water under a film to a
/// room at 100 K, with another face to a room at 0 K that it cannot reach, lands on 100 K, at the
/// limit it has with that face's room at 300 K, to the bit.
#[test]
fn only_what_a_cell_meets_binds_it() {
    let start = 274.15;
    let heavy = Substance::bulk("heavy", Density::from_si(1e4)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(100.0),
        specific_heat: SpecificHeat::j_per_kg_k(1e4),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: 0.0,
    });
    let mut beside = Solid3D::new(
        "beside",
        water(0.0),
        (2, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(start),
    )
    .losing_from(Face::XMin, film(10.0, 0.0));
    beside.fill(heavy, |i, _, _| i == 0);
    beside.set_temperature(0, 0, 0, Temperature::from_si(50.0));
    let dt = limit(&beside);
    step(&mut beside, dt).expect("water beside a cold solid at its limit");
    let landed = row_of(&beside);
    println!("water beside a solid at 50 K, stepped at its limit: {landed:?}");
    assert!(
        lands_on(landed[1], 50.0, start),
        "the water landed at {} K",
        landed[1]
    );

    let sealed = |room: f64| {
        Solid3D::new(
            "sealed",
            water(0.0),
            (1, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(start),
        )
        .losing_from(Face::XMin, film(1e3, 100.0))
        .losing_from(Face::XMax, film(0.0, room))
    };
    let dt = limit(&sealed(0.0));
    assert_eq!(
        dt.to_bits(),
        limit(&sealed(300.0)).to_bits(),
        "a face that carries nothing bound the cell"
    );
    let mut b = sealed(0.0);
    step(&mut b, dt).expect("water under a film at its limit");
    let landed = b.temperature_at(0, 0, 0).to_si();
    println!("water under a film to 100 K and a face to 0 K it cannot reach: {landed} K");
    assert!(
        lands_on(landed, 100.0, start),
        "the water landed at {landed} K"
    );

    // Nor is the far end of a clearance it cannot reach: a perfect insulator at 1 K, which the
    // pair's half cell carries nothing through, across one from water under a film to 100 K.
    let insulator =
        Substance::bulk("insulator", Density::from_si(1000.0)).with_thermal(ThermalProps {
            conductivity: ThermalConductivity::w_per_m_k(0.0),
            specific_heat: SpecificHeat::j_per_kg_k(1000.0),
            expansion: ThermalExpansion::ppm_per_k(10.0),
            emissivity: 0.9,
        });
    let mut across = Solid3D::new(
        "across",
        water(0.9),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(start),
    )
    .empty(|i, _, _| i == 1)
    .losing_from(Face::XMin, film(1e3, 100.0));
    across.fill(insulator, |i, _, _| i == 2);
    across.set_temperature(2, 0, 0, Temperature::from_si(1.0));
    let dt = limit(&across);
    step(&mut across, dt).expect("water across from an insulator at its limit");
    let landed = across.temperature_at(0, 0, 0).to_si();
    println!(
        "water under a film to 100 K, across a clearance from an insulator at 1 K: {landed} K"
    );
    assert!(
        lands_on(landed, 100.0, start),
        "the water across from an insulator landed at {landed} K"
    );
}

/// **A rod of water freezing in a room at 20 K never goes below it.** Twelve cells in a row from
/// 274.15 K to 634.15 K, every face under a film of 1e5 W/m²·K that also radiates, stepped at its
/// limit, asked again each step, until every cell has frozen and gone below 116.70 K: every cell
/// after every step between the room and the warmest it started, at least one cell crossing from
/// liquid to solid in one step, which is the step the limit was wrong for, and the coldest landing
/// on the room, where the limit puts the cell that sets it.
///
/// A rod and a film that strong because a cell undershoots its coldest input only when that input
/// is most of what it meets: four fifths of its row and more, from 274 K to a room at 20 K. An end
/// cell of a rod is nine tenths room, and of six cells from 274.15 K at the slope's limit came to
/// −31.98 K; a cell of a thicker block under a weaker film, a third room, lands past the mean of
/// what it meets and not past the room, and a test of it passed before this change. A film that
/// strong freezes the rod through in three steps.
#[test]
fn a_rod_freezing_in_a_cryogenic_room_stays_above_it() {
    let (nx, ny, nz) = (12, 1, 1);
    let dx: f64 = 1e-3;
    let room = 20.0;
    let mut b = Solid3D::new(
        "freezing",
        water(0.96),
        (nx, ny, nz),
        Length::from_si(dx),
        Temperature::from_si(274.15),
    );
    for face in FACES {
        let (a, c) = match face {
            Face::XMin | Face::XMax => (ny, nz),
            Face::YMin | Face::YMax => (nx, nz),
            Face::ZMin | Face::ZMax => (nx, ny),
        };
        b = b.losing_from(
            face,
            Environment {
                ambient: Temperature::from_si(room),
                convection_w_per_m2_k: 1e5,
                area: Area::from_si((a * c) as f64 * dx * dx),
            },
        );
    }
    let mut warmest: f64 = 0.0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let t = 274.15 + 60.0 * ((i + 2 * j + 3 * k) % 7) as f64;
                b.set_temperature(i, j, k, Temperature::from_si(t));
                warmest = warmest.max(t);
            }
        }
    }
    let cells = |b: &Solid3D| {
        let mut out = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    out.push((
                        b.temperature_at(i, j, k).to_si(),
                        b.melted_fraction_at(i, j, k),
                    ));
                }
            }
        }
        out
    };
    let (mut steps, mut crossings, mut coldest) = (0, 0, f64::INFINITY);
    loop {
        let before = cells(&b);
        let dt = limit(&b);
        step(&mut b, dt).unwrap_or_else(|e| panic!("step {steps}: {e}"));
        steps += 1;
        let after = cells(&b);
        for (n, ((_, was), (t, now))) in before.iter().zip(&after).enumerate() {
            coldest = coldest.min(*t);
            assert!(
                (*t >= room || lands_on(*t, room, 273.15)) && *t <= warmest,
                "step {steps}: cell {n} at {t} K, outside {room} K to {warmest} K"
            );
            if *was == 1.0 && *now == 0.0 {
                crossings += 1;
            }
        }
        if after.iter().all(|(t, phi)| *phi == 0.0 && *t < 116.0) || steps == 100_000 {
            break;
        }
    }
    println!(
        "frozen through in {steps} steps, {crossings} cells liquid to solid in one step, \
         coldest {coldest} K"
    );
    assert!(steps < 100_000, "the block did not freeze through");
    assert!(
        crossings > 0,
        "no cell crossed in one step, so the bound was never at stake"
    );
    assert!(
        lands_on(coldest, room, 273.15),
        "the coldest cell came to {coldest} K, not the room"
    );
}
