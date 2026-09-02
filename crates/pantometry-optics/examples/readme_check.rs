// Compile-check the two documented examples exactly as written.
//
// The name is where they used to be. `README.md` was cut back to what this is and how to run it,
// and the snippets moved to `EXAMPLES.md` with the rest of the explaining; the name stays because
// it is in the gate, in CI and in two documents, and a rename buys nothing a comment does not.
use pantometry_optics::{fresnel_reflectance, Material, Spectrum, SurfaceFinish};
use pantometry_units::Length;

fn optics_example() {
    let bk7 = Material::from_catalog("N-BK7").unwrap();
    let n = bk7.index(Length::nm(587.56));
    let bare = fresnel_reflectance(1.0, n, 1.0);

    let coated = SurfaceFinish::broadband_ar().reflectance_at(1.0, n, 1.0, Length::nm(550.0));
    assert!(coated < bare / 10.0);

    let tungsten = Spectrum::blackbody(3200.0);
    assert!(tungsten.at(Length::nm(450.0)) < 0.45 * tungsten.at(Length::nm(650.0)));
}

fn coupling_example() {
    use pantometry_units::{
        Area, HeatCapacity, Irradiance, Mass, Power, SpecificHeat, Temperature,
    };

    let area: Area = Length::mm(10.0) * Length::mm(10.0);
    let absorbed: Power = Irradiance::mw_per_cm2(50.0) * area * 0.02;

    let capacity: HeatCapacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
    let rise: Temperature = (absorbed * pantometry_units::Time::s(1.0)) / capacity;
    assert!((rise.to_si() - 0.000_582_7).abs() < 1e-7, "{rise:?}");
}

fn main() {
    optics_example();
    coupling_example();
    println!("README examples ok");
}
