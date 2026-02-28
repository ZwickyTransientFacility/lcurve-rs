use lcurve::model::Model;
use lcurve::grid;
use lcurve::brightness;
use lcurve_roche::Star;
use lcurve_subs::planck;

fn main() {
    let mdl = Model::from_file("test_data/test_model.dat").unwrap();

    println!("=== Model Parameters ===");
    println!("q={}, iangle={}, r1={}, r2={}", mdl.q.value, mdl.iangle.value, mdl.r1.value, mdl.r2.value);
    println!("t1={}, t2={}, wavelength={}", mdl.t1.value, mdl.t2.value, mdl.wavelength);
    println!("roche1={}, roche2={}", mdl.roche1, mdl.roche2);
    println!("nlat1f={}, nlat2f={}", mdl.nlat1f, mdl.nlat2f);
    println!("npole={}, gravity_dark1={}, gravity_dark2={}", mdl.npole, mdl.gravity_dark1.value, mdl.gravity_dark2.value);
    println!("gdark_bolom1={}, gdark_bolom2={}", mdl.gdark_bolom1, mdl.gdark_bolom2);
    println!("absorb={}", mdl.absorb.value);

    let (r1, r2) = mdl.get_r1r2();
    println!("get_r1r2 -> r1={}, r2={}", r1, r2);

    let rl1 = lcurve_roche::lagrange::xl1(mdl.q.value).unwrap();
    let xl12v = lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value).unwrap();
    println!("xl1={}, 1-xl1={}", rl1, 1.0 - rl1);
    println!("xl12={}, 1-xl12={}", xl12v, 1.0 - xl12v);

    println!("\n=== Planck ===");
    println!("planck(550, 15000) = {:.15e}", planck(550.0, 15000.0));
    println!("planck(550, 3500)  = {:.15e}", planck(550.0, 3500.0));

    println!("\n=== Star 1 Fine Grid ===");
    let mut star1f = grid::set_star_grid(&mdl, Star::Primary, true).unwrap();
    println!("  N = {}", star1f.len());
    let a1: f64 = star1f.iter().map(|p| p.area as f64).sum();
    println!("  total area = {:.15e}", a1);
    println!("  4*pi*r1^2  = {:.15e}", 4.0 * std::f64::consts::PI * r1 * r1);

    println!("\n=== Star 2 Fine Grid ===");
    let mut star2f = grid::set_star_grid(&mdl, Star::Secondary, true).unwrap();
    println!("  N = {}", star2f.len());
    let a2: f64 = star2f.iter().map(|p| p.area as f64).sum();
    println!("  total area = {:.15e}", a2);

    println!("\n=== Setting continuum ===");
    brightness::set_star_continuum(&mdl, &mut star1f, &mut star2f).unwrap();

    let f1: f64 = star1f.iter().map(|p| p.flux as f64).sum();
    let f2: f64 = star2f.iter().map(|p| p.flux as f64).sum();
    println!("  star1 total B*A = {:.15e}", f1);
    println!("  star2 total B*A = {:.15e}", f2);
    println!("  sum             = {:.15e}", f1 + f2);

    // Check star1 flux: should be planck(550,15000) * 4*pi*r1^2
    let expected_f1 = planck(550.0, 15000.0) * 4.0 * std::f64::consts::PI * r1 * r1;
    println!("  star1 expected (planck*4piR^2) = {:.15e}", expected_f1);

    // Check a few elements
    println!("\n=== Star1 elem[0] ===");
    let p = &star1f[0];
    println!("  posn=({:.10e},{:.10e},{:.10e})", p.posn.x, p.posn.y, p.posn.z);
    println!("  dirn=({:.10e},{:.10e},{:.10e})", p.dirn.x, p.dirn.y, p.dirn.z);
    println!("  area={:.10e} grav={:.10e} flux={:.10e}", p.area, p.gravity, p.flux);

    println!("\n=== Star2 elem[0] ===");
    let p = &star2f[0];
    println!("  posn=({:.10e},{:.10e},{:.10e})", p.posn.x, p.posn.y, p.posn.z);
    println!("  dirn=({:.10e},{:.10e},{:.10e})", p.dirn.x, p.dirn.y, p.dirn.z);
    println!("  area={:.10e} grav={:.10e} flux={:.10e}", p.area, p.gravity, p.flux);
}
