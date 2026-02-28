use lcurve::model::Model;
use lcurve::types::{Datum, Data};
use lcurve::orchestration::light_curve_comp;

/// Regression test: compare Rust light curve output against C++ reference data.
///
/// The C++ reference was generated with Tom Marsh's lroche on test_model.dat
/// with 500 evenly-spaced points from time -0.2 to 1.2.
#[test]
fn test_regression_vs_cpp() {
    let mdl = Model::from_file("test_data/test_model.dat")
        .expect("Failed to read test model");

    // Generate same time grid as C++ reference
    let ntime = 500;
    let time1 = -0.2;
    let time2 = 1.2;
    let expose = 0.001;

    let data: Data = (0..ntime)
        .map(|i| {
            let time = time1 + (time2 - time1) * i as f64 / (ntime - 1) as f64;
            Datum {
                time,
                expose,
                flux: 0.0,
                ferr: 0.0,
                weight: 1.0,
                ndiv: 1,
            }
        })
        .collect();

    let result = light_curve_comp(&mdl, &data, false, false)
        .expect("light_curve_comp failed");

    // Read C++ reference data
    let cpp_text = std::fs::read_to_string("test_data/test_lc.dat")
        .expect("Failed to read C++ reference data");

    let cpp_fluxes: Vec<f64> = cpp_text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            fields[2].parse::<f64>().expect("Failed to parse C++ flux")
        })
        .collect();

    assert_eq!(
        result.calc.len(),
        cpp_fluxes.len(),
        "Point count mismatch: Rust={} vs C++={}",
        result.calc.len(),
        cpp_fluxes.len()
    );

    let mut max_rel_err = 0.0f64;
    let mut max_rel_idx = 0;
    let mut sum_rel_err = 0.0;

    for (i, (&rust_flux, &cpp_flux)) in result.calc.iter().zip(cpp_fluxes.iter()).enumerate() {
        let rel_err = if cpp_flux != 0.0 {
            ((rust_flux - cpp_flux) / cpp_flux).abs()
        } else {
            (rust_flux - cpp_flux).abs()
        };

        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_rel_idx = i;
        }
        sum_rel_err += rel_err;
    }

    let mean_rel_err = sum_rel_err / cpp_fluxes.len() as f64;

    eprintln!("Regression test results:");
    eprintln!("  Max relative error: {:.6e} at point {}", max_rel_err, max_rel_idx);
    eprintln!("  Mean relative error: {:.6e}", mean_rel_err);

    // f32 intermediate precision gives ~1e-7 relative error at best;
    // we see ~4e-9 mean, ~4e-9 max. Use 1e-6 as a safe threshold.
    assert!(
        max_rel_err < 1e-6,
        "Max relative error {:.6e} exceeds threshold 1e-6 at point {}",
        max_rel_err, max_rel_idx
    );
}
