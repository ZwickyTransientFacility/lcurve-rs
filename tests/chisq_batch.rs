use lcurve::model::Model;
use lcurve::types::{Datum, Data};
use lcurve::orchestration::{light_curve_comp, chisq_batch};

/// Build synthetic data (no real observations — weight=1, ferr=1, flux=0).
fn make_data(ntime: usize) -> Data {
    let time1 = -0.2;
    let time2 = 1.2;
    (0..ntime)
        .map(|i| {
            let time = time1 + (time2 - time1) * i as f64 / (ntime - 1) as f64;
            Datum {
                time,
                expose: 0.001,
                flux: 0.0,
                ferr: 1.0,
                weight: 1.0,
                ndiv: 1,
            }
        })
        .collect()
}

/// chisq_batch with a single parameter set must match a sequential
/// light_curve_comp call with the same parameters.
#[test]
fn batch_single_matches_sequential() {
    let mdl = Model::from_file("test_data/test_model.dat")
        .expect("Failed to read test model");
    let data = make_data(100);

    let names = vec!["q", "iangle"];
    let values = vec![mdl.q.value, mdl.iangle.value];

    let batch = chisq_batch(&mdl, &data, &names, &values, false);
    let seq = light_curve_comp(&mdl, &data, false, true)
        .expect("light_curve_comp failed");

    assert_eq!(batch.len(), 1);
    assert!(
        (batch[0] - seq.chisq).abs() < 1e-10 * seq.chisq.abs().max(1.0),
        "batch={} vs sequential={}",
        batch[0],
        seq.chisq
    );
}

/// chisq_batch with N parameter sets must produce the same results
/// as N individual light_curve_comp calls.
#[test]
fn batch_matches_n_sequential() {
    let base = Model::from_file("test_data/test_model.dat")
        .expect("Failed to read test model");
    let data = make_data(100);

    let names = vec!["q", "iangle"];
    let q_base = base.q.value;
    let ia_base = base.iangle.value;

    // 5 parameter sets with small perturbations around the truth
    let offsets = [-0.02, -0.01, 0.0, 0.01, 0.02];
    let n = offsets.len();
    let mut flat_values = Vec::with_capacity(n * 2);
    for &dq in &offsets {
        flat_values.push(q_base + dq);
        flat_values.push(ia_base);
    }

    let batch = chisq_batch(&base, &data, &names, &flat_values, false);
    assert_eq!(batch.len(), n);

    for (i, &dq) in offsets.iter().enumerate() {
        let mut mdl = base.clone();
        mdl.set_param_value("q", q_base + dq);

        let seq = light_curve_comp(&mdl, &data, false, true)
            .expect("light_curve_comp failed");

        let tol = 1e-10 * seq.chisq.abs().max(1.0);
        assert!(
            (batch[i] - seq.chisq).abs() < tol,
            "set {}: batch={} vs sequential={} (diff={})",
            i,
            batch[i],
            seq.chisq,
            (batch[i] - seq.chisq).abs()
        );
    }
}

/// Unphysical parameters should produce NaN, not panic.
#[test]
fn batch_unphysical_gives_nan() {
    let base = Model::from_file("test_data/test_model.dat")
        .expect("Failed to read test model");
    let data = make_data(50);

    // r2 > Roche lobe should cause an error → NaN
    let names = vec!["r2"];
    let values = vec![0.99]; // very large — almost certainly larger than Roche lobe

    let batch = chisq_batch(&base, &data, &names, &values, false);
    assert_eq!(batch.len(), 1);
    assert!(
        batch[0].is_nan(),
        "Expected NaN for unphysical r2, got {}",
        batch[0]
    );
}

/// chisq_batch with scale=true should also match sequential.
#[test]
fn batch_with_scale_matches_sequential() {
    let base = Model::from_file("test_data/test_model.dat")
        .expect("Failed to read test model");
    let data = make_data(100);

    let names = vec!["q"];
    let values = vec![base.q.value];

    let batch = chisq_batch(&base, &data, &names, &values, true);
    let seq = light_curve_comp(&base, &data, true, true)
        .expect("light_curve_comp failed");

    assert_eq!(batch.len(), 1);
    assert!(
        (batch[0] - seq.chisq).abs() < 1e-10 * seq.chisq.abs().max(1.0),
        "scaled: batch={} vs sequential={}",
        batch[0],
        seq.chisq
    );
}
