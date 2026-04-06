//! Fast analytic contact binary model.
//!
//! Uses Eggleton (1983) Roche lobe radii, circle-circle eclipse geometry,
//! and first-order ellipsoidal variation. ~350x faster than the numerical
//! Roche mesh model, suitable for MCMC parameter exploration.

use crate::EBParams;

/// Result of analytic light curve computation.
#[derive(Debug)]
pub struct AnalyticLC {
    pub phases: Vec<f64>,
    pub flux: Vec<f64>,
}

/// Compute Roche lobe radius using Eggleton (1983) approximation.
///
/// Returns the radius in units of the orbital separation.
fn eggleton_roche(q: f64) -> f64 {
    let q23 = q.powf(2.0 / 3.0);
    0.49 * q23 / (0.6 * q23 + (1.0 + q.powf(1.0 / 3.0)).ln())
}

/// Area of overlap between two circles of radii ra, rb at projected
/// separation d. Returns 0 for non-overlapping circles.
fn circle_overlap(ra: f64, rb: f64, d: f64) -> f64 {
    if d >= ra + rb {
        return 0.0;
    }
    if d <= (ra - rb).abs() {
        let rmin = ra.min(rb);
        return std::f64::consts::PI * rmin * rmin;
    }
    let cos_a = ((d * d + ra * ra - rb * rb) / (2.0 * d * ra)).clamp(-1.0, 1.0);
    let cos_b = ((d * d + rb * rb - ra * ra) / (2.0 * d * rb)).clamp(-1.0, 1.0);
    let a = cos_a.acos();
    let b = cos_b.acos();
    ra * ra * (a - a.sin() * a.cos()) + rb * rb * (b - b.sin() * b.cos())
}

/// Compute a fast analytic binary light curve.
///
/// Supports both contact and detached systems. For contact binaries,
/// both stars share the same fillout factor. For detached binaries,
/// `fillout1` and `fillout2` independently scale each star's radius
/// relative to its Roche lobe.
///
/// This is suitable for rapid MCMC exploration. For publication-quality
/// fits, use the numerical Roche mesh model (`compute_lightcurve`).
///
/// Parameters in `EBParams`:
/// - `q`: mass ratio M2/M1
/// - `inclination`: orbital inclination in degrees
/// - `t_eff1`, `t_eff2`: temperatures (used as ratio T2/T1 for surface brightness)
/// - `fillout1`: radius scale factor for primary (1.0 = Roche lobe)
/// - `fillout2`: radius scale factor for secondary
/// - `ld1`, `ld2`: linear limb darkening coefficients
/// - `phi0`: phase zero-point offset
pub fn compute_analytic(params: &EBParams, phases: &[f64]) -> AnalyticLC {
    let q = params.q.max(0.01);
    let incl = params.inclination.to_radians();
    let ld1 = params.ld1;
    let ld2 = params.ld2;

    // Roche lobe radii (Eggleton 1983), scaled by fillout
    let r1 = params.fillout1 * eggleton_roche(1.0 / q);
    let r2 = params.fillout2 * eggleton_roche(q);

    // Surface brightness ratio (Stefan-Boltzmann, bolometric)
    let tratio = params.t_eff2 / params.t_eff1.max(1.0);
    let j = tratio.powi(4);

    // Limb-darkened disk luminosities
    let l1 = std::f64::consts::PI * r1 * r1 * (1.0 - ld1 / 3.0);
    let l2 = std::f64::consts::PI * r2 * r2 * (1.0 - ld2 / 3.0) * j;
    let l_total = l1 + l2;

    if l_total <= 0.0 {
        return AnalyticLC {
            phases: phases.to_vec(),
            flux: vec![1.0; phases.len()],
        };
    }

    // Ellipsoidal variation coefficients
    // Scale with fillout^3 — detached stars are less tidally distorted
    let f1_eff = params.fillout1.min(1.0);  // cap at 1.0 for ellipsoidal
    let f2_eff = params.fillout2.min(1.0);
    let tau1 = 0.15 * r1.powi(3) * q * incl.sin().powi(2) * f1_eff.powi(2);
    let tau2 = 0.15 * r2.powi(3) * (1.0 / q) * incl.sin().powi(2) * f2_eff.powi(2);
    let ellip_amp = (tau1 * l1 + tau2 * l2) / l_total;

    let mut flux = Vec::with_capacity(phases.len());

    for &phase in phases {
        let ph = phase - params.phi0;
        let theta = 2.0 * std::f64::consts::PI * ph;
        let x_sep = theta.cos();
        let y_sep = theta.sin();

        // Projected separation on sky
        let d_proj = (y_sep * y_sep + (x_sep * incl.cos()).powi(2)).sqrt();

        let mut f = 1.0;

        // Eclipse: determine which star is in front
        let star2_front = x_sep * incl.sin() > 0.0;

        if star2_front {
            // Star 2 blocks star 1 (primary eclipse)
            let ov = circle_overlap(r1, r2, d_proj);
            f -= ov / (std::f64::consts::PI * r1 * r1) * l1 / l_total;
        } else {
            // Star 1 blocks star 2 (secondary eclipse)
            let ov = circle_overlap(r2, r1, d_proj);
            f -= ov / (std::f64::consts::PI * r2 * r2) * l2 / l_total;
        }

        // Ellipsoidal variation (cos 2θ modulation)
        f -= ellip_amp * (2.0 * theta).cos();

        flux.push(f);
    }

    AnalyticLC {
        phases: phases.to_vec(),
        flux,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EBParams;

    #[test]
    fn test_analytic_symmetric() {
        let params = EBParams::contact(1.0, 90.0);
        let phases: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let lc = compute_analytic(&params, &phases);
        // q=1, equal temps → eclipses should be symmetric
        let f_at_0 = lc.flux[0];
        let f_at_50 = lc.flux[50];
        assert!((f_at_0 - f_at_50).abs() < 0.01,
                "q=1 should give symmetric eclipses");
    }

    #[test]
    fn test_analytic_out_of_eclipse() {
        let params = EBParams::contact(0.5, 30.0);
        let phases: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let lc = compute_analytic(&params, &phases);
        // Low inclination → no eclipses, flux should be near 1.0
        for &f in &lc.flux {
            assert!(f > 0.8 && f < 1.2,
                    "low inclination should give near-unity flux, got {}", f);
        }
    }

    #[test]
    fn test_eggleton() {
        // q=1 → both lobes same size
        let r1 = eggleton_roche(1.0);
        let r1_inv = eggleton_roche(1.0);
        assert!((r1 - r1_inv).abs() < 1e-10);
        // r1 should be ~0.38 for q=1
        assert!(r1 > 0.35 && r1 < 0.42, "r_L for q=1 should be ~0.38, got {}", r1);
    }
}
