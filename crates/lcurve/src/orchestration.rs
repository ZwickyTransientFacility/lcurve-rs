use lcurve_subs::{sqr, TWOPI, PI, DAY, C};
use lcurve_roche::Star;
use rayon::prelude::*;
use crate::model::Model;
use crate::types::{Data, Ginterp, Point};
use crate::LcurveError;

/// Result of a light curve computation.
pub struct LcResult {
    /// Computed flux values
    pub calc: Vec<f64>,
    /// White dwarf contribution at phase 0.5
    pub wdwarf: f64,
    /// Weighted chi-squared
    pub chisq: f64,
    /// Weighted number of OK points
    pub wnok: f64,
    /// Flux-weighted log10(g) for star 1 (CGS)
    pub logg1: f64,
    /// Flux-weighted log10(g) for star 2 (CGS)
    pub logg2: f64,
    /// Volume-averaged radius of star 1
    pub rv1: f64,
    /// Volume-averaged radius of star 2
    pub rv2: f64,
    /// Scale factors [star1, disc, edge, spot, star2]
    pub sfac: Vec<f64>,
}

/// Compute a complete light curve from a model and data.
///
/// * `mdl` - the model parameters
/// * `data` - the data defining the times
/// * `scale` - compute scale factors by SVD or least squares
/// * `rdata` - true if real data (compute chi-squared)
pub fn light_curve_comp(
    mdl: &Model, data: &Data, scale: bool, rdata: bool,
) -> Result<LcResult, LcurveError> {
    let n = data.len();
    let mut calc = vec![0.0; n];

    let (mut r1, mut r2) = mdl.get_r1r2();
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value)?;
    if r2 < 0.0 { r2 = rl2; } else if r2 > rl2 {
        return Err(LcurveError::Generic("light_curve_comp: secondary larger than Roche lobe".into()));
    }

    let ldc1 = mdl.get_ldc1();
    let ldc2 = mdl.get_ldc2();

    // Gravitational lensing radius
    let rlens1 = if mdl.glens1 {
        let gm = (1000.0 * mdl.velocity_scale.value).powi(3) * mdl.tperiod * DAY / TWOPI;
        let a = (gm / sqr(TWOPI / DAY / mdl.tperiod)).powf(1.0 / 3.0);
        4.0 * gm / (1.0 + mdl.q.value) / a / sqr(C)
    } else {
        0.0
    };

    // Generate star grids — fine
    let mut star1f = crate::grid::set_star_grid(mdl, Star::Primary, true)?;
    let mut star2f = crate::grid::set_star_grid(mdl, Star::Secondary, true)?;
    crate::brightness::set_star_continuum(mdl, &mut star1f, &mut star2f)?;

    // Coarse grids
    let mut star1c = if mdl.nlat1f == mdl.nlat1c {
        star1f.clone()
    } else {
        crate::grid::set_star_grid(mdl, Star::Primary, false)?
    };

    let copy2 = mdl.nlat2f == mdl.nlat2c
        && (!mdl.npole || r1 >= r2 || (mdl.nlatfill == 0 && mdl.nlngfill == 0));

    let mut star2c = if copy2 {
        star2f.clone()
    } else {
        crate::grid::set_star_grid(mdl, Star::Secondary, false)?
    };

    if mdl.nlat1c != mdl.nlat1f || !copy2 {
        crate::brightness::set_star_continuum(mdl, &mut star1c, &mut star2c)?;
    }

    // Grid interpolation
    let mut gint = Ginterp {
        phase1: mdl.phase1,
        phase2: mdl.phase2,
        scale11: 1.0, scale12: 1.0,
        scale21: 1.0, scale22: 1.0,
    };

    if mdl.nlat1c != mdl.nlat1f {
        let ff = crate::flux::comp_star1(mdl.iangle.value, &ldc1, 0.9999999999 * mdl.phase1,
                                          0.0, 1, mdl.q.value, mdl.beam_factor1.value,
                                          mdl.velocity_scale.value as f32, &gint, &star1f, &star1c);
        let fc = crate::flux::comp_star1(mdl.iangle.value, &ldc1, 1.0000000001 * mdl.phase1,
                                          0.0, 1, mdl.q.value, mdl.beam_factor1.value,
                                          mdl.velocity_scale.value as f32, &gint, &star1f, &star1c);
        gint.scale11 = ff / fc;

        let ff = crate::flux::comp_star1(mdl.iangle.value, &ldc1, 1.0 - 0.9999999999 * mdl.phase1,
                                          0.0, 1, mdl.q.value, mdl.beam_factor1.value,
                                          mdl.velocity_scale.value as f32, &gint, &star1f, &star1c);
        let fc = crate::flux::comp_star1(mdl.iangle.value, &ldc1, 1.0 - 1.0000000001 * mdl.phase1,
                                          0.0, 1, mdl.q.value, mdl.beam_factor1.value,
                                          mdl.velocity_scale.value as f32, &gint, &star1f, &star1c);
        gint.scale12 = ff / fc;
    }

    if !copy2 {
        let ff = crate::flux::comp_star2(mdl.iangle.value, &ldc2, 1.0 - 1.0000000001 * mdl.phase2,
                                          0.0, 1, mdl.q.value, mdl.beam_factor2.value,
                                          mdl.velocity_scale.value as f32, mdl.glens1, rlens1,
                                          &gint, &star2f, &star2c);
        let fc = crate::flux::comp_star2(mdl.iangle.value, &ldc2, 1.0 - 0.9999999999 * mdl.phase2,
                                          0.0, 1, mdl.q.value, mdl.beam_factor2.value,
                                          mdl.velocity_scale.value as f32, mdl.glens1, rlens1,
                                          &gint, &star2f, &star2c);
        gint.scale21 = ff / fc;

        let ff = crate::flux::comp_star2(mdl.iangle.value, &ldc2, 1.0000000001 * mdl.phase2,
                                          0.0, 1, mdl.q.value, mdl.beam_factor2.value,
                                          mdl.velocity_scale.value as f32, mdl.glens1, rlens1,
                                          &gint, &star2f, &star2c);
        let fc = crate::flux::comp_star2(mdl.iangle.value, &ldc2, 0.9999999999 * mdl.phase2,
                                          0.0, 1, mdl.q.value, mdl.beam_factor2.value,
                                          mdl.velocity_scale.value as f32, mdl.glens1, rlens1,
                                          &gint, &star2f, &star2c);
        gint.scale22 = ff / fc;
    }

    // Disc and spot
    let mut disc: Vec<Point> = Vec::new();
    let mut edge: Vec<Point> = Vec::new();
    let mut spot: Vec<Point> = Vec::new();

    if mdl.add_disc {
        disc = crate::grid::set_disc_grid(mdl)?;
        edge = crate::grid::set_disc_edge(mdl, true)?;

        let rdisc1 = if mdl.rdisc1.value > 0.0 { mdl.rdisc1.value } else { r1 };
        let rdisc2 = if mdl.rdisc2.value > 0.0 { mdl.rdisc2.value } else { mdl.radius_spot.value };

        if mdl.opaque {
            // Apply disc eclipse to star grids
            for grids in [&mut star1f, &mut star1c, &mut star2f, &mut star2c] {
                grids.par_iter_mut().for_each(|pt| {
                    let eclipses = lcurve_roche::disc_eclipse::disc_eclipse(
                        mdl.iangle.value, rdisc1, rdisc2,
                        mdl.beta_disc.value, mdl.height_disc.value, &pt.posn);
                    for e in eclipses {
                        pt.eclipse.push(e);
                    }
                });
            }
        }

        crate::brightness::set_disc_continuum(
            rdisc2, mdl.temp_disc.value, mdl.texp_disc.value,
            mdl.wavelength, &mut disc);

        crate::brightness::set_edge_continuum(
            mdl.temp_edge.value, r2, mdl.t2.value.abs(),
            mdl.absorb_edge.value, mdl.wavelength, &mut edge);
    }

    if mdl.add_spot {
        spot = crate::grid::set_bright_spot_grid(mdl)?;
    }

    // Polynomial fudge factors
    let (mut xmin, mut xmax) = (data[0].time, data[0].time);
    for d in data.iter().skip(1) {
        if d.time < xmin { xmin = d.time; }
        if d.time > xmax { xmax = d.time; }
    }
    let middle = (xmin + xmax) / 2.0;
    let range = (xmax - xmin) / 2.0;

    // Compute light curve
    let ncol = if mdl.t2.value > 0.0 { 5 } else { 4 };
    let mut fcomp = vec![vec![0.0; ncol]; n];

    if mdl.iscale {
        fcomp = (0..n).into_par_iter().map(|np| {
            // Compute phase with quadratic correction
            let mut phase = (data[np].time - mdl.t0.value) / mdl.period.value;
            for _ in 0..4 {
                phase -= (mdl.t0.value + phase * (mdl.period.value + mdl.pdot.value * phase) - data[np].time)
                    / (mdl.period.value + 2.0 * mdl.pdot.value * phase);
            }
            // Roemer delay
            phase += mdl.deltat.value / mdl.period.value / 2.0 * ((2.0 * PI * phase).cos() - 1.0);

            let expose = data[np].expose / mdl.period.value;
            let frac = (data[np].time - middle) / range;
            let slfac = 1.0 + frac * (mdl.slope.value + frac * (mdl.quad.value + frac * mdl.cube.value));

            let mut row = vec![0.0; ncol];
            row[0] = slfac * crate::flux::comp_star1(
                mdl.iangle.value, &ldc1, phase, expose, data[np].ndiv,
                mdl.q.value, mdl.beam_factor1.value, mdl.velocity_scale.value as f32,
                &gint, &star1f, &star1c);

            row[1] = slfac * crate::flux::comp_disc(
                mdl.iangle.value, mdl.lin_limb_disc.value, mdl.quad_limb_disc.value,
                phase, expose, data[np].ndiv, mdl.q.value,
                mdl.velocity_scale.value as f32, &disc);

            row[2] = slfac * crate::flux::comp_edge(
                mdl.iangle.value, mdl.lin_limb_disc.value, mdl.quad_limb_disc.value,
                phase, expose, data[np].ndiv, mdl.q.value,
                mdl.velocity_scale.value as f32, &edge);

            row[3] = slfac * crate::flux::comp_spot(
                mdl.iangle.value, phase, expose, data[np].ndiv,
                mdl.q.value, mdl.velocity_scale.value as f32, &spot);

            if mdl.t2.value > 0.0 {
                row[4] = slfac * crate::flux::comp_star2(
                    mdl.iangle.value, &ldc2, phase, expose, data[np].ndiv,
                    mdl.q.value, mdl.beam_factor2.value, mdl.velocity_scale.value as f32,
                    mdl.glens1, rlens1, &gint, &star2f, &star2c);
            }
            row
        }).collect();
    } else {
        calc = (0..n).into_par_iter().map(|np| {
            // Compute phase with quadratic correction
            let mut phase = (data[np].time - mdl.t0.value) / mdl.period.value;
            for _ in 0..4 {
                phase -= (mdl.t0.value + phase * (mdl.period.value + mdl.pdot.value * phase) - data[np].time)
                    / (mdl.period.value + 2.0 * mdl.pdot.value * phase);
            }
            // Roemer delay
            phase += mdl.deltat.value / mdl.period.value / 2.0 * ((2.0 * PI * phase).cos() - 1.0);

            let expose = data[np].expose / mdl.period.value;
            let frac = (data[np].time - middle) / range;
            let slfac = 1.0 + frac * (mdl.slope.value + frac * (mdl.quad.value + frac * mdl.cube.value));

            slfac * crate::flux::comp_light(
                mdl.iangle.value, &ldc1, &ldc2,
                mdl.lin_limb_disc.value, mdl.quad_limb_disc.value,
                phase, expose, data[np].ndiv,
                mdl.q.value, mdl.beam_factor1.value, mdl.beam_factor2.value,
                mdl.spin1.value, mdl.spin2.value, mdl.velocity_scale.value as f32,
                mdl.glens1, rlens1, &gint,
                &star1f, &star2f, &star1c, &star2c, &disc, &edge, &spot)
                + mdl.third.value
        }).collect();
    }

    // White dwarf contribution at phase 0.5
    let mut wdwarf = crate::flux::comp_star1(
        mdl.iangle.value, &ldc1, 0.5, 0.0, 1, mdl.q.value,
        mdl.beam_factor1.value, mdl.velocity_scale.value as f32,
        &gint, &star1f, &star1c);

    let mut sfac = vec![1.0; 5];
    let (mut chisq, mut wnok) = (0.0, 0.0);

    if scale {
        if mdl.iscale {
            // SVD fit for individual scale factors
            let ncomp = ncol;
            sfac = svd_scale(data, &fcomp, ncomp);
            wdwarf *= sfac[0];
            if ncol == 4 {
                // Expand to 5 elements with star2 = 0
                let old = sfac.clone();
                sfac = vec![old[0], old[1], old[2], old[3], 0.0];
            }
            for np in 0..n {
                calc[np] = sfac[0] * fcomp[np][0] + sfac[1] * fcomp[np][1]
                    + sfac[2] * fcomp[np][2] + sfac[3] * fcomp[np][3];
                if mdl.t2.value > 0.0 {
                    calc[np] += sfac[4] * fcomp[np][4];
                }
            }
            (chisq, wnok) = compute_chisq(data, &calc);
        } else {
            let ssfac = re_scale(data, &mut calc);
            wdwarf *= ssfac.0;
            chisq = ssfac.1;
            wnok = ssfac.2;
            sfac = vec![ssfac.0; 5];
        }
    } else {
        wdwarf *= sfac[0];
        if mdl.iscale {
            for np in 0..n {
                calc[np] = sfac[0] * fcomp[np][0] + sfac[1] * fcomp[np][1]
                    + sfac[2] * fcomp[np][2] + sfac[3] * fcomp[np][3];
                if mdl.t2.value > 0.0 {
                    calc[np] += sfac[4] * fcomp[np][4];
                }
            }
        } else {
            for c in calc.iter_mut() {
                *c *= sfac[0];
            }
        }
        if rdata {
            let r = compute_chisq(data, &calc);
            chisq = r.0;
            wnok = r.1;
        }
    }

    let logg1 = crate::flux::comp_gravity1(mdl, &star1f);
    let logg2 = crate::flux::comp_gravity2(mdl, &star2f);
    let rv1 = crate::flux::comp_radius1(&star1f);
    let rv2 = crate::flux::comp_radius2(&star2f);

    Ok(LcResult { calc, wdwarf, chisq, wnok, logg1, logg2, rv1, rv2, sfac })
}

/// Evaluate chi-squared for a batch of parameter sets in parallel.
///
/// * `base` - template model (cloned for each parameter set)
/// * `data` - observed data points (read once, shared across all evaluations)
/// * `param_names` - names of parameters being varied
/// * `param_values` - flat row-major array of length `N * ndim`
/// * `scale` - whether to autoscale each model to the data
///
/// Returns a `Vec<f64>` of length N.  Sets that cause errors produce `f64::NAN`.
pub fn chisq_batch(
    base: &Model,
    data: &Data,
    param_names: &[&str],
    param_values: &[f64],
    scale: bool,
) -> Vec<f64> {
    let ndim = param_names.len();
    assert!(
        ndim > 0 && param_values.len() % ndim == 0,
        "param_values length must be a multiple of param_names length"
    );
    let n = param_values.len() / ndim;

    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut mdl = base.clone();
            let row = &param_values[i * ndim..(i + 1) * ndim];
            for (j, name) in param_names.iter().enumerate() {
                mdl.set_param_value(name, row[j]);
            }
            match light_curve_comp(&mdl, data, scale, true) {
                Ok(res) => res.chisq,
                Err(_) => f64::NAN,
            }
        })
        .collect()
}

/// Re-scale fit to minimize chi-squared with a single global factor.
/// Returns (scale, chisq, wnok).
fn re_scale(data: &Data, fit: &mut [f64]) -> (f64, f64, f64) {
    let mut sdy = 0.0;
    let mut syy = 0.0;
    let mut wnok = 0.0;
    for (i, d) in data.iter().enumerate() {
        if d.weight > 0.0 {
            let wgt = d.weight / sqr(d.ferr);
            sdy += wgt * d.flux * fit[i];
            syy += wgt * fit[i] * fit[i];
            wnok += d.weight;
        }
    }
    if wnok > 0.0 && syy > 0.0 {
        let scale = sdy / syy;
        for f in fit.iter_mut() {
            *f *= scale;
        }
        let mut chisq = 0.0;
        for (i, d) in data.iter().enumerate() {
            if d.weight > 0.0 {
                chisq += d.weight * sqr((d.flux - fit[i]) / d.ferr);
            }
        }
        (scale, chisq, wnok)
    } else {
        (1.0, 0.0, 0.0)
    }
}

/// Compute chi-squared from data and fit.
fn compute_chisq(data: &Data, calc: &[f64]) -> (f64, f64) {
    let mut chisq = 0.0;
    let mut wnok = 0.0;
    for (i, d) in data.iter().enumerate() {
        if d.weight > 0.0 {
            wnok += d.weight;
            chisq += d.weight * sqr((d.flux - calc[i]) / d.ferr);
        }
    }
    (chisq, wnok)
}

/// Simple SVD-based scale factor determination.
/// Fits data = sum_j sfac[j] * fcomp[i][j] using least squares.
fn svd_scale(data: &Data, fcomp: &[Vec<f64>], ncol: usize) -> Vec<f64> {
    use nalgebra::{DMatrix, DVector};

    let n = data.len();
    let mut nok = 0;
    for d in data {
        if d.weight > 0.0 { nok += 1; }
    }
    if nok < ncol || nok == 0 {
        return vec![1.0; ncol];
    }

    let mut a = DMatrix::zeros(nok, ncol);
    let mut b = DVector::zeros(nok);
    let mut row = 0;
    for (i, d) in data.iter().enumerate() {
        if d.weight > 0.0 {
            let w = d.weight.sqrt() / d.ferr;
            for j in 0..ncol {
                a[(row, j)] = fcomp[i][j] * w;
            }
            b[row] = d.flux * w;
            row += 1;
        }
    }

    let svd = a.svd(true, true);
    match svd.solve(&b, 1e-12) {
        Ok(x) => x.iter().copied().collect(),
        Err(_) => vec![1.0; ncol],
    }
}
