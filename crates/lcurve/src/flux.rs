use lcurve_subs::{sqr, Vec3, TWOPI};
use lcurve_roche::misc::set_earth_cosi_sini as set_earth;
use crate::types::{Ginterp, LDC, Point};

/// Compute total flux at a given phase with exposure smearing.
#[allow(clippy::too_many_arguments)]
pub fn comp_light(
    iangle: f64, ldc1: &LDC, ldc2: &LDC,
    lin_limb_disc: f64, quad_limb_disc: f64,
    phase: f64, expose: f64, ndiv: i32, q: f64,
    beam_factor1: f64, beam_factor2: f64,
    spin1: f64, spin2: f64, vscale: f32,
    glens1: bool, rlens1: f64,
    gint: &Ginterp,
    star1f: &[Point], star2f: &[Point],
    star1c: &[Point], star2c: &[Point],
    disc: &[Point], edge: &[Point], spot: &[Point],
) -> f64 {
    let xcofm = q / (1.0 + q);
    let cosi = (TWOPI * iangle / 360.0).cos();
    let sini = (TWOPI * iangle / 360.0).sin();
    let vfac = vscale as f64 / (lcurve_subs::C / 1e3);

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);

        let ptype = gint.grid_type(phi);
        let star1 = if ptype == 1 { star1f } else { star1c };
        let star2 = if ptype == 3 { star2f } else { star2c };

        let mut ssum = 0.0;

        // Star 1
        for pt in star1 {
            if pt.visible(phi) {
                let mu = Vec3::dot(&earth, &pt.dirn);
                if ldc1.see(mu) {
                    if beam_factor1 != 0.0 {
                        let vx = -vfac * spin1 * pt.posn.y;
                        let vy = vfac * (spin1 * pt.posn.x - xcofm);
                        let vr = -(earth.x * vx + earth.y * vy);
                        let vn = pt.dirn.x * vx + pt.dirn.y * vy;
                        let mud = mu - mu * vr - vn;
                        ssum += mu * pt.flux as f64 * (1.0 - beam_factor1 * vr) * ldc1.imu(mud);
                    } else {
                        ssum += mu * pt.flux as f64 * ldc1.imu(mu);
                    }
                }
            }
        }
        ssum *= gint.scale1(phi);

        // Star 2
        let mut ssum2 = 0.0;
        for pt in star2 {
            if pt.visible(phi) {
                let mu = Vec3::dot(&earth, &pt.dirn);
                if ldc2.see(mu) {
                    let mag = if glens1 {
                        lensing_mag(&pt.posn, &earth, rlens1)
                    } else {
                        1.0
                    };

                    if beam_factor2 != 0.0 {
                        let vx = -vfac * spin2 * pt.posn.y;
                        let vy = vfac * (spin2 * (pt.posn.x - 1.0) + 1.0 - xcofm);
                        let vr = -(earth.x * vx + earth.y * vy);
                        let vn = pt.dirn.x * vx + pt.dirn.y * vy;
                        let mud = mu - mu * vr - vn;
                        ssum2 += mu * mag * pt.flux as f64 * (1.0 - beam_factor2 * vr) * ldc2.imu(mud);
                    } else {
                        ssum2 += mu * mag * pt.flux as f64 * ldc2.imu(mu);
                    }
                }
            }
        }
        ssum += gint.scale2(phi) * ssum2;

        // Disc
        for pt in disc {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                let ommu = 1.0 - mu;
                ssum += mu * pt.flux as f64 * (1.0 - ommu * (lin_limb_disc + quad_limb_disc * ommu));
            }
        }

        // Edge
        for pt in edge {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                let ommu = 1.0 - mu;
                ssum += mu * pt.flux as f64 * (1.0 - ommu * (lin_limb_disc + quad_limb_disc * ommu));
            }
        }

        // Spot
        for pt in spot {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                ssum += mu * pt.flux as f64;
            }
        }

        sum += wgt * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux from star 1 only.
#[allow(clippy::too_many_arguments)]
pub fn comp_star1(
    iangle: f64, ldc1: &LDC, phase: f64, expose: f64, ndiv: i32,
    q: f64, beam_factor1: f64, vscale: f32,
    gint: &Ginterp, star1f: &[Point], star1c: &[Point],
) -> f64 {
    let xcofm = q / (1.0 + q);
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();
    let vfac = vscale as f64 / (lcurve_subs::C / 1e3);

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);
        let star1 = if gint.grid_type(phi) == 1 { star1f } else { star1c };

        let mut ssum = 0.0;
        for pt in star1 {
            if pt.visible(phi) {
                let mu = Vec3::dot(&earth, &pt.dirn);
                if ldc1.see(mu) {
                    if beam_factor1 != 0.0 {
                        let vx = -vfac * pt.posn.y;
                        let vy = vfac * (pt.posn.x - xcofm);
                        let vr = -(earth.x * vx + earth.y * vy);
                        let vn = pt.dirn.x * vx + pt.dirn.y * vy;
                        let mud = mu - mu * vr - vn;
                        ssum += mu * pt.flux as f64 * (1.0 - beam_factor1 * vr) * ldc1.imu(mud);
                    } else {
                        let vr = vfac * (earth.x * pt.posn.y - earth.y * (pt.posn.x - xcofm));
                        ssum += mu * pt.flux as f64 * (1.0 - beam_factor1 * vr) * ldc1.imu(mu);
                    }
                }
            }
        }
        sum += wgt * gint.scale1(phase) * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux from star 2 only.
#[allow(clippy::too_many_arguments)]
pub fn comp_star2(
    iangle: f64, ldc2: &LDC, phase: f64, expose: f64, ndiv: i32,
    q: f64, beam_factor2: f64, vscale: f32,
    glens1: bool, rlens1: f64,
    gint: &Ginterp, star2f: &[Point], star2c: &[Point],
) -> f64 {
    let xcofm = q / (1.0 + q);
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();
    let vfac = vscale as f64 / (lcurve_subs::C / 1e3);

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);
        let star2 = if gint.grid_type(phi) == 3 { star2f } else { star2c };

        let mut ssum = 0.0;
        for pt in star2 {
            if pt.visible(phi) {
                let mu = Vec3::dot(&earth, &pt.dirn);
                if ldc2.see(mu) {
                    let mag = if glens1 { lensing_mag(&pt.posn, &earth, rlens1) } else { 1.0 };

                    if beam_factor2 != 0.0 {
                        let vx = -vfac * pt.posn.y;
                        let vy = vfac * (pt.posn.x - xcofm);
                        let vr = -(earth.x * vx + earth.y * vy);
                        let vn = pt.dirn.x * vx + pt.dirn.y * vy;
                        let mud = mu - mu * vr - vn;
                        ssum += mu * mag * pt.flux as f64 * (1.0 - beam_factor2 * vr) * ldc2.imu(mud);
                    } else {
                        let vr = vfac * (earth.x * pt.posn.y - earth.y * (pt.posn.x - xcofm));
                        ssum += mu * mag * pt.flux as f64 * (1.0 - beam_factor2 * vr) * ldc2.imu(mu);
                    }
                }
            }
        }
        sum += wgt * gint.scale2(phi) * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux from disc.
pub fn comp_disc(
    iangle: f64, lin_limb_disc: f64, quad_limb_disc: f64,
    phase: f64, expose: f64, ndiv: i32, _q: f64, _vscale: f32,
    disc: &[Point],
) -> f64 {
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);
        let mut ssum = 0.0;
        for pt in disc {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                let ommu = 1.0 - mu;
                ssum += mu * pt.flux as f64 * (1.0 - ommu * (lin_limb_disc + quad_limb_disc * ommu));
            }
        }
        sum += wgt * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux from disc edge.
pub fn comp_edge(
    iangle: f64, lin_limb_disc: f64, quad_limb_disc: f64,
    phase: f64, expose: f64, ndiv: i32, _q: f64, _vscale: f32,
    edge: &[Point],
) -> f64 {
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);
        let mut ssum = 0.0;
        for pt in edge {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                let ommu = 1.0 - mu;
                ssum += mu * pt.flux as f64 * (1.0 - ommu * (lin_limb_disc + quad_limb_disc * ommu));
            }
        }
        sum += wgt * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux from bright spot.
pub fn comp_spot(
    iangle: f64, phase: f64, expose: f64, ndiv: i32, _q: f64, _vscale: f32,
    spot: &[Point],
) -> f64 {
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();

    let mut sum = 0.0;
    for nd in 0..ndiv {
        let (phi, wgt) = exposure_weight(phase, expose, nd, ndiv);
        let earth = set_earth(cosi, sini, phi);
        let mut ssum = 0.0;
        for pt in spot {
            let mu = Vec3::dot(&earth, &pt.dirn);
            if mu > 0.0 && pt.visible(phi) {
                ssum += mu * pt.flux as f64;
            }
        }
        sum += wgt * ssum;
    }
    sum / 1.max(ndiv - 1) as f64
}

/// Compute flux-weighted gravity for star 1 (returns log10(g) in CGS).
pub fn comp_gravity1(mdl: &crate::model::Model, star1: &[Point]) -> f64 {
    let gm = (1000.0 * mdl.velocity_scale.value).powi(3) * mdl.tperiod * lcurve_subs::DAY / TWOPI;
    let a = (gm / sqr(TWOPI / lcurve_subs::DAY / mdl.tperiod)).powf(1.0 / 3.0);
    let gscale = 100.0 * gm / sqr(a);

    let (r1, _) = mdl.get_r1r2();
    let rl1 = lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin1.value).unwrap_or(0.5);

    let gref = if mdl.roche1 {
        let acc = mdl.delta_phase / 10.0;
        let ffac1 = r1 / rl1;
        let (rref1, pref1) = lcurve_roche::surface::ref_sphere(
            mdl.q.value, lcurve_roche::Star::Primary, mdl.spin1.value, ffac1).unwrap();
        let dirn = Vec3::new(-1.0, 0.0, 0.0);
        let (_, _, _, g) = lcurve_roche::surface::face(
            mdl.q.value, lcurve_roche::Star::Primary, mdl.spin1.value, &dirn, rref1, pref1, acc).unwrap();
        g * gscale
    } else {
        gscale / (1.0 + mdl.q.value) / sqr(r1)
    };

    let (mut sumfg, mut sumf) = (0.0, 0.0);
    for pt in star1 {
        sumfg += pt.flux as f64 * pt.gravity as f64;
        sumf += pt.flux as f64;
    }
    if gref > 0.0 && sumfg > 0.0 && sumf > 0.0 {
        (gref * sumfg / sumf).log10()
    } else {
        0.0
    }
}

/// Compute flux-weighted gravity for star 2 (returns log10(g) in CGS).
pub fn comp_gravity2(mdl: &crate::model::Model, star2: &[Point]) -> f64 {
    let gm = (1000.0 * mdl.velocity_scale.value).powi(3) * mdl.tperiod * lcurve_subs::DAY / TWOPI;
    let a = (gm / sqr(TWOPI / lcurve_subs::DAY / mdl.tperiod)).powf(1.0 / 3.0);
    let gscale = 100.0 * gm / sqr(a);

    let (_, r2) = mdl.get_r1r2();
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value).unwrap_or(0.5);

    let gref = if mdl.roche2 {
        let acc = mdl.delta_phase / 10.0;
        let ffac2 = r2 / rl2;
        let (rref2, pref2) = lcurve_roche::surface::ref_sphere(
            mdl.q.value, lcurve_roche::Star::Secondary, mdl.spin2.value, ffac2).unwrap();
        let dirn = Vec3::new(1.0, 0.0, 0.0);
        let (_, _, _, g) = lcurve_roche::surface::face(
            mdl.q.value, lcurve_roche::Star::Secondary, mdl.spin2.value, &dirn, rref2, pref2, acc).unwrap();
        g * gscale
    } else {
        gscale * mdl.q.value / (1.0 + mdl.q.value) / sqr(r2)
    };

    let (mut sumfg, mut sumf) = (0.0, 0.0);
    for pt in star2 {
        sumfg += pt.flux as f64 * pt.gravity as f64;
        sumf += pt.flux as f64;
    }
    if gref > 0.0 && sumfg > 0.0 && sumf > 0.0 {
        (gref * sumfg / sumf).log10()
    } else {
        0.0
    }
}

/// Compute volume-averaged radius for star 1.
pub fn comp_radius1(star1: &[Point]) -> f64 {
    let cofm1 = Vec3::ZERO;
    let (mut sumsa, mut sumvol) = (0.0, 0.0);
    for pt in star1 {
        let vec = pt.posn - cofm1;
        let r = vec.length();
        let rcosa = Vec3::dot(&pt.dirn, &vec);
        sumsa += pt.area as f64 * rcosa / r.powi(3);
        sumvol += pt.area as f64 * rcosa;
    }
    (sumvol / sumsa).powf(1.0 / 3.0)
}

/// Compute volume-averaged radius for star 2.
pub fn comp_radius2(star2: &[Point]) -> f64 {
    let cofm2 = Vec3::new(1.0, 0.0, 0.0);
    let (mut sumsa, mut sumvol) = (0.0, 0.0);
    for pt in star2 {
        let vec = pt.posn - cofm2;
        let r = vec.length();
        let rcosa = Vec3::dot(&pt.dirn, &vec);
        sumsa += pt.area as f64 * rcosa / r.powi(3);
        sumvol += pt.area as f64 * rcosa;
    }
    (sumvol / sumsa).powf(1.0 / 3.0)
}

/// Compute exposure sub-division phase and weight.
#[inline]
fn exposure_weight(phase: f64, expose: f64, nd: i32, ndiv: i32) -> (f64, f64) {
    if ndiv == 1 {
        (phase, 1.0)
    } else {
        let phi = phase + expose * (nd as f64 - (ndiv - 1) as f64 / 2.0) / (ndiv - 1) as f64;
        let wgt = if nd == 0 || nd == ndiv - 1 { 0.5 } else { 1.0 };
        (phi, wgt)
    }
}

/// Compute gravitational lensing magnification.
fn lensing_mag(posn: &Vec3, earth: &Vec3, rlens1: f64) -> f64 {
    let s = *posn;
    let d = -Vec3::dot(&s, earth);
    if d <= 0.0 {
        return 1.0;
    }
    let p = (s + d * *earth).length();
    let ph = p / 2.0;
    let phsq = ph * ph;
    let rd = rlens1 * d;
    let pd = if phsq > 25.0 * rd {
        p + rd / p
    } else {
        ph + (phsq + rd).sqrt()
    };
    pd * pd / (pd - ph) / ph / 4.0
}
