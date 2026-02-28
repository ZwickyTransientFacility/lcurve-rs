use lcurve_subs::{sqr, Vec3, PI, TWOPI};
use lcurve_roche::{self, Star};
use crate::model::Model;
use crate::types::{EclipseVec, Point};
use crate::LcurveError;

/// Compute number of surface elements given grid parameters.
pub fn numface(nlat: i32, infill: bool, thelo: f64, thehi: f64, nlatfill: i32, nlngfill: i32) -> usize {
    let mut nface = 0usize;
    let dtheta = PI / nlat as f64;

    if infill {
        let nl1 = (thelo / dtheta).ceil() as i32;
        for i in 0..nl1 {
            let theta = thelo * (i as f64 + 0.5) / nl1 as f64;
            let dphi = dtheta / theta.sin();
            nface += 16.max((TWOPI / dphi) as usize);
        }
        let nl2 = (((1 + nlatfill) as f64 * (thehi - thelo)) / dtheta).ceil() as i32;
        for i in 0..nl2 {
            let theta = thelo + (thehi - thelo) * (i as f64 + 0.5) / nl2 as f64;
            let dphi = dtheta / theta.sin() / (1 + nlngfill) as f64;
            nface += 8.max((PI / dphi) as usize);
        }
        let nl3 = ((thehi - thelo) / dtheta).ceil() as i32;
        for i in 0..nl3 {
            let theta = thelo + (thehi - thelo) * (i as f64 + 0.5) / nl3 as f64;
            let dphi = dtheta / theta.sin();
            nface += 8.max((PI / dphi) as usize);
        }
        let nl4 = ((PI - thehi) / dtheta).ceil() as i32;
        for i in 0..nl4 {
            let theta = thehi + (PI - thehi) * (i as f64 + 0.5) / nl4 as f64;
            let dphi = dtheta / theta.sin();
            nface += 16.max((TWOPI / dphi) as usize);
        }
    } else {
        let nlat = (PI / dtheta).ceil() as i32;
        for i in 0..nlat {
            let theta = PI * (i as f64 + 0.5) / nlat as f64;
            let dphi = dtheta / theta.sin();
            nface += 16.max((TWOPI / dphi) as usize);
        }
    }
    nface
}

/// Convenience wrapper for star eclipse computation.
pub fn star_eclipse(
    q: f64, r: f64, spin: f64, ffac: f64, iangle: f64, posn: &Vec3,
    delta: f64, roche: bool, star: Star, eclipses: &mut EclipseVec,
) -> Result<(), LcurveError> {
    let ri = iangle.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();
    let cofm = match star {
        Star::Primary => Vec3::ZERO,
        Star::Secondary => Vec3::new(1.0, 0.0, 0.0),
    };

    if roche {
        if let Some((ingress, egress)) = lcurve_roche::eclipse::ingress_egress(
            q, star, spin, ffac, iangle, delta, posn,
        )? {
            eclipses.push((ingress, egress));
        }
    } else if let Some((phi1, phi2, _lam1, _lam2)) =
        lcurve_roche::eclipse::sphere_eclipse_phase(cosi, sini, posn, &cofm, r)
    {
        eclipses.push((phi1, phi2));
    }
    Ok(())
}

/// Set up grid elements for a star (primary or secondary).
pub fn set_star_grid(
    mdl: &Model, which_star: Star, fine: bool,
) -> Result<Vec<Point>, LcurveError> {
    let (mut r1, mut r2) = mdl.get_r1r2();

    let eclipse = match which_star {
        Star::Primary => mdl.eclipse1,
        Star::Secondary => mdl.eclipse2,
    };
    let nlat = match which_star {
        Star::Primary => if fine { mdl.nlat1f } else { mdl.nlat1c },
        Star::Secondary => if fine { mdl.nlat2f } else { mdl.nlat2c },
    };
    let nlatfill = match which_star {
        Star::Primary => 0,
        Star::Secondary => if fine { mdl.nlatfill } else { 0 },
    };
    let nlngfill = match which_star {
        Star::Primary => 0,
        Star::Secondary => if fine { mdl.nlngfill } else { 0 },
    };

    let rl1 = lcurve_roche::lagrange::xl11(mdl.q.value, mdl.spin1.value)?;
    if r1 < 0.0 { r1 = rl1; } else if r1 > rl1 {
        return Err(LcurveError::Generic("set_star_grid: primary star larger than Roche lobe".into()));
    }
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value)?;
    if r2 < 0.0 { r2 = rl2; } else if r2 > rl2 {
        return Err(LcurveError::Generic("set_star_grid: secondary star larger than Roche lobe".into()));
    }

    let ffac1 = r1 / rl1;
    let (rref1, pref1) = lcurve_roche::surface::ref_sphere(mdl.q.value, Star::Primary, mdl.spin1.value, ffac1)?;
    let ffac2 = r2 / rl2;
    let (rref2, pref2) = lcurve_roche::surface::ref_sphere(mdl.q.value, Star::Secondary, mdl.spin2.value, ffac2)?;

    // Determine infill latitude range (only for secondary fine grid)
    let infill = mdl.npole && matches!(which_star, Star::Secondary)
        && (nlatfill > 0 || nlngfill > 0) && r2 > r1;

    let (mut thelo, mut thehi) = (0.0, 0.0);
    if infill {
        let rangle = mdl.iangle.value.to_radians();
        let cosi = rangle.cos();
        let ratio = (cosi + r1) / r2;
        if ratio >= 1.0 {
            thehi = PI / 2.0 + rangle;
        } else {
            let mut llo_val = 0.0f64;
            let mut lhi_val = PI / 2.0;
            while lhi_val > llo_val + 1e-7 {
                let lmid = (llo_val + lhi_val) / 2.0;
                let xy = envelope(rangle, lmid, r1);
                if xy.0 * xy.0 + xy.1 * xy.1 < r2 * r2 {
                    llo_val = lmid;
                } else {
                    lhi_val = lmid;
                }
            }
            let sini = rangle.sin();
            let xy = envelope(rangle, (llo_val + lhi_val) / 2.0, r1);
            thehi = (PI / 2.0 - mdl.llo.to_radians())
                .max((PI / 2.0 + rangle)
                    .min((xy.1 * sini / r2).acos() + mdl.lfudge.to_radians()));
        }
        let ratio2 = (cosi - r1) / r2;
        if ratio2 >= 1.0 {
            return set_star_grid_inner(mdl, which_star, r1, r2, rref1, pref1, rref2, pref2,
                                       ffac1, ffac2, eclipse, nlat, false, 0.0, 0.0, 0, 0);
        } else if ratio2 <= -1.0 {
            thelo = 0.0;
        } else {
            thelo = (PI / 2.0 - mdl.lhi.to_radians())
                .min((PI / 2.0 - ratio2.acos() + rangle - mdl.lfudge.to_radians()).max(0.0));
        }
    }

    set_star_grid_inner(mdl, which_star, r1, r2, rref1, pref1, rref2, pref2,
                         ffac1, ffac2, eclipse, nlat, infill, thelo, thehi, nlatfill, nlngfill)
}

fn set_star_grid_inner(
    mdl: &Model, which_star: Star, r1: f64, r2: f64,
    rref1: f64, pref1: f64, rref2: f64, pref2: f64,
    ffac1: f64, ffac2: f64, eclipse: bool, nlat: i32,
    infill: bool, thelo: f64, thehi: f64, nlatfill: i32, nlngfill: i32,
) -> Result<Vec<Point>, LcurveError> {
    let acc = mdl.delta_phase / 10.0;
    let dtheta = PI / nlat as f64;

    // Compute reference gravity
    let gref = match which_star {
        Star::Primary if mdl.roche1 => {
            let dirn = Vec3::new(-1.0, 0.0, 0.0);
            let (_, _, _, g) = lcurve_roche::surface::face(
                mdl.q.value, Star::Primary, mdl.spin1.value, &dirn, rref1, pref1, acc)?;
            g
        }
        Star::Secondary if mdl.roche2 => {
            let dirn = Vec3::new(1.0, 0.0, 0.0);
            let (_, _, _, g) = lcurve_roche::surface::face(
                mdl.q.value, Star::Secondary, mdl.spin2.value, &dirn, rref2, pref2, acc)?;
            g
        }
        _ => 1.0,
    };

    let mut star = Vec::new();

    if infill {
        add_faces(&mut star, 0.0, thelo, dtheta, 0, 0, mdl, which_star,
                  r1, r2, rref1, rref2, pref1, pref2, ffac1, ffac2, eclipse, gref, acc)?;
        add_faces(&mut star, thelo, thehi, dtheta, nlatfill, nlngfill, mdl, which_star,
                  r1, r2, rref1, rref2, pref1, pref2, ffac1, ffac2, eclipse, gref, acc)?;
        add_faces(&mut star, thehi, PI, dtheta, 0, 0, mdl, which_star,
                  r1, r2, rref1, rref2, pref1, pref2, ffac1, ffac2, eclipse, gref, acc)?;
    } else {
        add_faces(&mut star, 0.0, PI, dtheta, 0, 0, mdl, which_star,
                  r1, r2, rref1, rref2, pref1, pref2, ffac1, ffac2, eclipse, gref, acc)?;
    }

    Ok(star)
}

#[allow(clippy::too_many_arguments)]
fn add_faces(
    star: &mut Vec<Point>,
    tlo: f64, thi: f64, dtheta: f64,
    nlatfill: i32, nlngfill: i32,
    mdl: &Model, which_star: Star,
    r1: f64, r2: f64, rref1: f64, rref2: f64, pref1: f64, pref2: f64,
    ffac1: f64, ffac2: f64, eclipse: bool, gref: f64, acc: f64,
) -> Result<(), LcurveError> {
    let cofm1 = Vec3::ZERO;
    let cofm2 = Vec3::new(1.0, 0.0, 0.0);
    let ri = mdl.iangle.value.to_radians();
    let cosi = ri.cos();
    let sini = ri.sin();

    let infill = nlatfill > 0 || nlngfill > 0;
    let (nlat, nlat1, nlat2);
    if infill {
        nlat1 = ((1 + nlatfill) as f64 * (thi - tlo) / dtheta).ceil() as i32;
        nlat2 = ((thi - tlo) / dtheta).ceil() as i32;
        nlat = nlat1 + nlat2;
    } else {
        nlat1 = 0;
        nlat2 = 0;
        nlat = ((thi - tlo) / dtheta).ceil() as i32;
    }

    for nt in 0..nlat {
        let (theta, nphi, nl, phi1, phi2);
        if infill {
            if nt < nlat1 {
                theta = tlo + (thi - tlo) * (nt as f64 + 0.5) / nlat1 as f64;
                let sint = theta.sin();
                nphi = 8.max((PI * sint * (1 + nlngfill) as f64 / dtheta) as i32);
                nl = nlat1;
                match which_star {
                    Star::Primary => { phi1 = -PI / 2.0; phi2 = PI / 2.0; }
                    Star::Secondary => { phi1 = PI / 2.0; phi2 = 3.0 * PI / 2.0; }
                }
            } else {
                theta = tlo + (thi - tlo) * ((nt - nlat1) as f64 + 0.5) / nlat2 as f64;
                let sint = theta.sin();
                nphi = 8.max((PI * sint / dtheta) as i32);
                nl = nlat2;
                match which_star {
                    Star::Primary => { phi1 = PI / 2.0; phi2 = 3.0 * PI / 2.0; }
                    Star::Secondary => { phi1 = -PI / 2.0; phi2 = PI / 2.0; }
                }
            }
        } else {
            theta = tlo + (thi - tlo) * (nt as f64 + 0.5) / nlat as f64;
            let sint = theta.sin();
            nphi = 16.max((TWOPI * sint / dtheta) as i32);
            nl = nlat;
            phi1 = 0.0;
            phi2 = TWOPI;
        }

        let sint = theta.sin();
        let cost = theta.cos();

        for np in 0..nphi {
            let phi = phi1 + (phi2 - phi1) * (np as f64 + 0.5) / nphi as f64;
            let sinp = phi.sin();
            let cosp = phi.cos();

            let dirn = if mdl.npole {
                Vec3::new(sint * cosp, sint * sinp, cost)
            } else {
                Vec3::new(cost, sint * cosp, sint * sinp)
            };

            let (posn, dvec, rad, gravity);
            match which_star {
                Star::Primary if mdl.roche1 => {
                    let r = lcurve_roche::surface::face(
                        mdl.q.value, Star::Primary, mdl.spin1.value, &dirn, rref1, pref1, acc)?;
                    posn = r.0; dvec = r.1; rad = r.2; gravity = r.3;
                }
                Star::Secondary if mdl.roche2 => {
                    let r = lcurve_roche::surface::face(
                        mdl.q.value, Star::Secondary, mdl.spin2.value, &dirn, rref2, pref2, acc)?;
                    posn = r.0; dvec = r.1; rad = r.2; gravity = r.3;
                }
                Star::Primary => {
                    rad = r1;
                    posn = cofm1 + rad * dirn;
                    dvec = dirn;
                    gravity = 1.0;
                }
                Star::Secondary => {
                    rad = r2;
                    posn = cofm2 + rad * dirn;
                    dvec = dirn;
                    gravity = 1.0;
                }
            }

            let area = ((phi2 - phi1) / nphi as f64 * rad * sint)
                * ((thi - tlo) / nl as f64 * rad)
                / Vec3::dot(&dirn, &dvec);

            let mut eclipses = Vec::new();
            if eclipse {
                match which_star {
                    Star::Primary => {
                        star_eclipse(mdl.q.value, r2, mdl.spin2.value, ffac2,
                                     mdl.iangle.value, &posn, mdl.delta_phase,
                                     mdl.roche2, Star::Secondary, &mut eclipses)?;
                    }
                    Star::Secondary => {
                        star_eclipse(mdl.q.value, r1, mdl.spin1.value, ffac1,
                                     mdl.iangle.value, &posn, mdl.delta_phase,
                                     mdl.roche1, Star::Primary, &mut eclipses)?;
                    }
                }
            }

            star.push(Point::new(posn, dvec, area, gravity / gref, eclipses));
        }
    }
    Ok(())
}

fn envelope(rangle: f64, lambda: f64, r1: f64) -> (f64, f64) {
    let sini = rangle.sin();
    let cosi = rangle.cos();
    let sinl = lambda.sin();
    let cosl = lambda.cos();
    let norm = (cosi * cosi + sini * sini * cosl * cosl).sqrt();
    let x = sinl + r1 * cosi * sinl / norm;
    let y = -cosi * cosl - r1 * cosl / norm;
    (x, y)
}

/// Set up disc grid elements (top surface).
pub fn set_disc_grid(mdl: &Model) -> Result<Vec<Point>, LcurveError> {
    const EFAC: f64 = 1.0000001;

    let (mut r1, mut r2) = mdl.get_r1r2();
    let rl1 = lcurve_roche::lagrange::xl11(mdl.q.value, mdl.spin1.value)?;
    if r1 < 0.0 { r1 = rl1; }
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value)?;
    if r2 < 0.0 { r2 = rl2; }

    let rdisc1 = if mdl.rdisc1.value > 0.0 { mdl.rdisc1.value } else { r1 };
    let rdisc2 = if mdl.rdisc2.value > 0.0 { mdl.rdisc2.value } else { mdl.radius_spot.value };

    let ffac1 = r1 / rl1;
    let ffac2 = r2 / rl2;

    let drad = (rdisc2 - rdisc1) / mdl.nrad as f64;
    let drrad = rdisc2 / mdl.nrad as f64;

    let mut disc = Vec::new();

    for i in 0..mdl.nrad {
        let rad = rdisc1 + (rdisc2 - rdisc1) * (i as f64 + 0.5) / mdl.nrad as f64;
        let ntheta = 8.max((TWOPI * rad / drrad).ceil() as i32);
        let h = EFAC * mdl.height_disc.value * rad.powf(mdl.beta_disc.value);
        let tanp = mdl.beta_disc.value * h / rad;
        let cosp = 1.0 / (1.0 + sqr(tanp)).sqrt();
        let sinp = tanp * cosp;
        let area = TWOPI / ntheta as f64 * rad * drad;

        for j in 0..ntheta {
            let theta = TWOPI * j as f64 / ntheta as f64;
            let sint = theta.sin();
            let cost = theta.cos();
            let posn = Vec3::new(rad * cost, rad * sint, h);
            let dirn = Vec3::new(-cost * sinp, -sint * sinp, cosp);

            let mut eclipses = Vec::new();
            if mdl.opaque {
                eclipses = lcurve_roche::disc_eclipse::disc_eclipse(
                    mdl.iangle.value, rdisc1, rdisc2,
                    mdl.beta_disc.value, mdl.height_disc.value, &posn);
            }
            if mdl.eclipse1 {
                star_eclipse(mdl.q.value, r1, mdl.spin1.value, ffac1,
                             mdl.iangle.value, &posn, mdl.delta_phase,
                             mdl.roche1, Star::Primary, &mut eclipses)?;
            }
            if mdl.eclipse2 {
                star_eclipse(mdl.q.value, r2, mdl.spin2.value, ffac2,
                             mdl.iangle.value, &posn, mdl.delta_phase,
                             mdl.roche2, Star::Secondary, &mut eclipses)?;
            }
            disc.push(Point::new(posn, dirn, area, 1.0, eclipses));
        }
    }
    Ok(disc)
}

/// Set up disc edge elements (outer or inner rim).
pub fn set_disc_edge(mdl: &Model, outer: bool) -> Result<Vec<Point>, LcurveError> {
    const EFAC: f64 = 1.0000001;

    let (mut r1, mut r2) = mdl.get_r1r2();
    let rl1 = lcurve_roche::lagrange::xl11(mdl.q.value, mdl.spin1.value)?;
    if r1 < 0.0 { r1 = rl1; }
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value)?;
    if r2 < 0.0 { r2 = rl2; }

    let rdisc1 = if mdl.rdisc1.value > 0.0 { mdl.rdisc1.value } else { r1 };
    let rdisc2 = if mdl.rdisc2.value > 0.0 { mdl.rdisc2.value } else { mdl.radius_spot.value };
    let size = rdisc2 / mdl.nrad as f64;

    let ffac1 = r1 / rl1;
    let ffac2 = r2 / rl2;

    let rad_val = if outer { rdisc2 } else { rdisc1 };
    let h = mdl.height_disc.value * rad_val.powf(mdl.beta_disc.value);
    let nout = (2.0 * h / size) as i32 + 2;
    let ntheta = 8.max((TWOPI * rad_val / size).ceil() as i32);

    let area = TWOPI / ntheta as f64 * rad_val * 2.0 * h / (nout - 1) as f64;
    let mut edge = Vec::new();

    for i in 0..ntheta {
        let theta = TWOPI * i as f64 / ntheta as f64;
        let sint = theta.sin();
        let cost = theta.cos();

        // Upper rim element
        let posn = if outer {
            Vec3::new(EFAC * rad_val * cost, EFAC * rad_val * sint, EFAC * h)
        } else {
            Vec3::new(rad_val * cost / EFAC, rad_val * sint / EFAC, EFAC * h)
        };
        let dirn = Vec3::new(cost, sint, 0.0);

        let mut eclipses = Vec::new();
        if mdl.opaque {
            eclipses = lcurve_roche::disc_eclipse::disc_eclipse(
                mdl.iangle.value, rdisc1, rdisc2,
                mdl.beta_disc.value, mdl.height_disc.value, &posn);
        }
        if mdl.eclipse1 {
            let _ = star_eclipse(mdl.q.value, r1, mdl.spin1.value, ffac1,
                                 mdl.iangle.value, &posn, mdl.delta_phase,
                                 mdl.roche1, Star::Primary, &mut eclipses);
        }
        if mdl.eclipse2 {
            star_eclipse(mdl.q.value, r2, mdl.spin2.value, ffac2,
                         mdl.iangle.value, &posn, mdl.delta_phase,
                         mdl.roche2, Star::Secondary, &mut eclipses)?;
        }
        edge.push(Point::new(posn, dirn, area / 2.0, 1.0, eclipses));

        // Lower elements
        for j in 0..(nout - 1) {
            let z = -h + 2.0 * h * j as f64 / (nout - 1) as f64;
            let (posn, dirn) = if outer {
                (Vec3::new(EFAC * rad_val * cost, EFAC * rad_val * sint, z),
                 Vec3::new(cost, sint, 0.0))
            } else {
                (Vec3::new(rad_val * cost / EFAC, rad_val * sint / EFAC, z),
                 Vec3::new(-cost, -sint, 0.0))
            };

            let mut eclipses = Vec::new();
            if mdl.opaque {
                eclipses = lcurve_roche::disc_eclipse::disc_eclipse(
                    mdl.iangle.value, rdisc1, rdisc2,
                    mdl.beta_disc.value, mdl.height_disc.value, &posn);
            }
            if mdl.eclipse1 {
                let _ = star_eclipse(mdl.q.value, r1, mdl.spin1.value, ffac1,
                                     mdl.iangle.value, &posn, mdl.delta_phase,
                                     mdl.roche1, Star::Primary, &mut eclipses);
            }
            if mdl.eclipse2 {
                star_eclipse(mdl.q.value, r2, mdl.spin2.value, ffac2,
                             mdl.iangle.value, &posn, mdl.delta_phase,
                             mdl.roche2, Star::Secondary, &mut eclipses)?;
            }
            edge.push(Point::new(posn, dirn, area, 1.0, eclipses));
        }
    }
    Ok(edge)
}

/// Set up bright spot grid elements.
pub fn set_bright_spot_grid(mdl: &Model) -> Result<Vec<Point>, LcurveError> {
    let (mut r1, mut r2) = mdl.get_r1r2();
    let rl1 = lcurve_roche::lagrange::xl11(mdl.q.value, mdl.spin1.value)?;
    if r1 < 0.0 { r1 = rl1; } else if r1 > rl1 {
        return Err(LcurveError::Generic("set_bright_spot_grid: primary larger than Roche lobe".into()));
    }
    let rl2 = 1.0 - lcurve_roche::lagrange::xl12(mdl.q.value, mdl.spin2.value)?;
    if r2 < 0.0 { r2 = rl2; } else if r2 > rl2 {
        return Err(LcurveError::Generic("set_bright_spot_grid: secondary larger than Roche lobe".into()));
    }

    let ffac1 = r1 / rl1;
    let ffac2 = r2 / rl2;

    // Locate spot position via gas stream
    let (mut bspot, mut v) = lcurve_roche::stream::strinit(mdl.q.value)?;
    lcurve_roche::stream::stradv(mdl.q.value, &mut bspot, &mut v,
                                  mdl.radius_spot.value, 1e-10, 1e-3)?;

    let theta = bspot.y.atan2(bspot.x) + TWOPI / 4.0 + mdl.angle_spot.value.to_radians();
    let alpha = mdl.yaw_spot.value.to_radians();
    let tilt = mdl.tilt_spot.value.to_radians();

    let bvec = Vec3::new(theta.cos(), theta.sin(), 0.0);
    let pvec = Vec3::new(0.0, 0.0, 1.0);
    let tvec = Vec3::new(
        tilt.sin() * (theta + alpha).sin(),
        -tilt.sin() * (theta + alpha).cos(),
        tilt.cos(),
    );

    let bmax = (mdl.expon_spot.value / mdl.epow_spot.value).powf(1.0 / mdl.epow_spot.value);
    let sfac = 20.0 + bmax;
    let area = sfac * mdl.length_spot.value * mdl.height_spot.value / (mdl.nspot - 1) as f64;
    let bright_ref = lcurve_subs::planck(mdl.wavelength, mdl.temp_spot.value);

    let mut spot = Vec::with_capacity(2 * mdl.nspot as usize);

    for i in 0..mdl.nspot {
        let dist = sfac * i as f64 / (mdl.nspot - 1) as f64;
        let posn = bspot + mdl.length_spot.value * (dist - bmax) * bvec;

        let mut eclipses = Vec::new();
        if mdl.eclipse1 {
            star_eclipse(mdl.q.value, r1, mdl.spin1.value, ffac1,
                         mdl.iangle.value, &posn, mdl.delta_phase,
                         mdl.roche1, Star::Primary, &mut eclipses)?;
        }
        if mdl.eclipse2 {
            star_eclipse(mdl.q.value, r2, mdl.spin2.value, ffac2,
                         mdl.iangle.value, &posn, mdl.delta_phase,
                         mdl.roche2, Star::Secondary, &mut eclipses)?;
        }

        let bright = bright_ref * (dist / bmax).powf(mdl.expon_spot.value)
            * (mdl.expon_spot.value / mdl.epow_spot.value - dist.powf(mdl.epow_spot.value)).exp();

        // Tilted strip
        let mut pt_tilt = Point::new(posn, tvec, area, 1.0, eclipses.clone());
        pt_tilt.flux = (bright * (1.0 - mdl.cfrac_spot.value) * area) as f32;
        spot.push(pt_tilt);

        // Parallel (upward) strip
        let mut pt_par = Point::new(posn, pvec, area, 1.0, eclipses);
        pt_par.flux = (bright * mdl.cfrac_spot.value * area) as f32;
        spot.push(pt_par);
    }

    Ok(spot)
}
