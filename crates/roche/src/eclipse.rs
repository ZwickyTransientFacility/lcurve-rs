use lcurve_subs::{sqr, Vec3, TWOPI};
use crate::{RocheError, Star};

/// Test if point `r` is occulted by Roche-lobe-filling secondary.
///
/// Uses sphere pre-screening then steps along photon path.
/// Returns true if eclipsed.
pub fn blink(q: f64, r: &Vec3, e: &Vec3, acc: f64) -> Result<bool, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("blink: q={} <= 0", q)));
    }
    if acc <= 0.0 {
        return Err(RocheError::Generic(format!("blink: acc={} <= 0", acc)));
    }

    let xcm = q / (1.0 + q);
    let c1 = 2.0 / (1.0 + q);
    let c2 = 2.0 * xcm;

    let rl1 = crate::lagrange::xl1(q)?;

    // Critical potential at L1
    let r1_l1 = rl1;
    let r2_l1 = 1.0 - rl1;
    let xc_l1 = rl1 - xcm;
    let crit = c1 / r1_l1 + c2 / r2_l1 + xc_l1 * xc_l1;

    let rsphere = 1.0 - rl1;
    let pp = rsphere * rsphere;
    let step = rsphere * acc;

    // Sphere test: does photon path cross the sphere centered on secondary?
    let xt = r.x - 1.0;
    let b = e.x * xt + e.y * r.y + e.z * r.z;
    let c_val = xt * xt + r.y * r.y + r.z * r.z - pp;

    let fac = b * b - c_val;
    if fac <= 0.0 {
        return Ok(false);
    }
    let fac = fac.sqrt();

    let par2 = -b + fac;
    if par2 <= 0.0 {
        return Ok(false);
    }

    let par1 = (-b - fac).max(0.0);

    // Start at closest approach
    let par = if b < 0.0 { -b } else { 0.0 };
    let x1 = r.x + par * e.x;
    let x2 = r.y + par * e.y;
    let x3 = r.z + par * e.z;

    let xm = x1 - 1.0;
    let rr = x2 * x2 + x3 * x3;
    let rs2 = xm * xm + rr;

    if rs2 <= 0.0 {
        return Ok(true);
    }

    let r1 = (x1 * x1 + rr).sqrt();
    let r2 = rs2.sqrt();
    let xc = x1 - xcm;
    let yy = x2 * x2;
    let mut c_pot = c1 / r1 + c2 / r2 + xc * xc + yy;

    if c_pot > crit {
        return Ok(true);
    }

    // Determine step direction from derivative
    let a_val = x1 * e.x + x2 * e.y;
    let b_val = a_val + x3 * e.z;
    let rs1 = x1 * x1 + rr;
    let deriv = -c1 * b_val / (rs1 * r1) - c2 * (b_val - e.x) / (rs2 * r2)
        + 2.0 * (a_val - xcm * e.x);

    let (p1, p2) = if deriv > 0.0 {
        (par, par2)
    } else {
        (par, par1)
    };

    let nstep = ((p2 - p1).abs() / step + 0.5).floor() as usize;
    let nstep = nstep.max(2);
    let dp = (p2 - p1) / nstep as f64;
    let mut cmax = c_pot - 1.0;

    for i in 1..=nstep {
        if c_pot <= cmax {
            break;
        }
        cmax = c_pot;
        let p = p1 + dp * i as f64;
        let x1 = r.x + p * e.x;
        let x2 = r.y + p * e.y;
        let x3 = r.z + p * e.z;

        let xm = x1 - 1.0;
        let rr = x2 * x2 + x3 * x3;
        let r2 = (xm * xm + rr).sqrt();
        if r2 <= 0.0 {
            return Ok(true);
        }
        let r1 = (x1 * x1 + rr).sqrt();
        let xc = x1 - xcm;
        let yy = x2 * x2;
        c_pot = c1 / r1 + c2 / r2 + xc * xc + yy;
        if c_pot > crit {
            return Ok(true);
        }
    }

    Ok(false)
}

/// Test if point `p` is eclipsed by Roche-distorted star (general, asynchronous).
///
/// Searches along line-of-sight for minimum potential below surface value.
pub fn fblink(
    q: f64,
    star: Star,
    spin: f64,
    ffac: f64,
    acc: f64,
    earth: &Vec3,
    p: &Vec3,
) -> Result<bool, RocheError> {
    let (rref, pref) = crate::surface::ref_sphere(q, star, spin, ffac)?;

    let cofm = match star {
        Star::Primary => Vec3::ZERO,
        Star::Secondary => Vec3::new(1.0, 0.0, 0.0),
    };

    let (lam1, lam2) = match sphere_eclipse_los(earth, p, &cofm, rref) {
        Some((l1, l2)) => (l1, l2),
        None => return Ok(false),
    };

    if lam1 == 0.0 {
        return Ok(true);
    }

    // 1D search along lambda for minimum potential
    let func = |lam: f64| crate::potential::rpot_along_los(q, star, spin, earth, p, lam);

    let mut nstep = 1usize;
    let mut step_size = lam2 - lam1;
    let mut f1 = 0.0f64;
    let mut f2 = 0.0f64;
    let mut flam = 1.0;
    let mut lam = lam1;

    while step_size > acc {
        lam = lam1 + step_size / 2.0;
        let mut bracketted = false;
        for n in 0..nstep {
            flam = func(lam + n as f64 * step_size);
            if flam <= pref {
                return Ok(true);
            }
            if nstep == 1 {
                f1 = func(lam1);
                f2 = func(lam2);
            }
            if flam < f1 && flam < f2 {
                lam = lam + n as f64 * step_size;
                bracketted = true;
                break;
            }
        }
        if bracketted {
            break;
        }
        step_size /= 2.0;
        nstep *= 2;
    }

    if flam < f1 && flam < f2 {
        // Bracketed minimum — refine with dbrent
        let mut func_mut = |x: f64| crate::potential::rpot_along_los(q, star, spin, earth, p, x);
        let mut dfunc_mut = |x: f64| crate::potential::drpot_along_los(q, star, spin, earth, p, x);
        let result = lcurve_subs::numerical::dbrent::dbrent(
            lam1, lam, lam2, &mut func_mut, &mut dfunc_mut, acc, true, pref,
        );
        match result {
            Ok((fmin, _xmin)) => Ok(fmin < pref),
            Err(_) => Ok(false),
        }
    } else {
        Ok(false)
    }
}

/// Test if sphere eclipses a point, returning multiplier range along line-of-sight.
///
/// Returns `Some((lam1, lam2))` if sphere eclipses the point, `None` otherwise.
pub fn sphere_eclipse_los(
    earth: &Vec3,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
) -> Option<(f64, f64)> {
    let d = *r - *c;
    let bquad = Vec3::dot(earth, &d);
    if bquad >= 0.0 {
        return None;
    }

    let cquad = d.sqr() - sqr(rsphere);
    let fac = sqr(bquad) - cquad;
    if fac <= 0.0 {
        return None;
    }

    let fac = fac.sqrt();
    let lam2 = -bquad + fac;
    let lam1 = (cquad / lam2).max(0.0);
    Some((lam1, lam2))
}

/// Full sphere eclipse test with phase range computation.
///
/// Returns `Some((phi1, phi2, lam1, lam2))` or `None`.
pub fn sphere_eclipse_phase(
    cosi: f64,
    sini: f64,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
) -> Option<(f64, f64, f64, f64)> {
    let d = *r - *c;
    let pdist = (sqr(d.x) + sqr(d.y)).sqrt();
    let bquad = d.z * cosi - pdist * sini;
    if bquad >= 0.0 {
        return None;
    }

    let cquad = d.sqr() - sqr(rsphere);
    let fac = sqr(bquad) - cquad;
    if fac <= 0.0 {
        return None;
    }

    let fac = fac.sqrt();
    let lam2 = -bquad + fac;
    let lam1 = (cquad / lam2).max(0.0);

    let (phi1, phi2) = if cquad < 0.0 {
        (0.0, 1.0)
    } else {
        let delta = ((cosi * d.z + cquad.sqrt()) / (sini * pdist)).acos();
        let phi = d.y.atan2(-d.x);
        let mut p1 = (phi - delta) / TWOPI;
        p1 -= p1.floor();
        let p2 = p1 + 2.0 * delta / TWOPI;
        (p1, p2)
    };

    Some((phi1, phi2, lam1, lam2))
}

/// Compute ingress and egress phases for eclipse of point `r` by Roche-distorted star.
///
/// Returns `Some((ingress, egress))` if eclipsed, `None` otherwise.
pub fn ingress_egress(
    q: f64,
    star: Star,
    spin: f64,
    ffac: f64,
    iangle: f64,
    delta: f64,
    r: &Vec3,
) -> Result<Option<(f64, f64)>, RocheError> {
    let (rref, pref) = crate::surface::ref_sphere(q, star, spin, ffac)?;

    let irad = iangle.to_radians();
    let cosi = irad.cos();
    let sini = irad.sin();

    let cofm = match star {
        Star::Primary => Vec3::ZERO,
        Star::Secondary => Vec3::new(1.0, 0.0, 0.0),
    };

    // Test for sphere eclipse
    let (phi1, phi2, lam1, lam2) =
        match sphere_eclipse_phase(cosi, sini, r, &cofm, rref) {
            Some(v) => v,
            None => return Ok(None),
        };

    let acc = 2.0 * (2.0 * TWOPI * (lam2 - lam1) * delta).sqrt();

    // Check if potential minimum lies below surface
    if let Some(phi) = pot_min_search(q, cosi, sini, star, spin, ffac, r, phi1, phi2, rref, pref, acc)? {
        // Binary chop for ingress
        let mut pin = phi;
        let mut pout = phi1;
        while (pin - pout).abs() > delta {
            let pmid = (pin + pout) / 2.0;
            let earth = crate::misc::set_earth_cosi_sini(cosi, sini, pmid);
            if fblink(q, star, spin, ffac, acc, &earth, r)? {
                pin = pmid;
            } else {
                pout = pmid;
            }
        }
        let mut ingress = (pin + pout) / 2.0;
        ingress -= ingress.floor();

        // Binary chop for egress
        pin = phi;
        pout = phi2;
        while (pin - pout).abs() > delta {
            let pmid = (pin + pout) / 2.0;
            let earth = crate::misc::set_earth_cosi_sini(cosi, sini, pmid);
            if fblink(q, star, spin, ffac, acc, &earth, r)? {
                pin = pmid;
            } else {
                pout = pmid;
            }
        }
        let mut egress = (pin + pout) / 2.0;
        egress -= egress.floor();
        if egress < ingress {
            egress += 1.0;
        }

        Ok(Some((ingress, egress)))
    } else {
        Ok(None)
    }
}

/// Search for eclipse along line-of-sight cone using fblink.
///
/// Uses a grid search over phases and fblink (with proper 1D line-of-sight
/// minimisation) at each phase to detect eclipses reliably, including
/// grazing eclipses where the potential minimum is not at the midpoint.
///
/// Returns Some(phase) where eclipse occurs, or None.
fn pot_min_search(
    q: f64,
    cosi: f64,
    sini: f64,
    star: Star,
    spin: f64,
    ffac: f64,
    r: &Vec3,
    phi1: f64,
    phi2: f64,
    _rref: f64,
    _pref: f64,
    acc: f64,
) -> Result<Option<f64>, RocheError> {
    let nsteps = ((phi2 - phi1) / (acc / TWOPI)).ceil() as usize;
    let nsteps = nsteps.max(10).min(1000);
    let dphi = (phi2 - phi1) / nsteps as f64;

    for i in 0..=nsteps {
        let phi = phi1 + i as f64 * dphi;
        let earth = crate::misc::set_earth_cosi_sini(cosi, sini, phi);
        if fblink(q, star, spin, ffac, acc, &earth, r)? {
            return Ok(Some(phi));
        }
    }
    Ok(None)
}
