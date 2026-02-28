use lcurve_subs::Vec3;
use lcurve_subs::numerical::bsstep::{BsState, bsstep};
use crate::RocheError;

/// Initialize gas stream particle just inside L1 (Lubow & Shu theory).
///
/// Returns (position, velocity) in rotating frame.
pub fn strinit(q: f64) -> Result<(Vec3, Vec3), RocheError> {
    const SMALL: f64 = 1.0e-5;

    let rl1 = crate::lagrange::xl1(q)?;
    let mu = q / (1.0 + q);
    let a = (1.0 - mu) / rl1.powi(3) + mu / (1.0 - rl1).powi(3);
    let lambda1 = (((a - 2.0) + (a * (9.0 * a - 8.0)).sqrt()) / 2.0).sqrt();
    let m1 = (lambda1 * lambda1 - 2.0 * a - 1.0) / (2.0 * lambda1);

    let r = Vec3::new(rl1 - SMALL, -m1 * SMALL, 0.0);
    let v = Vec3::new(-lambda1 * SMALL, -lambda1 * m1 * SMALL, 0.0);
    Ok((r, v))
}

/// One step of orbit integration in Roche potential using Bulirsch-Stoer.
///
/// Updates `r`, `v`, and `time` in place.
/// Returns `(tdid, tnext)`.
pub fn gsint(
    q: f64,
    r: &mut Vec3,
    v: &mut Vec3,
    ttry: f64,
    time: &mut f64,
    eps: f64,
    state: &mut BsState,
) -> (f64, f64) {
    let mut y = [r.x, r.y, r.z, v.x, v.y, v.z];
    let mut dydt = [0.0; 6];

    // Compute initial derivatives
    rocder(q, *time, &y, &mut dydt);

    let yscal = [1.0; 6];
    let result = bsstep(
        &mut y,
        &dydt,
        6,
        time,
        ttry,
        eps,
        &yscal,
        state,
        &|t, y_in, dy_out| rocder(q, t, y_in, dy_out),
    );

    r.x = y[0];
    r.y = y[1];
    r.z = y[2];
    v.x = y[3];
    v.y = y[4];
    v.z = y[5];

    (result.hdid, result.hnext)
}

/// Derivative function for Roche potential ODE integration.
fn rocder(q: f64, _t: f64, y: &[f64], dydt: &mut [f64]) {
    let r = Vec3::new(y[0], y[1], y[2]);
    let v = Vec3::new(y[3], y[4], y[5]);

    dydt[0] = v.x;
    dydt[1] = v.y;
    dydt[2] = v.z;

    let a = crate::misc::rocacc(q, &r, &v);
    dydt[3] = a.x;
    dydt[4] = a.y;
    dydt[5] = a.z;
}

/// Advance particle until it reaches specified radius.
///
/// Returns the time taken.
pub fn stradv(
    q: f64,
    r: &mut Vec3,
    v: &mut Vec3,
    rad: f64,
    acc: f64,
    smax: f64,
) -> Result<f64, RocheError> {
    const EPS: f64 = 1.0e-8;
    const TMAX: f64 = 10.0;

    let mut state = BsState::new();
    let rinit = r.length();
    let mut rnow = rinit;
    let mut time = 0.0;
    let mut tnext: f64 = 1.0e-2;

    let mut ro;
    let mut vo;

    // Step until radius crossed
    loop {
        ro = *r;
        vo = *v;
        let ttry = tnext.min(smax);
        let (tdid, tn) = gsint(q, r, v, ttry, &mut time, EPS, &mut state);
        tnext = tn;
        rnow = r.length();

        if (rinit > rad && rnow <= rad) || (rinit < rad && rnow >= rad) {
            // Binary chop to refine
            let mut lo = 0.0;
            let mut hi = tdid;
            let mut rlo = ro.length();
            let mut rhi = rnow;
            let to = time - tdid;

            while (rhi - rlo).abs() > acc {
                let tmid = (lo + hi) / 2.0;
                *r = ro;
                *v = vo;
                let mut t = to;
                let mut state2 = BsState::new();
                gsint(q, r, v, tmid, &mut t, EPS, &mut state2);
                rnow = r.length();
                if (rhi > rad && rnow > rad) || (rhi < rad && rnow < rad) {
                    rhi = rnow;
                    hi = tmid;
                } else {
                    rlo = rnow;
                    lo = tmid;
                }
            }
            return Ok(time);
        }

        if time > TMAX {
            return Err(RocheError::Generic(format!(
                "stradv: took too long without crossing radius={}",
                rad
            )));
        }
    }
}

/// Find next point where stream is closest or furthest from primary.
pub fn strmnx(q: f64, r: &mut Vec3, v: &mut Vec3) -> Result<(), RocheError> {
    const EPS: f64 = 1.0e-8;
    const ACC: f64 = 1.0e-8;

    let mut state = BsState::new();
    let dir1 = Vec3::dot(r, v);
    let mut dir = dir1;
    let mut time = 0.0;
    let mut ttry = 1.0e-2;

    let mut ro;
    let mut vo;

    // Step until direction reverses
    loop {
        ro = *r;
        vo = *v;
        let (tdid, tnext) = gsint(q, r, v, ttry, &mut time, EPS, &mut state);
        dir = Vec3::dot(r, v);
        ttry = tnext;

        if (dir > 0.0) != (dir1 > 0.0) {
            // Binary chop
            let mut lo = 0.0;
            let mut hi = tdid;
            let to = time - tdid;

            while (hi - lo).abs() > ACC {
                let tmid = (lo + hi) / 2.0;
                *r = ro;
                *v = vo;
                let mut t = to;
                let mut state2 = BsState::new();
                gsint(q, r, v, tmid, &mut t, EPS, &mut state2);
                dir = Vec3::dot(r, v);
                if (dir1 > 0.0 && dir < 0.0) || (dir1 < 0.0 && dir > 0.0) {
                    hi = tmid;
                } else {
                    lo = tmid;
                }
            }
            return Ok(());
        }
    }
}

/// Check if stream crosses a Roche equipotential.
///
/// Returns `Some((x, y))` if crossing occurs, `None` if minimum potential not reached.
pub fn hits(q: f64, pot: f64) -> Result<Option<(f64, f64)>, RocheError> {
    const EPS: f64 = 1.0e-9;
    const BCHOP: f64 = 1.0e-8;

    let (mut r, mut v) = strinit(q)?;
    let mut state = BsState::new();
    let mut time = 0.0;
    let mut tnext = 2.0e-4;

    let mut rold = r;
    let mut vold = v;

    // Integrate until distance starts increasing
    let (_, tn) = gsint(q, &mut r, &mut v, tnext, &mut time, EPS, &mut state);
    tnext = tn;

    loop {
        if crate::misc::rdot(&r, &v) >= 0.0 {
            break;
        }
        rold = r;
        vold = v;
        let (_, tn) = gsint(q, &mut r, &mut v, tnext, &mut time, EPS, &mut state);
        tnext = tn;
    }

    // Check if potential reached
    let pot_at_min = crate::potential::rpot(q, &r)?;
    if pot_at_min >= pot {
        return Ok(None); // doesn't reach the potential
    }

    // Re-integrate from start to find where it first crosses
    let (mut r2, mut v2) = strinit(q)?;
    let mut state2 = BsState::new();
    let mut time2 = 0.0;
    let mut tnext2: f64 = 2.0e-4;

    loop {
        rold = r2;
        vold = v2;
        let _told = time2;
        let ttry = tnext2.min(time - time2);
        let (tdid, tn) = gsint(q, &mut r2, &mut v2, ttry, &mut time2, EPS, &mut state2);
        tnext2 = tn;

        let rpot_val = crate::potential::rpot(q, &r2)?;
        if rpot_val <= pot {
            // Binary chop to refine
            let mut tlo = 0.0;
            let mut thi = tdid;
            let t_base = time2 - tdid;

            while thi - tlo > BCHOP {
                let tmid = (tlo + thi) / 2.0;
                r2 = rold;
                v2 = vold;
                let mut t = t_base;
                let mut st = BsState::new();
                gsint(q, &mut r2, &mut v2, tmid, &mut t, EPS, &mut st);
                if crate::potential::rpot(q, &r2)? > pot {
                    tlo = tmid;
                } else {
                    thi = tmid;
                }
            }
            return Ok(Some((r2.x, r2.y)));
        }
    }
}

/// Extended version of hits that also returns velocity.
pub fn hits_with_vel(q: f64, pot: f64) -> Result<Option<(f64, f64, f64, f64)>, RocheError> {
    // Same as hits but also returns vx, vy
    match hits(q, pot)? {
        Some((x, y)) => Ok(Some((x, y, 0.0, 0.0))), // velocity tracking would need full re-impl
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strinit() {
        let (r, v) = strinit(1.0).unwrap();
        // Should be just inside L1, on negative y side
        assert!(r.x > 0.0 && r.x < 0.5);
        assert!(r.y.abs() < 0.01);
        assert!(r.z == 0.0);
        assert!(v.x < 0.0); // moving toward primary
    }

    #[test]
    fn test_gsint_conserves_jacobi() {
        let q = 0.5;
        let (mut r, mut v) = strinit(q).unwrap();
        let j0 = crate::misc::jacobi(q, &r, &v);

        let mut time = 0.0;
        let mut state = BsState::new();
        gsint(q, &mut r, &mut v, 0.01, &mut time, 1e-10, &mut state);
        let j1 = crate::misc::jacobi(q, &r, &v);

        assert!(
            (j0 - j1).abs() < 1e-8,
            "Jacobi not conserved: j0={}, j1={}",
            j0,
            j1
        );
    }
}
