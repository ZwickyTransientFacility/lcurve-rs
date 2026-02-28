use lcurve_subs::Vec3;
use crate::{RocheError, Star};

/// Compute reference sphere radius and potential for a Roche-distorted star.
///
/// The reference sphere just encloses the star surface along the line of centres.
///
/// * `q` - mass ratio M2/M1
/// * `star` - which star
/// * `spin` - spin/orbital frequency ratio
/// * `ffac` - linear filling factor (0-1, where 1 = Roche-filling)
///
/// Returns `(rref, pref)` — reference radius and reference potential.
pub fn ref_sphere(q: f64, star: Star, spin: f64, ffac: f64) -> Result<(f64, f64), RocheError> {
    let rl1 = match star {
        Star::Primary => crate::lagrange::xl11(q, spin)?,
        Star::Secondary => 1.0 - crate::lagrange::xl12(q, spin)?,
    };

    // Linear filling factor: surface at ffac * rl1 distance
    let tref = rl1;

    // Compute potential at the actual surface radius
    let surf_r = ffac * tref;
    let p = match star {
        Star::Primary => Vec3::new(surf_r, 0.0, 0.0),
        Star::Secondary => Vec3::new(1.0 - surf_r, 0.0, 0.0),
    };

    let pref = match star {
        Star::Primary => crate::potential::rpot1(q, spin, &p)?,
        Star::Secondary => crate::potential::rpot2(q, spin, &p)?,
    };

    // Expand reference sphere by 0.1% but cap at Roche lobe radius
    let rref = (1.001 * ffac).min(1.0) * tref;

    Ok((rref, pref))
}

/// Find surface point on Roche-distorted star in given direction.
///
/// Uses binary chop to locate point on equipotential surface.
///
/// * `q` - mass ratio
/// * `star` - which star
/// * `spin` - spin/orbital frequency
/// * `dirn` - direction to search (from star center)
/// * `rref` - reference sphere radius
/// * `pref` - reference potential (defines the surface)
/// * `acc` - fractional accuracy
///
/// Returns `(position, normal, radius, gravity)`.
pub fn face(
    q: f64,
    star: Star,
    spin: f64,
    dirn: &Vec3,
    rref: f64,
    pref: f64,
    acc: f64,
) -> Result<(Vec3, Vec3, f64, f64), RocheError> {
    let cofm = match star {
        Star::Primary => Vec3::ZERO,
        Star::Secondary => Vec3::new(1.0, 0.0, 0.0),
    };

    let pot_fn = |radius: f64| -> Result<f64, RocheError> {
        let p = cofm + radius * *dirn;
        match star {
            Star::Primary => crate::potential::rpot1(q, spin, &p),
            Star::Secondary => crate::potential::rpot2(q, spin, &p),
        }
    };

    // Find r1, r2 such that r1 is below reference potential and r2 is above.
    // Search downward from rref in halving steps (matches C++ ref_sphere).
    let mut lo = rref / 2.0;
    let mut hi = rref;
    let mut tref = pref + 1.0;

    for _ in 0..30 {
        if tref <= pref {
            break;
        }
        lo = hi / 2.0;
        tref = pot_fn(lo)?;
        if tref > pref {
            hi = lo;
        }
    }

    // Binary chop with absolute accuracy (matches C++)
    for _ in 0..100 {
        if hi - lo <= acc {
            break;
        }
        let mid = (lo + hi) / 2.0;
        let pot_mid = pot_fn(mid)?;
        if pot_mid < pref {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let radius = (lo + hi) / 2.0;
    let pvec = cofm + radius * *dirn;

    // Compute gradient (normal direction and gravity)
    let grad = match star {
        Star::Primary => crate::potential::drpot1(q, spin, &pvec)?,
        Star::Secondary => crate::potential::drpot2(q, spin, &pvec)?,
    };

    let gravity = grad.length();
    let mut dvec = grad;
    let _ = dvec.unit(); // normalize to get outward normal

    Ok((pvec, dvec, radius, gravity))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ref_sphere() {
        let (rref, pref) = ref_sphere(1.0, Star::Secondary, 1.0, 1.0).unwrap();
        assert!(rref > 0.0);
        assert!(pref < 0.0);
    }

    #[test]
    fn test_face_primary() {
        let q = 0.5;
        let spin = 1.0;
        let ffac = 0.9;
        let (rref, pref) = ref_sphere(q, Star::Primary, spin, ffac).unwrap();
        let dirn = Vec3::new(0.0, 1.0, 0.0); // perpendicular to line of centres

        let (pvec, dvec, r, g) = face(q, Star::Primary, spin, &dirn, rref, pref, 1e-6).unwrap();
        assert!(r > 0.0);
        assert!(g > 0.0);
        assert!(pvec.y > 0.0);
    }
}
