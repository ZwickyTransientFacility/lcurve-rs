use lcurve_subs::{sqr, Vec3};
use crate::RocheError;

/// Synchronous Roche potential.
///
/// * `q` - mass ratio M2/M1 (must be > 0)
/// * `p` - position (units of binary separation)
pub fn rpot(q: f64, p: &Vec3) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("rpot: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let comp = 1.0 - mu;
    let x2y2 = sqr(p.x) + sqr(p.y);
    let r1sq = x2y2 + sqr(p.z);
    let r1 = r1sq.sqrt();
    let r2 = (r1sq + 1.0 - 2.0 * p.x).sqrt();
    Ok(-comp / r1 - mu / r2 - (x2y2 + mu * (mu - 2.0 * p.x)) / 2.0)
}

/// Asynchronous Roche potential for primary star.
pub fn rpot1(q: f64, spin: f64, p: &Vec3) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("rpot1: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let comp = 1.0 - mu;
    let x2y2 = sqr(p.x) + sqr(p.y);
    let r1sq = x2y2 + sqr(p.z);
    let r1 = r1sq.sqrt();
    let r2 = (r1sq + 1.0 - 2.0 * p.x).sqrt();
    Ok(-comp / r1 - mu / r2 - sqr(spin) * x2y2 / 2.0 + mu * p.x)
}

/// Asynchronous Roche potential for secondary star.
pub fn rpot2(q: f64, spin: f64, p: &Vec3) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("rpot2: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let comp = 1.0 - mu;
    let x2y2 = sqr(p.x) + sqr(p.y);
    let r1sq = x2y2 + sqr(p.z);
    let r1 = r1sq.sqrt();
    let r2 = (r1sq + 1.0 - 2.0 * p.x).sqrt();
    Ok(-comp / r1 - mu / r2 - sqr(spin) * (0.5 + 0.5 * x2y2 - p.x) - comp * p.x)
}

/// Gradient of synchronous Roche potential.
pub fn drpot(q: f64, p: &Vec3) -> Result<Vec3, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("drpot: q={} <= 0", q)));
    }
    let r1sq = p.sqr();
    let r2sq = r1sq + 1.0 - 2.0 * p.x;
    let r1 = r1sq.sqrt();
    let r2 = r2sq.sqrt();
    let mu = q / (1.0 + q);
    let mu1 = mu / r2 / r2sq;
    let comp = (1.0 - mu) / r1 / r1sq;
    Ok(Vec3::new(
        comp * p.x + mu1 * (p.x - 1.0) - p.x + mu,
        comp * p.y + mu1 * p.y - p.y,
        comp * p.z + mu1 * p.z,
    ))
}

/// Gradient of asynchronous Roche potential for primary star.
pub fn drpot1(q: f64, spin: f64, p: &Vec3) -> Result<Vec3, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("drpot1: q={} <= 0", q)));
    }
    let r1sq = p.sqr();
    let r2sq = r1sq + 1.0 - 2.0 * p.x;
    let r1 = r1sq.sqrt();
    let r2 = r2sq.sqrt();
    let mu = q / (1.0 + q);
    let mu1 = mu / r2 / r2sq;
    let comp = (1.0 - mu) / r1 / r1sq;
    let ssq = sqr(spin);
    Ok(Vec3::new(
        comp * p.x + mu1 * (p.x - 1.0) - ssq * p.x + mu,
        comp * p.y + mu1 * p.y - ssq * p.y,
        comp * p.z + mu1 * p.z,
    ))
}

/// Gradient of asynchronous Roche potential for secondary star.
pub fn drpot2(q: f64, spin: f64, p: &Vec3) -> Result<Vec3, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("drpot2: q={} <= 0", q)));
    }
    let r1sq = p.sqr();
    let r2sq = r1sq + 1.0 - 2.0 * p.x;
    let r1 = r1sq.sqrt();
    let r2 = r2sq.sqrt();
    let mu = q / (1.0 + q);
    let mu1 = mu / r2 / r2sq;
    let comp = (1.0 - mu) / r1 / r1sq;
    let ssq = sqr(spin);
    Ok(Vec3::new(
        comp * p.x + mu1 * (p.x - 1.0) - ssq * (p.x - 1.0) + mu - 1.0,
        comp * p.y + mu1 * p.y - ssq * p.y,
        comp * p.z + mu1 * p.z,
    ))
}

/// Compute potential value along the line-of-sight: point + lam * earth.
/// Used by fblink and pot_min for 1D minimization.
pub fn rpot_along_los(q: f64, star: crate::Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> f64 {
    let point = *p + lam * *earth;
    match star {
        crate::Star::Primary => rpot1(q, spin, &point).unwrap_or(0.0),
        crate::Star::Secondary => rpot2(q, spin, &point).unwrap_or(0.0),
    }
}

/// Compute potential gradient along the line-of-sight direction.
pub fn drpot_along_los(q: f64, star: crate::Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> f64 {
    let point = *p + lam * *earth;
    let grad = match star {
        crate::Star::Primary => drpot1(q, spin, &point).unwrap_or(Vec3::ZERO),
        crate::Star::Secondary => drpot2(q, spin, &point).unwrap_or(Vec3::ZERO),
    };
    Vec3::dot(&grad, earth)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rpot_symmetry() {
        let p1 = rpot(1.0, &Vec3::new(0.5, 0.0, 0.0)).unwrap();
        assert!(p1.is_finite());
    }

    #[test]
    fn test_rpot_negative() {
        let p = rpot(1.0, &Vec3::new(0.1, 0.0, 0.0)).unwrap();
        assert!(p < 0.0);
    }

    #[test]
    fn test_drpot_at_l1() {
        let q = 1.0;
        let xl1 = crate::lagrange::xl1(q).unwrap();
        let grad = drpot(q, &Vec3::new(xl1, 0.0, 0.0)).unwrap();
        assert!(grad.x.abs() < 1e-10, "grad.x={}", grad.x);
        assert!(grad.y.abs() < 1e-15);
        assert!(grad.z.abs() < 1e-15);
    }
}
