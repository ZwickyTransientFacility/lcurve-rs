use lcurve_subs::{sqr, Vec3, TWOPI};

/// Acceleration in rotating Roche potential (includes Coriolis).
pub fn rocacc(q: f64, r: &Vec3, v: &Vec3) -> Vec3 {
    let f1 = 1.0 / (1.0 + q);
    let f2 = f1 * q;
    let yzsq = sqr(r.y) + sqr(r.z);
    let r1sq = sqr(r.x) + yzsq;
    let r2sq = sqr(r.x - 1.0) + yzsq;
    let fm1 = f1 / (r1sq * r1sq.sqrt());
    let fm2 = f2 / (r2sq * r2sq.sqrt());
    let fm3 = fm1 + fm2;
    Vec3::new(
        -fm3 * r.x + fm2 + 2.0 * v.y + r.x - f2,
        -fm3 * r.y - 2.0 * v.x + r.y,
        -fm3 * r.z,
    )
}

/// Jacobi constant (energy integral in the rotating frame).
pub fn jacobi(q: f64, r: &Vec3, v: &Vec3) -> f64 {
    let f1 = 1.0 / (1.0 + q);
    let f2 = f1 * q;
    let yzsq = sqr(r.y) + sqr(r.z);
    (sqr(v.x) + sqr(v.y) + sqr(v.z) - sqr(r.y) - sqr(r.x - f2)) / 2.0
        - f1 / (sqr(r.x) + yzsq).sqrt()
        - f2 / (sqr(r.x - 1.0) + yzsq).sqrt()
}

/// dr/dt = r·v / |r|.
pub fn rdot(r: &Vec3, v: &Vec3) -> f64 {
    Vec3::dot(r, v) / r.length()
}

/// Unit vector toward Earth given inclination (degrees) and orbital phase.
pub fn set_earth(iangle: f64, phase: f64) -> Vec3 {
    let irad = iangle.to_radians();
    set_earth_cosi_sini(irad.cos(), irad.sin(), phase)
}

/// Unit vector toward Earth given pre-computed cos(i), sin(i) and phase.
pub fn set_earth_cosi_sini(cosi: f64, sini: f64, phase: f64) -> Vec3 {
    let phi = TWOPI * phase;
    Vec3::new(sini * phi.cos(), -sini * phi.sin(), cosi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_earth_overhead() {
        let e = set_earth(0.0, 0.0);
        assert!((e.z - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_set_earth_edgeon() {
        // At phase=0, i=90: observer in +x direction
        let e = set_earth(90.0, 0.0);
        assert!((e.x - 1.0).abs() < 1e-10);
        assert!(e.y.abs() < 1e-10);
        assert!(e.z.abs() < 1e-10);

        // At phase=0.25, i=90: observer in -y direction
        let e = set_earth(90.0, 0.25);
        assert!(e.x.abs() < 1e-10);
        assert!((e.y - (-1.0)).abs() < 1e-10);
    }
}
