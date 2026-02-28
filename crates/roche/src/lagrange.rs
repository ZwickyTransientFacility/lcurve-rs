use crate::RocheError;

const NMAX: usize = 1000;
const EPS: f64 = 1.0e-12;

/// L1 Lagrangian point: x-distance from primary, normalized by separation.
/// L1 is between the two stars (0 < xl1 < 1).
pub fn xl1(q: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("xl1: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let a1 = -1.0 + mu;
    let a2 = 2.0 - 2.0 * mu;
    let a3 = -1.0 + mu;
    let a4 = 1.0 + 2.0 * mu;
    let a5 = -2.0 - mu;
    let a6 = 1.0f64;
    let (d1, d2, d3, d4, d5) = (a2, 2.0 * a3, 3.0 * a4, 4.0 * a5, 5.0 * a6);

    let mut x = 1.0 / (1.0 + q);
    for _ in 0..NMAX {
        let xold = x;
        let f = x * (x * (x * (x * (x * a6 + a5) + a4) + a3) + a2) + a1;
        let df = x * (x * (x * (x * d5 + d4) + d3) + d2) + d1;
        x -= f / df;
        if (x - xold).abs() <= EPS * x.abs() {
            return Ok(x);
        }
    }
    Err(RocheError::Generic("xl1: exceeded maximum iterations".into()))
}

/// L2 Lagrangian point: beyond secondary (xl2 > 1).
pub fn xl2(q: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("xl2: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let a1 = -1.0 + mu;
    let a2 = 2.0 - 2.0 * mu;
    let a3 = -1.0 - mu;
    let a4 = 1.0 + 2.0 * mu;
    let a5 = -2.0 - mu;
    let a6 = 1.0f64;
    let (d1, d2, d3, d4, d5) = (a2, 2.0 * a3, 3.0 * a4, 4.0 * a5, 5.0 * a6);

    let mut x = 1.5;
    for _ in 0..NMAX {
        let xold = x;
        let f = x * (x * (x * (x * (x * a6 + a5) + a4) + a3) + a2) + a1;
        let df = x * (x * (x * (x * d5 + d4) + d3) + d2) + d1;
        x -= f / df;
        if (x - xold).abs() <= EPS * x.abs() {
            return Ok(x);
        }
    }
    Err(RocheError::Generic("xl2: exceeded maximum iterations".into()))
}

/// L3 Lagrangian point: opposite secondary (xl3 < 0).
pub fn xl3(q: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("xl3: q={} <= 0", q)));
    }
    let mu = q / (1.0 + q);
    let a1 = 1.0 - mu;
    let a2 = -2.0 + 2.0 * mu;
    let a3 = 1.0 - mu;
    let a4 = 1.0 + 2.0 * mu;
    let a5 = -2.0 - mu;
    let a6 = 1.0f64;
    let (d1, d2, d3, d4, d5) = (a2, 2.0 * a3, 3.0 * a4, 4.0 * a5, 5.0 * a6);

    let mut x = -1.0;
    for _ in 0..NMAX {
        let xold = x;
        let f = x * (x * (x * (x * (x * a6 + a5) + a4) + a3) + a2) + a1;
        let df = x * (x * (x * (x * d5 + d4) + d3) + d2) + d1;
        x -= f / df;
        if (x - xold).abs() <= EPS * x.abs() {
            return Ok(x);
        }
    }
    Err(RocheError::Generic("xl3: exceeded maximum iterations".into()))
}

/// L1 with asynchronous rotation of primary.
pub fn xl11(q: f64, spin: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("xl11: q={} <= 0", q)));
    }
    let ssq = spin * spin;
    let mu = q / (1.0 + q);
    let a1 = -1.0 + mu;
    let a2 = 2.0 - 2.0 * mu;
    let a3 = -1.0 + mu;
    let a4 = ssq + 2.0 * mu;
    let a5 = -2.0 * ssq - mu;
    let a6 = ssq;
    let (d1, d2, d3, d4, d5) = (a2, 2.0 * a3, 3.0 * a4, 4.0 * a5, 5.0 * a6);

    let mut x = 1.0 / (1.0 + q);
    for _ in 0..NMAX {
        let xold = x;
        let f = x * (x * (x * (x * (x * a6 + a5) + a4) + a3) + a2) + a1;
        let df = x * (x * (x * (x * d5 + d4) + d3) + d2) + d1;
        x -= f / df;
        if (x - xold).abs() <= EPS * x.abs() {
            return Ok(x);
        }
    }
    Err(RocheError::Generic("xl11: exceeded maximum iterations".into()))
}

/// L1 with asynchronous rotation of secondary.
pub fn xl12(q: f64, spin: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("xl12: q={} <= 0", q)));
    }
    let ssq = spin * spin;
    let mu = q / (1.0 + q);
    let a1 = -1.0 + mu;
    let a2 = 2.0 - 2.0 * mu;
    let a3 = -ssq + mu;
    let a4 = 3.0 * ssq + 2.0 * mu - 2.0;
    let a5 = 1.0 - mu - 3.0 * ssq;
    let a6 = ssq;
    let (d1, d2, d3, d4, d5) = (a2, 2.0 * a3, 3.0 * a4, 4.0 * a5, 5.0 * a6);

    let mut x = 1.0 / (1.0 + q);
    for _ in 0..NMAX {
        let xold = x;
        let f = x * (x * (x * (x * (x * a6 + a5) + a4) + a3) + a2) + a1;
        let df = x * (x * (x * (x * d5 + d4) + d3) + d2) + d1;
        x -= f / df;
        if (x - xold).abs() <= EPS * x.abs() {
            return Ok(x);
        }
    }
    Err(RocheError::Generic("xl12: exceeded maximum iterations".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xl1_equal_mass() {
        let x = xl1(1.0).unwrap();
        assert!((x - 0.5).abs() < 1e-10, "xl1(1.0)={}", x);
    }

    #[test]
    fn test_xl1_range() {
        let x = xl1(0.1).unwrap();
        assert!(x > 0.0 && x < 1.0);
    }

    #[test]
    fn test_xl2_beyond_secondary() {
        let x = xl2(1.0).unwrap();
        assert!(x > 1.0);
    }

    #[test]
    fn test_xl3_opposite_side() {
        let x = xl3(1.0).unwrap();
        assert!(x < 0.0);
    }

    #[test]
    fn test_xl11_sync_equals_xl1() {
        let x1 = xl1(0.5).unwrap();
        let x11 = xl11(0.5, 1.0).unwrap();
        assert!((x1 - x11).abs() < 1e-10);
    }

    #[test]
    fn test_xl12_sync_equals_xl1() {
        let x1 = xl1(0.5).unwrap();
        let x12 = xl12(0.5, 1.0).unwrap();
        assert!((x1 - x12).abs() < 1e-10);
    }
}
