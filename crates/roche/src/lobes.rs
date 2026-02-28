use crate::RocheError;

/// Eggleton's formula for volume-averaged Roche lobe radius (R_L/a).
pub fn rlobe_eggleton(q: f64) -> Result<f64, RocheError> {
    if q <= 0.0 {
        return Err(RocheError::Generic(format!("rlobe_eggleton: q={} <= 0", q)));
    }
    let q3 = q.powf(1.0 / 3.0);
    Ok(0.49 * q3 * q3 / (0.6 * q3 * q3 + (1.0 + q3).ln()))
}

/// d log(R_L) / d log(M2) assuming M1+M2 = constant.
pub fn zeta_rlobe_eggleton(q: f64) -> f64 {
    let q1 = q.powf(1.0 / 3.0);
    let loneq = (1.0 + q1).ln();
    (1.0 + q) / 3.0 * (2.0 * loneq - q1 / (1.0 + q1)) / (0.6 * q1 * q1 + loneq)
}

/// d zeta / d q.
pub fn dzetadq_rlobe_eggleton(q: f64) -> f64 {
    let q1 = q.powf(1.0 / 3.0);
    let q2 = q1 * q1;
    let opq1 = 1.0 + q1;
    let loneq = opq1.ln();
    let denom = 0.6 * q2 + loneq;
    let numer = 2.0 * loneq - q1 / opq1;
    numer / denom / 3.0
        + (1.0 + q) / 3.0
            * ((1.0 + 2.0 * q1) / 3.0 / (q1 * opq1).powi(2)
                - numer * (0.4 / q1 + 1.0 / (3.0 * q2 * opq1)) / denom)
            / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eggleton_q1() {
        let r = rlobe_eggleton(1.0).unwrap();
        assert!((r - 0.3789).abs() < 0.001, "rlobe(1.0)={}", r);
    }

    #[test]
    fn test_zeta_finite() {
        // zeta = d log(rl) / d log(m2) with m1+m2=const
        // If M = m1+m2 = const and q = m2/m1, then m2 = M*q/(1+q)
        // d log(rl)/d log(m2) = (m2/rl) * drl/dm2
        // = (m2/rl) * (drl/dq) * (dq/dm2) = (q/rl) * drl/dq * (1+q)^2/M * M/(1+q) = (q/rl)*drl/dq*(1+q)
        let q = 0.5;
        let dq = 1e-7;
        let rl_lo = rlobe_eggleton(q - dq).unwrap();
        let rl_hi = rlobe_eggleton(q + dq).unwrap();
        let rl = rlobe_eggleton(q).unwrap();
        let drl_dq = (rl_hi - rl_lo) / (2.0 * dq);
        let zeta_fd = q * (1.0 + q) / rl * drl_dq;
        let zeta = zeta_rlobe_eggleton(q);
        assert!((zeta_fd - zeta).abs() < 1e-4, "fd={}, zeta={}", zeta_fd, zeta);
    }
}
