use crate::SubsError;

const MAXIT: usize = 100;

/// Find root of a function using combined Newton-Raphson / bisection.
///
/// `func(x)` must return `(f, df)` — function value and derivative.
/// The root must be bracketed: `f(x1)` and `f(x2)` must have opposite signs.
///
/// * `x1` - left bracket
/// * `x2` - right bracket
/// * `xacc` - required accuracy in x
pub fn rtsafe<F>(func: &F, x1: f64, x2: f64, xacc: f64) -> Result<f64, SubsError>
where
    F: Fn(f64) -> (f64, f64),
{
    let (fl, _) = func(x1);
    let (fh, _) = func(x2);

    if (fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0) {
        return Err(SubsError::RootNotBracketed { x1, x2, f1: fl, f2: fh });
    }

    if fl == 0.0 {
        return Ok(x1);
    }
    if fh == 0.0 {
        return Ok(x2);
    }

    let (mut xl, mut xh) = if fl < 0.0 { (x1, x2) } else { (x2, x1) };

    let mut rts = 0.5 * (x1 + x2);
    let mut dxold = (x2 - x1).abs();
    let mut dx = dxold;
    let (mut f, mut df) = func(rts);

    for _ in 0..MAXIT {
        // Use bisection if Newton step would leave bracket or converge too slowly
        if ((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0
            || (2.0 * f).abs() > (dxold * df).abs()
        {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if xl == rts {
                return Ok(rts);
            }
        } else {
            dxold = dx;
            dx = f / df;
            let temp = rts;
            rts -= dx;
            if temp == rts {
                return Ok(rts);
            }
        }

        if dx.abs() < xacc {
            return Ok(rts);
        }

        let result = func(rts);
        f = result.0;
        df = result.1;

        if f < 0.0 {
            xl = rts;
        } else {
            xh = rts;
        }
    }

    Err(SubsError::TooManyIterations("rtsafe".to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_sqrt2() {
        // f(x) = x^2 - 2, f'(x) = 2x
        let root = rtsafe(&|x| (x * x - 2.0, 2.0 * x), 1.0, 2.0, 1e-12).unwrap();
        assert!((root - std::f64::consts::SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn test_find_zero() {
        // f(x) = x, f'(x) = 1
        let root = rtsafe(&|x| (x, 1.0), -1.0, 1.0, 1e-12).unwrap();
        assert!(root.abs() < 1e-12);
    }

    #[test]
    fn test_not_bracketed() {
        let result = rtsafe(&|x| (x * x + 1.0, 2.0 * x), 1.0, 2.0, 1e-12);
        assert!(result.is_err());
    }
}
