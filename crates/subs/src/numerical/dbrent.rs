use crate::{sign, SubsError};

const ITMAX: usize = 100;

/// Brent's method for 1D minimization with derivatives.
///
/// Given a bracket `[ax, cx]` containing a minimum with `bx` as the
/// current best estimate, refines the minimum.
///
/// * `ax` - one extreme of x
/// * `bx` - mid x value (best guess)
/// * `cx` - other extreme of x
/// * `func` - function to minimize: `f(x) -> f64`
/// * `dfunc` - derivative function: `f(x) -> f64`
/// * `acc` - absolute accuracy in x
/// * `stopfast` - if true, return early when function drops below `fref`
/// * `fref` - reference function value for early stopping
///
/// Returns `(fmin, xmin)` — the function value and position at the minimum.
pub fn dbrent<F, D>(
    ax: f64,
    bx: f64,
    cx: f64,
    func: &mut F,
    dfunc: &mut D,
    acc: f64,
    stopfast: bool,
    fref: f64,
) -> Result<(f64, f64), SubsError>
where
    F: FnMut(f64) -> f64,
    D: FnMut(f64) -> f64,
{
    let mut a = ax.min(cx);
    let mut b = ax.max(cx);

    let mut x = bx;
    let mut w = bx;
    let mut v = bx;
    let mut fx = func(x);
    let mut fw = fx;
    let mut fv = fx;

    if stopfast && fx < fref {
        return Ok((fx, x));
    }

    let mut dx = dfunc(x);
    let mut dw = dx;
    let mut dv = dx;

    let mut e = 0.0f64;
    let mut d = 0.0f64;

    for _ in 0..ITMAX {
        let xm = 0.5 * (a + b);
        let tol1 = acc;
        let tol2 = 2.0 * tol1;

        if (x - xm).abs() <= tol2 - 0.5 * (b - a) {
            return Ok((fx, x));
        }

        if e.abs() > tol1 {
            let mut d1 = 2.0 * (b - a);
            let mut d2 = d1;
            if dw != dx {
                d1 = (w - x) * dx / (dx - dw);
            }
            if dv != dx {
                d2 = (v - x) * dx / (dx - dv);
            }

            let u1 = x + d1;
            let u2 = x + d2;
            let ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
            let ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
            let olde = e;
            e = d;

            if ok1 || ok2 {
                d = if ok1 && ok2 {
                    if d1.abs() < d2.abs() { d1 } else { d2 }
                } else if ok1 {
                    d1
                } else {
                    d2
                };

                if d.abs() <= (0.5 * olde).abs() {
                    let u = x + d;
                    if u - a < tol2 || b - u < tol2 {
                        d = sign(tol1, xm - x);
                    }
                } else {
                    e = if dx >= 0.0 { a - x } else { b - x };
                    d = 0.5 * e;
                }
            } else {
                e = if dx >= 0.0 { a - x } else { b - x };
                d = 0.5 * e;
            }
        } else {
            e = if dx >= 0.0 { a - x } else { b - x };
            d = 0.5 * e;
        }

        let (u, fu);
        if d.abs() >= tol1 {
            u = x + d;
            fu = func(u);
            if stopfast && fu < fref {
                return Ok((fu, u));
            }
        } else {
            u = x + sign(tol1, d);
            fu = func(u);
            if stopfast && fu < fref {
                return Ok((fu, u));
            }
            if fu > fx {
                return Ok((fx, x));
            }
        }

        let du = dfunc(u);

        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            dv = dw;
            w = x;
            fw = fx;
            dw = dx;
            x = u;
            fx = fu;
            dx = du;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                dv = dw;
                w = u;
                fw = fu;
                dw = du;
            } else if fu < fv || v == x || v == w {
                v = u;
                fv = fu;
                dv = du;
            }
        }
    }

    Err(SubsError::TooManyIterations("dbrent".to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quadratic_min() {
        // f(x) = (x-3)^2 + 1, min at x=3
        let mut func = |x: f64| (x - 3.0) * (x - 3.0) + 1.0;
        let mut dfunc = |x: f64| 2.0 * (x - 3.0);

        let (fmin, xmin) = dbrent(0.0, 2.0, 6.0, &mut func, &mut dfunc, 1e-10, false, 0.0).unwrap();

        assert!((xmin - 3.0).abs() < 1e-9, "xmin={}", xmin);
        assert!((fmin - 1.0).abs() < 1e-9, "fmin={}", fmin);
    }
}
