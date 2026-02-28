use crate::sqr;

/// Result of a Bulirsch-Stoer step.
pub struct BsStepResult {
    /// The step size that was actually completed.
    pub hdid: f64,
    /// The suggested next step size.
    pub hnext: f64,
    /// True if step size underflow (x+h == x).
    pub underflow: bool,
}

/// Persistent state for Bulirsch-Stoer integrator.
/// In the C++ code these are `static` variables; here we store them in a struct.
pub struct BsState {
    first: bool,
    kmax: usize,
    kopt: usize,
    epsold: f64,
    xnew: f64,
    a: [f64; 9],       // IMAXX = KMAXX+1 = 9
    alf: [[f64; 8]; 8], // KMAXX × KMAXX = 8×8
}

const KMAXX: usize = 8;
const IMAXX: usize = KMAXX + 1;
const NSEQ: [usize; IMAXX] = [2, 4, 6, 8, 10, 12, 14, 16, 18];

const SAFE1: f64 = 0.25;
const SAFE2: f64 = 0.7;
const REDMAX: f64 = 1.0e-5;
const REDMIN: f64 = 0.7;
const TINY: f64 = 1.0e-30;
const SCALMX: f64 = 0.1;

impl BsState {
    pub fn new() -> Self {
        BsState {
            first: true,
            kmax: 0,
            kopt: 0,
            epsold: -1.0,
            xnew: -1.0e29,
            a: [0.0; IMAXX],
            alf: [[0.0; KMAXX]; KMAXX],
        }
    }
}

impl Default for BsState {
    fn default() -> Self {
        Self::new()
    }
}

/// Modified midpoint method: advance y from xs to xs+htot using nstep substeps.
fn mmid(
    y: &[f64],
    dydx: &[f64],
    nv: usize,
    xs: f64,
    htot: f64,
    nstep: usize,
    yout: &mut [f64],
    derivs: &dyn Fn(f64, &[f64], &mut [f64]),
) {
    let h = htot / nstep as f64;
    let mut ym = vec![0.0; nv];
    let mut yn = vec![0.0; nv];

    for i in 0..nv {
        ym[i] = y[i];
        yn[i] = y[i] + h * dydx[i];
    }

    let mut x = xs + h;
    derivs(x, &yn, yout);

    let h2 = 2.0 * h;
    for _ in 1..nstep {
        for i in 0..nv {
            let swap = ym[i] + h2 * yout[i];
            ym[i] = yn[i];
            yn[i] = swap;
        }
        x += h;
        derivs(x, &yn, yout);
    }

    for i in 0..nv {
        yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i]);
    }
}

/// Polynomial extrapolation used by bsstep.
fn pzextr(
    iest: usize,
    xest: f64,
    yest: &[f64],
    yz: &mut [f64],
    dy: &mut [f64],
    nv: usize,
    d_tab: &mut Vec<Vec<f64>>,
    x_tab: &mut Vec<f64>,
) {
    // Ensure storage
    if x_tab.len() <= iest {
        x_tab.resize(iest + 1, 0.0);
    }
    for row in d_tab.iter_mut() {
        if row.len() <= iest {
            row.resize(iest + 1, 0.0);
        }
    }

    x_tab[iest] = xest;

    for j in 0..nv {
        dy[j] = yest[j];
        yz[j] = yest[j];
    }

    if iest == 0 {
        for j in 0..nv {
            d_tab[j][0] = yest[j];
        }
    } else {
        let mut c = vec![0.0; nv];
        for j in 0..nv {
            c[j] = yest[j];
        }
        for k1 in 0..iest {
            let delta = 1.0 / (x_tab[iest - k1 - 1] - xest);
            let f1 = xest * delta;
            let f2 = x_tab[iest - k1 - 1] * delta;
            for j in 0..nv {
                let q = d_tab[j][k1];
                d_tab[j][k1] = dy[j];
                let delta_val = c[j] - q;
                dy[j] = f1 * delta_val;
                c[j] = f2 * delta_val;
                yz[j] += dy[j];
            }
        }
        for j in 0..nv {
            d_tab[j][iest] = dy[j];
        }
    }
}

/// Bulirsch-Stoer adaptive ODE integration step.
///
/// Advances `y[0..nv]` from `xx` by approximately `htry`, maintaining accuracy `eps`.
///
/// `derivs(t, y, dydt)` computes the RHS of the ODE system.
///
/// Returns a `BsStepResult` with the actual step completed and suggested next step.
/// If `underflow` is true, the step size became too small (x+h == x).
pub fn bsstep(
    y: &mut [f64],
    dydx: &[f64],
    nv: usize,
    xx: &mut f64,
    htry: f64,
    eps: f64,
    yscal: &[f64],
    state: &mut BsState,
    derivs: &dyn Fn(f64, &[f64], &mut [f64]),
) -> BsStepResult {
    let mut yerr = vec![0.0; nv];
    let mut ysav = vec![0.0; nv];
    let mut yseq = vec![0.0; nv];
    let mut err = [0.0f64; KMAXX];

    // Extrapolation tableau storage
    let mut d_tab: Vec<Vec<f64>> = vec![vec![0.0; KMAXX]; nv];
    let mut x_tab: Vec<f64> = vec![0.0; KMAXX];

    // Recompute work coefficients if eps changed
    if eps != state.epsold {
        state.xnew = -1.0e29;
        let eps1 = SAFE1 * eps;
        state.a[0] = NSEQ[0] as f64 + 1.0;
        for k in 0..KMAXX {
            state.a[k + 1] = state.a[k] + NSEQ[k + 1] as f64;
        }
        for iq in 1..KMAXX {
            for k in 0..iq {
                state.alf[k][iq] = eps1.powf(
                    (state.a[k + 1] - state.a[iq + 1])
                        / ((state.a[iq + 1] - state.a[0] + 1.0) * (2 * k + 3) as f64),
                );
            }
        }
        state.epsold = eps;
        state.kopt = 1;
        for kopt in 1..KMAXX - 1 {
            if state.a[kopt + 1] > state.a[kopt] * state.alf[kopt - 1][kopt] {
                state.kopt = kopt;
                break;
            }
            state.kopt = kopt;
        }
        state.kmax = state.kopt;
    }

    let mut h = htry;
    for i in 0..nv {
        ysav[i] = y[i];
    }

    if *xx != state.xnew || h != state.xnew {
        // The C++ code compares h != hnext but hnext isn't stored separately.
        // We use the first flag to handle this.
        state.first = true;
        state.kopt = state.kmax;
    }

    let mut reduct = false;
    let mut km = 0usize;
    loop {
        let mut exitflag = false;
        for k in 0..=state.kmax {
            state.xnew = *xx + h;
            if state.xnew == *xx {
                return BsStepResult {
                    hdid: 0.0,
                    hnext: 0.0,
                    underflow: true,
                };
            }

            mmid(&ysav, dydx, nv, *xx, h, NSEQ[k], &mut yseq, derivs);
            let xest = sqr(h / NSEQ[k] as f64);
            pzextr(k, xest, &yseq, y, &mut yerr, nv, &mut d_tab, &mut x_tab);

            if k != 0 {
                let mut errmax = TINY;
                for i in 0..nv {
                    errmax = errmax.max((yerr[i] / yscal[i]).abs());
                }
                errmax /= eps;
                km = k - 1;
                err[km] = (errmax / SAFE1).powf(1.0 / (2 * km + 3) as f64);
            }

            if k != 0 && (k >= state.kopt.saturating_sub(1) || state.first) {
                if err[km] != 0.0 {
                    // errmax was computed above
                    let mut errmax = TINY;
                    for i in 0..nv {
                        errmax = errmax.max((yerr[i] / yscal[i]).abs());
                    }
                    errmax /= eps;

                    if errmax < 1.0 {
                        exitflag = true;
                        break;
                    }
                }

                if k == state.kmax || k == state.kopt + 1 {
                    break;
                } else if k == state.kopt && state.alf[state.kopt - 1][state.kopt] < err[km] {
                    break;
                } else if state.kopt == state.kmax
                    && state.alf[km][state.kmax - 1] < err[km]
                {
                    break;
                } else if state.alf[km][state.kopt] < err[km] {
                    break;
                }
            }
        }

        if exitflag {
            break;
        }

        // Reduce step
        let red = if err[km] > 0.0 {
            let r = SAFE2 / err[km];
            r.max(REDMAX).min(REDMIN)
        } else {
            REDMIN
        };
        h *= red;
        reduct = true;
    }

    *xx = state.xnew;
    let hdid = h;
    state.first = false;

    // Determine optimal order for next step
    let mut wrkmin = 1.0e35f64;
    let mut scale = 1.0;
    for kk in 0..=km {
        let fact = err[kk].max(SCALMX);
        let work = fact * state.a[kk + 1];
        if work < wrkmin {
            scale = fact;
            wrkmin = work;
            state.kopt = kk + 1;
        }
    }

    let mut hnext = h / scale;
    let k_final = km + 1; // k at break

    if state.kopt >= k_final && state.kopt != state.kmax && !reduct {
        let fact = (scale / state.alf[state.kopt - 1][state.kopt]).max(SCALMX);
        if state.a[state.kopt + 1] * fact <= wrkmin {
            hnext = h / fact;
            state.kopt += 1;
        }
    }

    BsStepResult {
        hdid,
        hnext,
        underflow: false,
    }
}

/// Stoermer's rule variant for 2nd order conservative ODEs.
fn stoerm(
    y: &[f64],
    d2y: &[f64],
    nvar: usize,
    xs: f64,
    htot: f64,
    nstep: usize,
    yout: &mut [f64],
    derivs: &dyn Fn(f64, &[f64], &mut [f64]),
) {
    let h = htot / nstep as f64;
    let halfh = 0.5 * h;
    let neqns = nvar / 2;
    let mut ytemp = vec![0.0; nvar];

    for i in 0..neqns {
        let n = neqns + i;
        ytemp[n] = h * (y[n] + halfh * d2y[i]);
        ytemp[i] = y[i] + ytemp[n];
    }

    let mut x = xs + h;
    derivs(x, &ytemp, yout);

    let h2 = h * h;
    for _ in 1..nstep {
        for i in 0..neqns {
            let n = neqns + i;
            ytemp[n] += h2 * yout[i];
            ytemp[i] += ytemp[n];
        }
        x += h;
        derivs(x, &ytemp, yout);
    }

    for i in 0..neqns {
        let n = neqns + i;
        yout[n] = ytemp[n] / h + halfh * yout[i];
        yout[i] = ytemp[i];
    }
}

/// Persistent state for Stoermer-variant Bulirsch-Stoer integrator.
pub struct BsStateSt {
    first: bool,
    kmax: usize,
    kopt: usize,
    epsold: f64,
    xnew: f64,
    a: [f64; 13],
    alf: [[f64; 12]; 12],
}

const KMAXX_ST: usize = 12;
const IMAXX_ST: usize = KMAXX_ST + 1;
const NSEQ_ST: [usize; IMAXX_ST] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

impl BsStateSt {
    pub fn new() -> Self {
        BsStateSt {
            first: true,
            kmax: 0,
            kopt: 0,
            epsold: -1.0,
            xnew: -1.0e29,
            a: [0.0; IMAXX_ST],
            alf: [[0.0; KMAXX_ST]; KMAXX_ST],
        }
    }
}

impl Default for BsStateSt {
    fn default() -> Self {
        Self::new()
    }
}

/// Bulirsch-Stoer step using Stoermer's rule for 2nd order conservative ODEs.
pub fn bsstepst(
    y: &mut [f64],
    dydx: &[f64],
    nv: usize,
    xx: &mut f64,
    htry: f64,
    eps: f64,
    yscal: &[f64],
    state: &mut BsStateSt,
    derivs: &dyn Fn(f64, &[f64], &mut [f64]),
) -> BsStepResult {
    let mut yerr = vec![0.0; nv];
    let mut ysav = vec![0.0; nv];
    let mut yseq = vec![0.0; nv];
    let mut err = [0.0f64; KMAXX_ST];
    let mut d_tab: Vec<Vec<f64>> = vec![vec![0.0; KMAXX_ST]; nv];
    let mut x_tab: Vec<f64> = vec![0.0; KMAXX_ST];

    if eps != state.epsold {
        state.xnew = -1.0e29;
        let eps1 = SAFE1 * eps;
        state.a[0] = NSEQ_ST[0] as f64 + 1.0;
        for k in 0..KMAXX_ST {
            state.a[k + 1] = state.a[k] + NSEQ_ST[k + 1] as f64;
        }
        for iq in 1..KMAXX_ST {
            for k in 0..iq {
                state.alf[k][iq] = eps1.powf(
                    (state.a[k + 1] - state.a[iq + 1])
                        / ((state.a[iq + 1] - state.a[0] + 1.0) * (2 * k + 3) as f64),
                );
            }
        }
        state.epsold = eps;
        state.kopt = 1;
        for kopt in 1..KMAXX_ST - 1 {
            if state.a[kopt + 1] > state.a[kopt] * state.alf[kopt - 1][kopt] {
                state.kopt = kopt;
                break;
            }
            state.kopt = kopt;
        }
        state.kmax = state.kopt;
    }

    let mut h = htry;
    for i in 0..nv {
        ysav[i] = y[i];
    }

    if *xx != state.xnew {
        state.first = true;
        state.kopt = state.kmax;
    }

    let mut reduct = false;
    let mut km = 0usize;

    loop {
        let mut exitflag = false;
        for k in 0..=state.kmax {
            state.xnew = *xx + h;
            if state.xnew == *xx {
                return BsStepResult {
                    hdid: 0.0,
                    hnext: 0.0,
                    underflow: true,
                };
            }

            stoerm(&ysav, dydx, nv, *xx, h, NSEQ_ST[k], &mut yseq, derivs);
            let xest = sqr(h / NSEQ_ST[k] as f64);
            pzextr(k, xest, &yseq, y, &mut yerr, nv, &mut d_tab, &mut x_tab);

            if k != 0 {
                let mut errmax = TINY;
                for i in 0..nv {
                    errmax = errmax.max((yerr[i] / yscal[i]).abs());
                }
                errmax /= eps;
                km = k - 1;
                err[km] = (errmax / SAFE1).powf(1.0 / (2 * km + 3) as f64);
            }

            if k != 0 && (k >= state.kopt.saturating_sub(1) || state.first) {
                if err[km] != 0.0 {
                    let mut errmax = TINY;
                    for i in 0..nv {
                        errmax = errmax.max((yerr[i] / yscal[i]).abs());
                    }
                    errmax /= eps;

                    if errmax < 1.0 {
                        exitflag = true;
                        break;
                    }
                }

                if k == state.kmax || k == state.kopt + 1 {
                    break;
                } else if k == state.kopt && state.alf[state.kopt - 1][state.kopt] < err[km] {
                    break;
                } else if state.kopt == state.kmax
                    && state.alf[km][state.kmax - 1] < err[km]
                {
                    break;
                } else if state.alf[km][state.kopt] < err[km] {
                    break;
                }
            }
        }

        if exitflag {
            break;
        }

        let red = if err[km] > 0.0 {
            let r = SAFE2 / err[km];
            r.max(REDMAX).min(REDMIN)
        } else {
            REDMIN
        };
        h *= red;
        reduct = true;
    }

    *xx = state.xnew;
    let hdid = h;
    state.first = false;

    let mut wrkmin = 1.0e35f64;
    let mut scale = 1.0;
    for kk in 0..=km {
        let fact = err[kk].max(SCALMX);
        let work = fact * state.a[kk + 1];
        if work < wrkmin {
            scale = fact;
            wrkmin = work;
            state.kopt = kk + 1;
        }
    }

    let mut hnext = h / scale;
    let k_final = km + 1;

    if state.kopt >= k_final && state.kopt != state.kmax && !reduct {
        let fact = (scale / state.alf[state.kopt - 1][state.kopt]).max(SCALMX);
        if state.a[state.kopt + 1] * fact <= wrkmin {
            hnext = h / fact;
            state.kopt += 1;
        }
    }

    BsStepResult {
        hdid,
        hnext,
        underflow: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_ode() {
        // dy/dx = y, y(0) = 1 => y(x) = e^x
        // Integrate from 0 to 1
        let nv = 1;
        let mut y = [1.0];
        let dydx = [1.0]; // initial derivative
        let mut xx = 0.0;
        let yscal = [1.0];
        let mut state = BsState::new();

        let derivs = |_t: f64, y: &[f64], dydt: &mut [f64]| {
            dydt[0] = y[0];
        };

        let result = bsstep(
            &mut y, &dydx, nv, &mut xx, 1.0, 1e-10, &yscal, &mut state, &derivs,
        );

        assert!(!result.underflow);
        // After integration, y should be close to e^(hdid)
        let expected = xx.exp();
        assert!(
            (y[0] - expected).abs() / expected < 1e-8,
            "y={}, expected={}, xx={}",
            y[0],
            expected,
            xx
        );
    }
}
