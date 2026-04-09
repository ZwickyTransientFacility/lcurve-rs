"""
Eclipsing binary parameter estimation using phoebe-rs + Eryn MCMC.

Fits mass ratio (q), inclination (i), temperature ratio, and fill factors
to observed eclipsing binary light curves.

Example
-------
>>> from lcurve_rs.eb_fitting import EBFitter
>>> fitter = EBFitter(phases, mags, errs)
>>> result = fitter.run_eryn(nwalkers=24, nsteps=3000, burn=500)
>>> print(f"q = {result.median('q'):.3f} +/- {result.std('q'):.3f}")
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# numpy 2.x compat: eryn uses np.in1d which was removed
if not hasattr(np, 'in1d'):
    np.in1d = np.isin

from lcurve_rs import EBParams, eb_lightcurve


@dataclass
class EBFitResult:
    """Result of EB parameter estimation."""

    param_names: list
    samples: np.ndarray        # shape (n_samples, n_params)
    log_likelihood: np.ndarray  # shape (n_samples,)
    raw_result: object = None

    def median(self, name: str) -> float:
        idx = self.param_names.index(name)
        return float(np.median(self.samples[:, idx]))

    def std(self, name: str) -> float:
        idx = self.param_names.index(name)
        return float(np.std(self.samples[:, idx]))

    def percentile(self, name: str, q: float) -> float:
        idx = self.param_names.index(name)
        return float(np.percentile(self.samples[:, idx], q))

    def summary(self) -> str:
        lines = []
        for i, name in enumerate(self.param_names):
            med = np.median(self.samples[:, i])
            lo = np.percentile(self.samples[:, i], 16)
            hi = np.percentile(self.samples[:, i], 84)
            lines.append(f"  {name:>12s} = {med:.4f}  (+{hi-med:.4f} / -{med-lo:.4f})")
        return "\n".join(lines)

    def best_params(self) -> dict:
        """Return parameter values at maximum likelihood."""
        idx = np.argmax(self.log_likelihood)
        return {name: float(self.samples[idx, i])
                for i, name in enumerate(self.param_names)}


class EBFitter:
    """Fit eclipsing binary parameters to an observed light curve.

    Parameters
    ----------
    phases : array
        Orbital phases (0 to 1).
    mags : array
        Observed magnitudes (relative or absolute).
    errs : array
        Magnitude uncertainties.
    eb_type : str
        'contact' or 'detached'.
    fit_params : list of str
        Which parameters to fit. Default: ['q', 'inclination', 'tratio'].
    n_grid : int
        Surface mesh resolution (default 25 for MCMC speed).
    """

    # Default parameter bounds
    PARAM_BOUNDS = {
        'q':           (0.05, 1.0),
        'inclination': (30.0, 90.0),
        'tratio':      (0.5, 1.0),   # T2/T1
        'fillout1':    (0.5, 1.0),
        'fillout2':    (0.5, 1.0),
        'l3':          (0.0, 0.5),
        'phi0':        (-0.05, 0.05),
    }

    def __init__(
        self,
        phases,
        mags,
        errs,
        eb_type: str = 'contact',
        fit_params: Optional[list] = None,
        n_grid: int = 25,
        band_data: Optional[dict] = None,
    ):
        """Initialize EB fitter.

        For single-band fitting, pass phases/mags/errs as arrays.
        For multi-band fitting, pass band_data as a dict:
            {band: (phases, mags, errs)} where band is 'g', 'r', or 'i'.
        If band_data is provided, phases/mags/errs are ignored.
        """
        self.eb_type = eb_type
        self.n_grid = n_grid

        if fit_params is None:
            fit_params = ['q', 'inclination', 'tratio']
        self.fit_params = fit_params
        self.ndim = len(fit_params)

        if band_data is not None:
            self.multiband = True
            self.band_data = {}
            for band, (p, m, e) in band_data.items():
                p = np.asarray(p, dtype=np.float64)
                m = np.asarray(m, dtype=np.float64)
                e = np.maximum(np.asarray(e, dtype=np.float64), 0.005)
                self.band_data[band] = (p, m, e, 1.0 / e**2)
        else:
            self.multiband = False
            self.phases = np.asarray(phases, dtype=np.float64)
            self.mags = np.asarray(mags, dtype=np.float64)
            self.errs = np.maximum(np.asarray(errs, dtype=np.float64), 0.005)
            self.w = 1.0 / self.errs**2

    def _build_params(self, theta: np.ndarray, passband: str = 'r') -> EBParams:
        """Construct EBParams from parameter vector."""
        values = {name: theta[i] for i, name in enumerate(self.fit_params)}

        q = values.get('q', 0.5)
        incl = values.get('inclination', 80.0)
        tratio = values.get('tratio', 0.9)
        t1 = 6000.0
        t2 = t1 * tratio

        if self.eb_type == 'contact':
            params = EBParams.contact(q, incl, t1=t1, t2=t2,
                                      n_grid=self.n_grid, passband=passband)
        else:
            f1 = values.get('fillout1', 0.5)
            f2 = values.get('fillout2', 0.5)
            params = EBParams.detached(q, incl, f1, f2, t1=t1, t2=t2,
                                       n_grid=self.n_grid, passband=passband)
        return params

    def _band_chi2(self, theta, phases, mags, w, passband='r'):
        """Compute chi2 for a single band."""
        try:
            params = self._build_params(theta, passband=passband)
            lc = eb_lightcurve(params, phases)
            model_mag = -2.5 * np.log10(np.array(lc['flux']))
        except Exception:
            return np.inf

        A = np.column_stack([model_mag, np.ones(len(phases))])
        Aw = A * w[:, None]
        try:
            coeffs = np.linalg.lstsq(Aw.T @ A, Aw.T @ mags, rcond=None)[0]
        except np.linalg.LinAlgError:
            return np.inf

        fitted = A @ coeffs
        return float(np.sum(w * (mags - fitted)**2))

    def log_likelihood(self, theta: np.ndarray) -> float:
        """Compute log-likelihood for a single parameter vector."""
        # Check bounds
        for i, name in enumerate(self.fit_params):
            lo, hi = self.PARAM_BOUNDS[name]
            if theta[i] < lo or theta[i] > hi:
                return -np.inf

        if self.multiband:
            chi2 = 0.0
            for band, (phases, mags, errs, w) in self.band_data.items():
                c = self._band_chi2(theta, phases, mags, w, passband=band)
                if not np.isfinite(c):
                    return -np.inf
                chi2 += c
            return -0.5 * chi2
        else:
            c = self._band_chi2(theta, self.phases, self.mags, self.w)
            if not np.isfinite(c):
                return -np.inf
            return -0.5 * c

    def log_likelihood_batch(self, theta_batch: np.ndarray) -> np.ndarray:
        """Vectorized log-likelihood for multiple parameter vectors.

        theta_batch: shape (n_walkers, 1, n_params) for eryn
        """
        if theta_batch.ndim == 3:
            # eryn format: (nwalkers, 1, ndim)
            theta_batch = theta_batch[:, 0, :]

        return np.array([self.log_likelihood(t) for t in theta_batch])

    def run_eryn(
        self,
        nwalkers: int = 0,
        nsteps: int = 3000,
        burn: int = 500,
        ntemps: int = 5,
        thin_by: int = 1,
        progress: bool = True,
        **kwargs,
    ) -> EBFitResult:
        """Run Eryn parallel-tempered MCMC.

        Parameters
        ----------
        nwalkers : int
            Number of walkers (default: 2*ndim + 2, minimum 24).
        nsteps : int
            Total steps after burn-in.
        burn : int
            Burn-in steps.
        ntemps : int
            Number of temperatures (default 5).

        Returns
        -------
        EBFitResult
        """
        from eryn.ensemble import EnsembleSampler
        from eryn.prior import ProbDistContainer

        if nwalkers <= 0:
            nwalkers = max(24, 2 * self.ndim + 2)

        # Build eryn priors
        class UniformAdapter:
            def __init__(self, lo, hi):
                self.lo, self.hi = lo, hi
                self.span = hi - lo
            def rvs(self, size=1):
                return np.random.uniform(self.lo, self.hi, size=size)
            def logpdf(self, x):
                x = np.asarray(x)
                out = np.where((x >= self.lo) & (x <= self.hi),
                               -np.log(self.span), -np.inf)
                return out

        eryn_priors = ProbDistContainer(
            {i: UniformAdapter(*self.PARAM_BOUNDS[name])
             for i, name in enumerate(self.fit_params)}
        )

        tempering_kwargs = {"ntemps": ntemps} if ntemps > 1 else {}

        # Initial positions
        actual_ntemps = ntemps if ntemps > 1 else 1
        coords = eryn_priors.rvs(size=(actual_ntemps, nwalkers, 1))

        sampler = EnsembleSampler(
            nwalkers,
            self.ndim,
            self.log_likelihood_batch,
            eryn_priors,
            tempering_kwargs=tempering_kwargs,
            vectorize=True,
            **kwargs,
        )
        sampler.run_mcmc(coords, nsteps, burn=burn, progress=progress,
                         thin_by=thin_by)

        # Extract cold chain
        chain = sampler.get_chain(temp_index=0)["model_0"]
        chain = chain.reshape(-1, self.ndim)
        logl = sampler.get_log_like(temp_index=0).reshape(-1)

        return EBFitResult(
            param_names=list(self.fit_params),
            samples=chain,
            log_likelihood=logl,
            raw_result=sampler,
        )
