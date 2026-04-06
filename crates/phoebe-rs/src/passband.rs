//! Passband-dependent surface brightness and limb darkening.
//!
//! Instead of bolometric T^4, computes bandpass-integrated Planck function
//! for ZTF g, ZTF r, Johnson V, or bolometric.

use std::f64::consts::PI;

/// Available passbands.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Passband {
    Bolometric,
    ZtfG,
    ZtfR,
    JohnsonV,
}

impl Passband {
    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "bolometric" | "bol" => Some(Passband::Bolometric),
            "ztf_g" | "ztf-g" | "g" => Some(Passband::ZtfG),
            "ztf_r" | "ztf-r" | "r" => Some(Passband::ZtfR),
            "johnson_v" | "johnson-v" | "v" => Some(Passband::JohnsonV),
            _ => None,
        }
    }

    /// Default linear limb darkening coefficient for a solar-type star.
    pub fn default_ld(&self) -> f64 {
        match self {
            Passband::Bolometric => 0.50,
            Passband::ZtfG => 0.65,
            Passband::ZtfR => 0.45,
            Passband::JohnsonV => 0.55,
        }
    }
}

/// Planck function B(λ, T) in SI units.
/// λ in metres, T in Kelvin.
fn planck(lam: f64, t: f64) -> f64 {
    const H: f64 = 6.626e-34;
    const C: f64 = 3.0e8;
    const K: f64 = 1.381e-23;

    let x = H * C / (lam * K * t);
    if x > 500.0 { return 0.0; }  // avoid overflow
    2.0 * H * C * C / (lam.powi(5)) / (x.exp() - 1.0)
}

/// Compute bandpass-integrated surface brightness for a given temperature.
///
/// Returns a value proportional to the flux per unit area in the passband.
/// The absolute scale doesn't matter since we normalise light curves.
pub fn surface_brightness(t_eff: f64, passband: Passband) -> f64 {
    match passband {
        Passband::Bolometric => {
            // Stefan-Boltzmann: sb ~ T^4
            t_eff.powi(4)
        },
        _ => {
            // Numerical integration of Planck × transmission
            let (lam_min, lam_max, n_pts) = (300e-9, 1100e-9, 200);
            let dlam = (lam_max - lam_min) / n_pts as f64;

            let mut flux = 0.0;
            for i in 0..n_pts {
                let lam = lam_min + (i as f64 + 0.5) * dlam;
                let lam_nm = lam * 1e9;
                let trans = transmission(lam_nm, passband);
                flux += planck(lam, t_eff) * trans * dlam;
            }
            flux
        }
    }
}

/// Approximate transmission function for a passband.
/// Input: wavelength in nm.
fn transmission(lam_nm: f64, passband: Passband) -> f64 {
    match passband {
        Passband::ZtfG => {
            // ZTF g-band: ~400-560 nm, peak ~480 nm
            if lam_nm < 400.0 || lam_nm > 560.0 { return 0.0; }
            (-0.5 * ((lam_nm - 480.0) / 40.0).powi(2)).exp()
        },
        Passband::ZtfR => {
            // ZTF r-band: ~560-730 nm, peak ~640 nm
            if lam_nm < 560.0 || lam_nm > 730.0 { return 0.0; }
            (-0.5 * ((lam_nm - 640.0) / 45.0).powi(2)).exp()
        },
        Passband::JohnsonV => {
            // Johnson V: ~480-640 nm, peak ~550 nm
            if lam_nm < 480.0 || lam_nm > 640.0 { return 0.0; }
            (-0.5 * ((lam_nm - 550.0) / 40.0).powi(2)).exp()
        },
        Passband::Bolometric => 1.0,
    }
}

/// Precomputed surface brightness lookup table for fast evaluation.
///
/// Stores SB values for temperatures 2000-10000 K in steps of 50 K.
pub struct SBTable {
    values: Vec<f64>,
    t_min: f64,
    dt: f64,
}

impl SBTable {
    /// Build lookup table for a passband.
    pub fn new(passband: Passband) -> Self {
        let t_min = 2000.0;
        let t_max = 10000.0;
        let dt = 50.0;
        let n = ((t_max - t_min) / dt) as usize + 1;

        let values: Vec<f64> = (0..n)
            .map(|i| surface_brightness(t_min + i as f64 * dt, passband))
            .collect();

        SBTable { values, t_min, dt }
    }

    /// Interpolate surface brightness at a given temperature.
    pub fn eval(&self, t_eff: f64) -> f64 {
        let idx_f = (t_eff - self.t_min) / self.dt;
        let idx = idx_f.floor() as usize;
        if idx >= self.values.len() - 1 {
            return *self.values.last().unwrap();
        }
        if idx_f < 0.0 {
            return self.values[0];
        }
        let frac = idx_f - idx as f64;
        self.values[idx] * (1.0 - frac) + self.values[idx + 1] * frac
    }
}
