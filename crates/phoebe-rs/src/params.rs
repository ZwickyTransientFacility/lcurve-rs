//! Binary system parameters for EB light curve synthesis.

use crate::EBType;
use crate::passband::Passband;

/// Parameters defining an eclipsing binary system.
#[derive(Debug, Clone)]
pub struct EBParams {
    /// Mass ratio q = M2/M1 (0 < q <= 1 by convention)
    pub q: f64,
    /// Orbital inclination in degrees (90 = edge-on)
    pub inclination: f64,
    /// Orbital period in days
    pub period: f64,
    /// Effective temperature of primary (K)
    pub t_eff1: f64,
    /// Effective temperature of secondary (K)
    pub t_eff2: f64,
    /// Roche lobe filling factor for primary (0-1, 1 = Roche-filling)
    pub fillout1: f64,
    /// Roche lobe filling factor for secondary (0-1)
    pub fillout2: f64,
    /// System type
    pub eb_type: EBType,
    /// Limb darkening coefficient for primary (linear law)
    pub ld1: f64,
    /// Limb darkening coefficient for secondary
    pub ld2: f64,
    /// Third light fraction (0 = no third light)
    pub l3: f64,
    /// Phase zero-point offset
    pub phi0: f64,
    /// Number of latitude strips for surface grid (higher = more accurate)
    pub n_grid: usize,
    /// Passband for flux computation
    pub passband: Passband,
}

impl Default for EBParams {
    fn default() -> Self {
        Self {
            q: 0.5,
            inclination: 85.0,
            period: 0.35,
            t_eff1: 6000.0,
            t_eff2: 5500.0,
            fillout1: 1.0,
            fillout2: 1.0,
            eb_type: EBType::Contact,
            ld1: 0.5,
            ld2: 0.5,
            l3: 0.0,
            phi0: 0.0,
            passband: Passband::Bolometric,
            n_grid: 40,
        }
    }
}

impl EBParams {
    /// Create a contact binary (W UMa) with given mass ratio and inclination.
    pub fn contact(q: f64, inclination: f64) -> Self {
        Self {
            q,
            inclination,
            fillout1: 1.0,
            fillout2: 1.0,
            eb_type: EBType::Contact,
            ..Default::default()
        }
    }

    /// Create a detached binary with given parameters.
    pub fn detached(q: f64, inclination: f64, r1_frac: f64, r2_frac: f64) -> Self {
        Self {
            q,
            inclination,
            fillout1: r1_frac,
            fillout2: r2_frac,
            eb_type: EBType::Detached,
            ..Default::default()
        }
    }

    /// Validate parameters.
    pub fn validate(&self) -> Result<(), crate::PhoebeError> {
        if self.q <= 0.0 || self.q > 10.0 {
            return Err(crate::PhoebeError::InvalidParams(
                format!("q={} out of range (0, 10]", self.q)));
        }
        if self.inclination < 0.0 || self.inclination > 90.0 {
            return Err(crate::PhoebeError::InvalidParams(
                format!("inclination={} out of range [0, 90]", self.inclination)));
        }
        if self.fillout1 < 0.0 || self.fillout1 > 1.5 {
            return Err(crate::PhoebeError::InvalidParams(
                format!("fillout1={} out of range [0, 1.5]", self.fillout1)));
        }
        Ok(())
    }
}
