use crate::constants;
use crate::sqr;

/// Factor: 2e27 * h * c  (wavelength in nm → frequency domain)
const FAC1: f64 = 2.0e27 * constants::H * constants::C;

/// Factor: 1e9 * h * c / k  (for exponential argument)
const FAC2: f64 = 1.0e9 * constants::H * constants::C / constants::K;

/// Planck function Bnu = (2 h nu^3 / c^2) / (exp(h nu/kT) - 1).
///
/// Output units: W/m^2/Hz/sr.
///
/// * `wave` - wavelength in nanometres
/// * `temp` - temperature in Kelvin
pub fn planck(wave: f64, temp: f64) -> f64 {
    let efac = FAC2 / (wave * temp);
    if efac > 40.0 {
        FAC1 * (-efac).exp() / (wave * sqr(wave))
    } else {
        FAC1 / (efac.exp() - 1.0) / (wave * sqr(wave))
    }
}

/// Logarithmic derivative of Planck function wrt wavelength:
/// d ln(Bnu) / d ln(lambda).
///
/// * `wave` - wavelength in nanometres
/// * `temp` - temperature in Kelvin
pub fn dplanck(wave: f64, temp: f64) -> f64 {
    let efac = FAC2 / (wave * temp);
    efac / (1.0 - (-efac).exp()) - 3.0
}

/// Logarithmic derivative of Planck function wrt temperature:
/// d ln(Bnu) / d ln(T).
///
/// * `wave` - wavelength in nanometres
/// * `temp` - temperature in Kelvin
pub fn dlpdlt(wave: f64, temp: f64) -> f64 {
    let efac = FAC2 / (wave * temp);
    efac / (1.0 - (-efac).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planck_positive() {
        // Planck function should be positive for any wave/temp
        let f = planck(500.0, 10000.0);
        assert!(f > 0.0);

        let f2 = planck(500.0, 5000.0);
        assert!(f2 > 0.0);

        // Hotter should be brighter at same wavelength (for Bnu)
        assert!(f > f2);
    }

    #[test]
    fn test_planck_high_temp_approx() {
        // For very high efac (low T, short wave), should still be positive
        let f = planck(100.0, 100.0);
        assert!(f >= 0.0);
    }

    #[test]
    fn test_dlpdlt_positive() {
        // dln(B)/dln(T) should always be positive (hotter = brighter)
        assert!(dlpdlt(500.0, 10000.0) > 0.0);
        assert!(dlpdlt(500.0, 3000.0) > 0.0);
    }
}
