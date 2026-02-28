/// Offset between JD and MJD
pub const MJD2JD: f64 = 2_400_000.5;

/// Speed of light, m/s (exact)
pub const C: f64 = 2.997_924_58e8;

/// Gravitational constant, MKS
pub const G: f64 = 6.673e-11;

/// Planck's constant, MKS
pub const H: f64 = 6.6262e-34;

/// Boltzmann's constant, MKS
pub const K: f64 = 1.3806e-23;

/// Charge on the electron (magnitude, C)
pub const E_CHARGE: f64 = 1.602_176_565e-19;

/// Mass of the electron, kg
pub const ME: f64 = 9.109_56e-31;

/// Mass of the proton, kg
pub const MP: f64 = 1.67e-27;

/// Stefan-Boltzmann constant, MKS
pub const SIGMA: f64 = 5.669_56e-8;

/// Thomson cross-section, MKS
pub const SIGMAT: f64 = 6.65e-29;

/// Astronomical unit, metres
pub const AU: f64 = 1.495_978_706_91e11;

/// Solar luminosity, Watts
pub const LSUN: f64 = 3.826e26;

/// Solar mass, kg
pub const MSUN: f64 = 1.989e30;

/// Gravitational parameter of the Sun, SI (m^3 s^-2)
pub const GMSUN: f64 = 1.327_124_420_99e20;

/// Gravitational parameter of the Sun, AU^3 YR^-2
pub const GMSUNA: f64 = 39.476_927_033_270_655;

/// Gauss' gravitational constant sqrt(G*MSUN), AU^(3/2) day^-1
pub const KGAUSS: f64 = 0.017_202_098_95;

/// G*MSUN, AU^3 day^-2 (Gauss^2)
pub const GMGAUSS: f64 = KGAUSS * KGAUSS;

/// Absolute visual magnitude of the Sun
pub const MVSUN: f64 = 4.75;

/// Parsec, metres
pub const PC: f64 = 3.085_678e16;

/// Solar radius, metres
pub const RSUN: f64 = 6.9599e8;

/// Effective temperature of the Sun, Kelvin
pub const TSUN: f64 = 5700.0;

/// Number of seconds in a day
pub const DAY: f64 = 86400.0;

/// Length of Julian year in seconds
pub const YEAR: f64 = 365.25 * DAY;

/// Integer number of seconds in a day
pub const IDAY: i64 = 86400;

/// Number of seconds in an hour
pub const HOUR: f64 = 3600.0;

/// Number of seconds in a minute
pub const MINUTE: f64 = 60.0;

/// Pi
pub const PI: f64 = std::f64::consts::PI;

/// 2*Pi
pub const TWOPI: f64 = 2.0 * PI;

/// Wavelength of Halpha, Angstroms
pub const HALPHA: f64 = 6562.76;

/// Wavelength of Hbeta, Angstroms
pub const HBETA: f64 = 4861.327;

/// Wavelength of Hgamma, Angstroms
pub const HGAMMA: f64 = 4340.465;

/// Wavelength of Hdelta, Angstroms
pub const HDELTA: f64 = 4340.465;
