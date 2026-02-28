use lcurve_subs::Vec3;

/// Eclipse information: vector of (ingress, egress) phase pairs.
pub type EclipseVec = Vec<(f64, f64)>;

/// A single surface element on a star, disc, or bright spot.
#[derive(Clone)]
pub struct Point {
    /// Position vector (units of binary separation)
    pub posn: Vec3,
    /// Outward-facing direction (unit vector)
    pub dirn: Vec3,
    /// Area of element (units of separation^2)
    pub area: f32,
    /// Gravity of element (relative to reference)
    pub gravity: f32,
    /// Ingress/egress phase pairs for eclipses
    pub eclipse: EclipseVec,
    /// Brightness * area
    pub flux: f32,
}

impl Default for Point {
    fn default() -> Self {
        Point {
            posn: Vec3::ZERO,
            dirn: Vec3::ZERO,
            area: 0.0,
            gravity: 1.0,
            eclipse: Vec::new(),
            flux: 0.0,
        }
    }
}

impl Point {
    pub fn new(posn: Vec3, dirn: Vec3, area: f64, gravity: f64, eclipse: EclipseVec) -> Self {
        Point {
            posn,
            dirn,
            area: area as f32,
            gravity: gravity as f32,
            eclipse,
            flux: 0.0,
        }
    }

    /// Returns true if the point is visible (not eclipsed) at the given phase.
    #[inline]
    pub fn visible(&self, phase: f64) -> bool {
        let phi = phase - phase.floor();
        for &(ing, eg) in &self.eclipse {
            if (phi >= ing && phi <= eg) || phi <= eg - 1.0 {
                return false;
            }
        }
        true
    }
}

/// Limb darkening type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LDCType {
    /// Polynomial: I(mu) = 1 - sum_i a_i (1-mu)^i
    Poly,
    /// Claret 4-coefficient: I(mu) = 1 - sum_i a_i (1 - mu^(i/2))
    Claret,
}

/// Limb darkening coefficients.
#[derive(Debug, Clone, Copy)]
pub struct LDC {
    ldc1: f64,
    ldc2: f64,
    ldc3: f64,
    ldc4: f64,
    mucrit: f64,
    ltype: LDCType,
}

impl Default for LDC {
    fn default() -> Self {
        LDC {
            ldc1: 0.0,
            ldc2: 0.0,
            ldc3: 0.0,
            ldc4: 0.0,
            mucrit: 0.0,
            ltype: LDCType::Poly,
        }
    }
}

impl LDC {
    pub fn new(ldc1: f64, ldc2: f64, ldc3: f64, ldc4: f64, mucrit: f64, ltype: LDCType) -> Self {
        LDC { ldc1, ldc2, ldc3, ldc4, mucrit, ltype }
    }

    /// Compute I(mu) — the limb-darkened intensity.
    #[inline]
    pub fn imu(&self, mu: f64) -> f64 {
        if mu <= 0.0 {
            return 0.0;
        }
        let mu = mu.min(1.0);
        let ommu = 1.0 - mu;
        let mut im = 1.0;
        match self.ltype {
            LDCType::Poly => {
                im -= ommu * (self.ldc1 + ommu * (self.ldc2 + ommu * (self.ldc3 + ommu * self.ldc4)));
            }
            LDCType::Claret => {
                im -= self.ldc1 + self.ldc2 + self.ldc3 + self.ldc4;
                let msq = mu.sqrt();
                im += msq * (self.ldc1 + msq * (self.ldc2 + msq * (self.ldc3 + msq * self.ldc4)));
            }
        }
        im
    }

    /// Returns true if mu exceeds the critical value.
    #[inline]
    pub fn see(&self, mu: f64) -> bool {
        mu > self.mucrit
    }
}

/// Grid interpolation parameters for switching between coarse and fine grids.
#[derive(Debug, Clone)]
pub struct Ginterp {
    /// Start phase of coarse grid (0 – 0.5)
    pub phase1: f64,
    /// End phase of coarse grid (0 – 0.5)
    pub phase2: f64,
    /// Scale factor star 1 at phase1
    pub scale11: f64,
    /// Scale factor star 1 at 1-phase1
    pub scale12: f64,
    /// Scale factor star 2 at -phase2
    pub scale21: f64,
    /// Scale factor star 2 at phase2
    pub scale22: f64,
}

impl Default for Ginterp {
    fn default() -> Self {
        Ginterp {
            phase1: 0.05,
            phase2: 0.45,
            scale11: 1.0,
            scale12: 1.0,
            scale21: 1.0,
            scale22: 1.0,
        }
    }
}

impl Ginterp {
    /// Returns scale factor for star 1 at a given phase.
    #[inline]
    pub fn scale1(&self, phase: f64) -> f64 {
        let pnorm = phase - phase.floor();
        if pnorm <= self.phase1 || pnorm >= 1.0 - self.phase1 {
            1.0
        } else {
            (self.scale11 * (1.0 - self.phase1 - pnorm) + self.scale12 * (pnorm - self.phase1))
                / (1.0 - 2.0 * self.phase1)
        }
    }

    /// Returns scale factor for star 2 at a given phase.
    #[inline]
    pub fn scale2(&self, phase: f64) -> f64 {
        let pnorm = phase - phase.floor();
        if pnorm >= self.phase2 && pnorm <= 1.0 - self.phase2 {
            1.0
        } else if pnorm < 0.5 {
            (self.scale22 * (self.phase2 - pnorm) + self.scale21 * (pnorm + self.phase2))
                / (2.0 * self.phase2)
        } else {
            (self.scale21 * (1.0 + self.phase2 - pnorm) + self.scale22 * (pnorm - 1.0 + self.phase2))
                / (2.0 * self.phase2)
        }
    }

    /// Returns integer type to represent grid situation.
    /// 1 = fine star1, coarse star2
    /// 2 = coarse both
    /// 3 = coarse star1, fine star2
    #[inline]
    pub fn grid_type(&self, phase: f64) -> i32 {
        let pnorm = phase - phase.floor();
        if pnorm <= self.phase1 || pnorm >= 1.0 - self.phase1 {
            1
        } else if (pnorm > self.phase1 && pnorm < self.phase2)
            || (pnorm > 1.0 - self.phase2 && pnorm < 1.0 - self.phase1)
        {
            2
        } else {
            3
        }
    }
}

/// A single data point of a light curve.
#[derive(Debug, Clone)]
pub struct Datum {
    /// Mid-exposure time
    pub time: f64,
    /// Exposure length (same units as time)
    pub expose: f64,
    /// Flux
    pub flux: f64,
    /// Uncertainty on flux
    pub ferr: f64,
    /// Weight factor
    pub weight: f64,
    /// Number of sub-divisions for exposure smearing
    pub ndiv: i32,
}

/// A collection of data points.
pub type Data = Vec<Datum>;

/// Read data from an ASCII file.
///
/// Lines starting with # or blank/whitespace are skipped.
/// Each data line has: time expose flux ferr weight ndiv
pub fn read_data(path: &str) -> Result<Data, crate::LcurveError> {
    let contents = std::fs::read_to_string(path)?;
    let mut data = Vec::new();
    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 6 {
            continue;
        }
        let time: f64 = parts[0].parse().map_err(|_| crate::LcurveError::Parse(format!("bad time: {}", parts[0])))?;
        let expose: f64 = parts[1].parse().map_err(|_| crate::LcurveError::Parse(format!("bad expose: {}", parts[1])))?;
        let flux: f64 = parts[2].parse().map_err(|_| crate::LcurveError::Parse(format!("bad flux: {}", parts[2])))?;
        let mut ferr: f64 = parts[3].parse().map_err(|_| crate::LcurveError::Parse(format!("bad ferr: {}", parts[3])))?;
        let mut weight: f64 = parts[4].parse().map_err(|_| crate::LcurveError::Parse(format!("bad weight: {}", parts[4])))?;
        let ndiv: i32 = parts[5].parse().map_err(|_| crate::LcurveError::Parse(format!("bad ndiv: {}", parts[5])))?;
        if ferr < 0.0 {
            ferr = -ferr;
            weight = 0.0;
        }
        data.push(Datum { time, expose, flux, ferr, weight, ndiv });
    }
    Ok(data)
}

/// Write data to an ASCII file.
pub fn write_data(path: &str, data: &[Datum]) -> Result<(), crate::LcurveError> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    for d in data {
        writeln!(f, "{:17.10e} {} {:10.6e} {} {} {}",
                 d.time, d.expose, d.flux, d.ferr, d.weight, d.ndiv)?;
    }
    Ok(())
}
