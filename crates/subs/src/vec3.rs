use crate::SubsError;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Cartesian 3-vector in double precision.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    /// Zero vector.
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

    /// Create from components.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Set all components.
    #[inline]
    pub fn set(&mut self, x: f64, y: f64, z: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    /// Squared magnitude.
    #[inline]
    pub fn sqr(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Magnitude (length).
    #[inline]
    pub fn length(&self) -> f64 {
        self.sqr().sqrt()
    }

    /// Normalize to unit vector. Returns error if zero vector.
    pub fn unit(&mut self) -> Result<(), SubsError> {
        let norm = self.sqr();
        if norm > 0.0 {
            let norm = norm.sqrt();
            self.x /= norm;
            self.y /= norm;
            self.z /= norm;
            Ok(())
        } else {
            Err(SubsError::NullVector)
        }
    }

    /// Return a unit vector in this direction (non-mutating).
    pub fn unitized(&self) -> Result<Vec3, SubsError> {
        let mut v = *self;
        v.unit()?;
        Ok(v)
    }

    /// Unit vector along X axis.
    #[inline]
    pub fn unit_x() -> Self {
        Vec3 { x: 1.0, y: 0.0, z: 0.0 }
    }

    /// Unit vector along Y axis.
    #[inline]
    pub fn unit_y() -> Self {
        Vec3 { x: 0.0, y: 1.0, z: 0.0 }
    }

    /// Unit vector along Z axis.
    #[inline]
    pub fn unit_z() -> Self {
        Vec3 { x: 0.0, y: 0.0, z: 1.0 }
    }

    /// Dot product.
    #[inline]
    pub fn dot(a: &Vec3, b: &Vec3) -> f64 {
        a.x * b.x + a.y * b.y + a.z * b.z
    }

    /// Cross product.
    #[inline]
    pub fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
        Vec3 {
            x: a.y * b.z - a.z * b.y,
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x,
        }
    }
}

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {}", self.x, self.y, self.z)
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3 { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, rhs: Vec3) -> Vec3 {
        Vec3 { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3 { x: -self.x, y: -self.y, z: -self.z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, rhs: f64) -> Vec3 {
        Vec3 { x: self.x * rhs, y: self.y * rhs, z: self.z * rhs }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, rhs: Vec3) -> Vec3 {
        Vec3 { x: self * rhs.x, y: self * rhs.y, z: self * rhs.z }
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, rhs: f64) -> Vec3 {
        Vec3 { x: self.x / rhs, y: self.y / rhs, z: self.z / rhs }
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Vec3) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Vec3) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl MulAssign<f64> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl DivAssign<f64> for Vec3 {
    #[inline]
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert!((Vec3::dot(&a, &b) - 32.0).abs() < 1e-15);
    }

    #[test]
    fn test_cross() {
        let a = Vec3::unit_x();
        let b = Vec3::unit_y();
        let c = Vec3::cross(&a, &b);
        assert!((c.x).abs() < 1e-15);
        assert!((c.y).abs() < 1e-15);
        assert!((c.z - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_length() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!((v.length() - 5.0).abs() < 1e-15);
    }

    #[test]
    fn test_unit() {
        let mut v = Vec3::new(3.0, 4.0, 0.0);
        v.unit().unwrap();
        assert!((v.length() - 1.0).abs() < 1e-15);
        assert!((v.x - 0.6).abs() < 1e-15);
        assert!((v.y - 0.8).abs() < 1e-15);
    }

    #[test]
    fn test_ops() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);

        let c = a + b;
        assert!((c.x - 5.0).abs() < 1e-15);

        let d = a - b;
        assert!((d.x - (-3.0)).abs() < 1e-15);

        let e = 2.0 * a;
        assert!((e.x - 2.0).abs() < 1e-15);
        assert!((e.y - 4.0).abs() < 1e-15);

        let f = -a;
        assert!((f.x - (-1.0)).abs() < 1e-15);
    }
}
