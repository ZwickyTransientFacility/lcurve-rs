pub mod vec3;
pub mod constants;
pub mod planck;
pub mod numerical;

pub use vec3::Vec3;
pub use constants::*;
pub use planck::{planck, dplanck, dlpdlt};
pub use numerical::rtsafe::rtsafe;
pub use numerical::bsstep::{bsstep, BsState};
pub use numerical::dbrent::dbrent;
pub use numerical::svdfit::svdfit;

/// Inline utility: square of a value
#[inline]
pub fn sqr(x: f64) -> f64 {
    x * x
}

/// Transfer sign of b to magnitude of a
#[inline]
pub fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

/// sqrt(a^2 + b^2) without overflow
#[inline]
pub fn pythag(a: f64, b: f64) -> f64 {
    a.hypot(b)
}

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SubsError {
    #[error("{0}")]
    Generic(String),
    #[error("root not bracketed: x1={x1}, x2={x2}, f1={f1}, f2={f2}")]
    RootNotBracketed { x1: f64, x2: f64, f1: f64, f2: f64 },
    #[error("too many iterations in {0}")]
    TooManyIterations(String),
    #[error("null vector in unit()")]
    NullVector,
    #[error("SVD error: {0}")]
    Svd(String),
}
