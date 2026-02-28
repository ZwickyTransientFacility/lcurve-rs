pub mod potential;
pub mod lagrange;
pub mod lobes;
pub mod eclipse;
pub mod disc_eclipse;
pub mod stream;
pub mod surface;
pub mod misc;

pub use lcurve_subs::Vec3;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum RocheError {
    #[error("{0}")]
    Generic(String),
    #[error("subs error: {0}")]
    Subs(#[from] lcurve_subs::SubsError),
}

/// Which star in the binary system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Star {
    Primary,
    Secondary,
}
