pub mod model;
pub mod types;
pub mod grid;
pub mod brightness;
pub mod flux;
pub mod orchestration;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum LcurveError {
    #[error("{0}")]
    Generic(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("roche error: {0}")]
    Roche(#[from] lcurve_roche::RocheError),
    #[error("subs error: {0}")]
    Subs(#[from] lcurve_subs::SubsError),
}
