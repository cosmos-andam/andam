//! Statistical analysis tools

pub mod fisher;
pub mod mcmc;

pub use fisher::FisherMatrix;
pub use mcmc::{Chain, MCMCSampler, Parameter};
