//! Statistical analysis tools

pub mod mcmc;
pub mod fisher;

pub use mcmc::{MCMCSampler, Chain, Parameter};
pub use fisher::FisherMatrix;
