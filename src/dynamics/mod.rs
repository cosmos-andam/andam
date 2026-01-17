//! Cosmic dynamics and evolution
//!
//! Implements Friedmann equations and universe evolution

pub mod friedmann;
pub mod components;

pub use friedmann::Universe;
pub use components::{Component, ComponentType};
