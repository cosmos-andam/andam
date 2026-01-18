//! Cosmic dynamics and evolution
//!
//! Implements Friedmann equations and universe evolution

pub mod components;
pub mod friedmann;

pub use components::{Component, ComponentType};
pub use friedmann::Universe;
