//! # Andam
//!
//! A comprehensive Rust library for cosmological calculations and visualizations.
//!
//! ## Features
//!
//! - Friedmann equation solver
//! - Multiple universe components (matter, radiation, dark energy)
//! - Cosmological distance measures
//! - CMB recombination physics and angular power spectrum
//! - Matter power spectrum
//! - Growth factor for structure formation
//! - Weak gravitational lensing
//! - Unit conversions
//! - 2D and interactive plotting capabilities
//!
//! ## Example
//!
//! ```rust
//! use andam::dynamics::Universe;
//! use andam::observations::luminosity_distance;
//!
//! let universe = Universe::benchmark();
//! let age = universe.age_today();
//! println!("Age of universe: {:.2} Gyr", age);
//!
//! let z = 1.0;
//! let d_l = luminosity_distance(z, &universe);
//! println!("Luminosity distance to z={}: {:.1} Mpc", z, d_l);
//! ```

pub mod advanced;
pub mod beyond_lcdm;
pub mod cmb;
pub mod constants;
pub mod dynamics;
pub mod early_universe;
pub mod observations;
pub mod perturbations;
pub mod statistics;
pub mod structure;
pub mod units;

#[cfg(feature = "plotting")]
pub mod visualization;

#[cfg(feature = "hdf5-storage")]
pub mod storage;

// Re-export commonly used items
pub use cmb::fluctuations::{angular_power_spectrum, dimensionless_power_spectrum};
pub use cmb::recombination::{ionization_fraction, recombination_redshift};
pub use constants::*;
pub use dynamics::Universe;
pub use early_universe::nucleosynthesis::{helium_abundance, run_bbn, BBNParameters};
pub use observations::{angular_diameter_distance, comoving_distance, luminosity_distance};
pub use perturbations::growth::{growth_factor, growth_rate};
pub use structure::power_spectrum::matter_power_spectrum;
pub use units::{Energy, Length, Mass, Time};
