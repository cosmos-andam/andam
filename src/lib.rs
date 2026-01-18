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

pub mod constants;
pub mod units;
pub mod dynamics;
pub mod observations;
pub mod cmb;
pub mod structure;
pub mod perturbations;
pub mod advanced;
pub mod early_universe;
pub mod visualization;
pub mod statistics;
pub mod beyond_lcdm;

#[cfg(feature = "hdf5-storage")]
pub mod storage;

// Re-export commonly used items
pub use constants::*;
pub use units::{Length, Mass, Time, Energy};
pub use dynamics::Universe;
pub use observations::{
    luminosity_distance,
    angular_diameter_distance,
    comoving_distance,
};
pub use cmb::recombination::{ionization_fraction, recombination_redshift};
pub use cmb::fluctuations::{angular_power_spectrum, dimensionless_power_spectrum};
pub use structure::power_spectrum::matter_power_spectrum;
pub use perturbations::growth::{growth_factor, growth_rate};
pub use early_universe::nucleosynthesis::{run_bbn, helium_abundance, BBNParameters};
