//! Cosmic Microwave Background physics

pub mod fluctuations;
pub mod polarization;
pub mod recombination;

pub use fluctuations::{
    acoustic_peak_positions, angular_power_spectrum, dimensionless_power_spectrum,
};
pub use polarization::{decompose_eb, PolarizationSpectrum, StokesParameters};
pub use recombination::{
    ionization_fraction, optical_depth, recombination_redshift, saha_equation,
};
