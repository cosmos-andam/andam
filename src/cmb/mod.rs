//! Cosmic Microwave Background physics

pub mod recombination;
pub mod fluctuations;
pub mod polarization;

pub use recombination::{
    saha_equation,
    ionization_fraction,
    recombination_redshift,
    optical_depth,
};
pub use fluctuations::{
    angular_power_spectrum,
    dimensionless_power_spectrum,
    acoustic_peak_positions,
};
pub use polarization::{
    StokesParameters,
    PolarizationSpectrum,
    decompose_eb,
};
