//! Observational cosmology tools

pub mod distances;

pub use distances::{
    luminosity_distance,
    angular_diameter_distance,
    comoving_distance,
    distance_modulus,
};
