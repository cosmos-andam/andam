//! Observational cosmology tools

pub mod distances;

pub use distances::{
    angular_diameter_distance, comoving_distance, distance_modulus, luminosity_distance,
};
