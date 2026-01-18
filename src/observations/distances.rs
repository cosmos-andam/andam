//! Cosmological distance measures
//!
//! Implements various distance definitions used in cosmology:
//! - Comoving distance
//! - Luminosity distance
//! - Angular diameter distance
//! - Distance modulus

use crate::constants::*;
use crate::dynamics::Universe;

/// Integrand for comoving distance: 1/H(z)
fn comoving_integrand(z: f64, universe: &Universe) -> f64 {
    1.0 / universe.hubble_z(z)
}

/// Comoving distance to redshift z \[Mpc\]
///
/// d_c = (c/H_0) ∫₀^z dz'/H(z')
pub fn comoving_distance(z: f64, universe: &Universe) -> f64 {
    let c_over_h0 = C / (universe.h0 * 1e3) * PARSEC * 1e-6; // Mpc

    // Simpson's rule integration
    let n = 1000;
    let dz = z / n as f64;
    let mut sum = comoving_integrand(0.0, universe) + comoving_integrand(z, universe);

    for i in 1..n {
        let zi = i as f64 * dz;
        let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += weight * comoving_integrand(zi, universe);
    }

    c_over_h0 * sum * dz / 3.0
}

/// Transverse comoving distance \[Mpc\]
///
/// Accounts for curvature:
/// - Flat: d_M = d_c
/// - Open: d_M = (c/H_0)/√|Ω_k| sinh(√|Ω_k| d_c H_0/c)
/// - Closed: d_M = (c/H_0)/√Ω_k sin(√Ω_k d_c H_0/c)
pub fn transverse_comoving_distance(z: f64, universe: &Universe) -> f64 {
    let d_c = comoving_distance(z, universe);
    let omega_k = 1.0 - universe.omega_total();
    let c_over_h0 = C / (universe.h0 * 1e3) * PARSEC * 1e-6;

    if omega_k.abs() < 1e-5 {
        // Flat universe
        d_c
    } else if omega_k > 0.0 {
        // Open universe
        let sqrt_ok = omega_k.sqrt();
        (c_over_h0 / sqrt_ok) * (sqrt_ok * d_c / c_over_h0).sinh()
    } else {
        // Closed universe
        let sqrt_ok = (-omega_k).sqrt();
        (c_over_h0 / sqrt_ok) * (sqrt_ok * d_c / c_over_h0).sin()
    }
}

/// Luminosity distance \[Mpc\]
///
/// d_L = (1+z) d_M
pub fn luminosity_distance(z: f64, universe: &Universe) -> f64 {
    (1.0 + z) * transverse_comoving_distance(z, universe)
}

/// Angular diameter distance \[Mpc\]
///
/// d_A = d_M / (1+z)
pub fn angular_diameter_distance(z: f64, universe: &Universe) -> f64 {
    transverse_comoving_distance(z, universe) / (1.0 + z)
}

/// Distance modulus (m - M)
///
/// μ = 5 log₁₀(d_L / 10 pc)
pub fn distance_modulus(z: f64, universe: &Universe) -> f64 {
    let d_l = luminosity_distance(z, universe);
    let d_l_pc = d_l * 1e6; // Convert Mpc to pc
    5.0 * (d_l_pc / 10.0).log10()
}

/// Comoving volume element dV_c/dz/dΩ [Mpc³/sr]
pub fn comoving_volume_element(z: f64, universe: &Universe) -> f64 {
    let d_h = C / (universe.h0 * 1e3) * PARSEC * 1e-6; // Hubble distance \[Mpc\]
    let d_m = transverse_comoving_distance(z, universe);
    let e_z = universe.hubble_z(z) / universe.h0;

    d_h * d_m * d_m / e_z
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_comoving_distance_z0() {
        let universe = Universe::benchmark();
        let d = comoving_distance(0.0, &universe);
        assert_relative_eq!(d, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_distance_relations() {
        let universe = Universe::benchmark();
        let z = 1.0;

        let d_l = luminosity_distance(z, &universe);
        let d_a = angular_diameter_distance(z, &universe);

        // d_L = (1+z)² d_A
        assert_relative_eq!(d_l, (1.0 + z).powi(2) * d_a, epsilon = 1e-6);
    }

    #[test]
    fn test_flat_universe_distances() {
        let universe = Universe::benchmark();
        let z = 0.5;

        let d_c = comoving_distance(z, &universe);
        let d_m = transverse_comoving_distance(z, &universe);

        // Should be equal for flat universe (within numerical precision)
        assert_relative_eq!(d_c, d_m, max_relative = 1e-6);
    }
}
