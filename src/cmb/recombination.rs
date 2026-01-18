//! CMB recombination physics
//!
//! Implements the Saha equation and ionization fraction calculation

use crate::constants::*;
use std::f64::consts::PI;

/// Saha equation for ionization equilibrium
///
/// n_e n_p / n_H = (m_e k T / 2π ℏ²)^(3/2) exp(-B/kT)
///
/// Returns ionization fraction X_e
pub fn saha_equation(temp_k: f64, n_h: f64) -> f64 {
    // Hydrogen binding energy [J]
    let binding_energy = 13.6 * EV_TO_J;

    // Prefactor
    let prefactor = ((M_E * K_B * temp_k) / (2.0 * PI * HBAR * HBAR)).powf(1.5);

    // Boltzmann factor
    let boltzmann = (-binding_energy / (K_B * temp_k)).exp();

    // Combined factor K = prefactor * boltzmann / n_h
    let k = prefactor * boltzmann / n_h;

    // Check for limiting cases to avoid numerical issues
    if k > 1e6 {
        return 1.0; // Fully ionized
    }
    if k < 1e-20 {
        return 0.0; // Fully neutral
    }

    // Solve quadratic equation for X_e
    // X_e² / (1 - X_e) = K
    // X_e² + K * X_e - K = 0
    let discriminant = k * k + 4.0 * k;
    let x_e = (-k + discriminant.sqrt()) / 2.0;

    x_e.clamp(0.0, 1.0)
}

/// Ionization fraction as a function of redshift
///
/// Uses Saha equation with corrections
pub fn ionization_fraction(z: f64, universe: &crate::dynamics::Universe) -> f64 {
    let temp = T_CMB * (1.0 + z);

    // Baryon number density today [m⁻³]
    let omega_b = 0.049; // Approximate baryon density parameter
    let rho_b = omega_b * universe.critical_density();
    let n_b = rho_b / M_P;

    // Scale with redshift
    let n_h = n_b * (1.0 + z).powi(3);

    saha_equation(temp, n_h)
}

/// Recombination redshift (when X_e ~ 0.5)
pub fn recombination_redshift(universe: &crate::dynamics::Universe) -> f64 {
    // Iterate to find z where X_e = 0.5
    let mut z_low = 500.0;
    let mut z_high = 2000.0;

    for _ in 0..50 {
        let z_mid = (z_low + z_high) / 2.0;
        let x_e = ionization_fraction(z_mid, universe);

        if x_e > 0.5 {
            // Too ionized, go to lower z (later time)
            z_high = z_mid;
        } else {
            // Too neutral, go to higher z (earlier time)
            z_low = z_mid;
        }
    }

    (z_low + z_high) / 2.0
}

/// Optical depth to Thomson scattering
///
/// τ = ∫ n_e σ_T c dt
pub fn optical_depth(z: f64, universe: &crate::dynamics::Universe) -> f64 {
    let n_steps = 1000;
    let dz = z / n_steps as f64;

    let mut tau = 0.0;

    for i in 0..n_steps {
        let zi = i as f64 * dz;
        let x_e = ionization_fraction(zi, universe);

        // Baryon number density
        let omega_b = 0.049;
        let rho_b = omega_b * universe.critical_density();
        let n_b = rho_b / M_P;
        let n_e = x_e * n_b * (1.0 + zi).powi(3);

        // dt/dz
        let h = universe.hubble_z(zi) * 1e3 / (PARSEC * 1e6); // SI units
        let dt_dz = 1.0 / ((1.0 + zi) * h);

        tau += n_e * SIGMA_T * C * dt_dz * dz;
    }

    tau
}

/// Visibility function g(z) = -dτ/dz exp(-τ)
///
/// Peaks at last scattering surface
pub fn visibility_function(z: f64, universe: &crate::dynamics::Universe) -> f64 {
    let dz = 1.0;
    let tau = optical_depth(z, universe);
    let tau_plus = optical_depth(z + dz, universe);
    let dtau_dz = (tau_plus - tau) / dz;

    -dtau_dz * (-tau).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_saha_high_temp() {
        let temp = 1e5; // High temperature
        let n_h = 1e6;
        let x_e = saha_equation(temp, n_h);

        // Should be nearly fully ionized
        assert!(x_e > 0.99);
    }

    #[test]
    fn test_saha_low_temp() {
        let temp = 1000.0; // Low temperature
        let n_h = 1e6;
        let x_e = saha_equation(temp, n_h);

        // Should be mostly neutral
        assert!(x_e < 0.1);
    }

    #[test]
    fn test_recombination_redshift() {
        let universe = crate::dynamics::Universe::benchmark();
        let z_rec = recombination_redshift(&universe);

        // Should be around z ~ 1100
        assert!(z_rec > 900.0 && z_rec < 1400.0);
    }
}
