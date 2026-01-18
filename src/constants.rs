//! Physical and cosmological constants
//!
//! All constants are in SI units unless otherwise specified

use std::f64::consts::PI;

/// Speed of light in vacuum [m/s]
pub const C: f64 = 2.99792458e8;

/// Gravitational constant [m³/(kg·s²)]
pub const G: f64 = 6.67430e-11;

/// Planck constant [J·s]
pub const H_PLANCK: f64 = 6.62607015e-34;

/// Reduced Planck constant (ℏ) [J·s]
pub const HBAR: f64 = 1.054571817e-34;

/// Boltzmann constant [J/K]
pub const K_B: f64 = 1.380649e-23;

/// Stefan-Boltzmann constant [W/(m²·K⁴)]
pub const SIGMA_SB: f64 = 5.670374419e-8;

/// Electron mass \[kg\]
pub const M_E: f64 = 9.1093837015e-31;

/// Proton mass \[kg\]
pub const M_P: f64 = 1.67262192369e-27;

/// Neutron mass \[kg\]
pub const M_N: f64 = 1.67492749804e-27;

/// Atomic mass unit \[kg\]
pub const AMU: f64 = 1.66053906660e-27;

/// Elementary charge \[C\]
pub const E_CHARGE: f64 = 1.602176634e-19;

/// Thomson scattering cross-section \[m²\]
pub const SIGMA_T: f64 = 6.6524587321e-29;

/// Fine structure constant (dimensionless)
pub const ALPHA: f64 = 7.2973525693e-3;

// Astronomical constants

/// Astronomical unit \[m\]
pub const AU: f64 = 1.495978707e11;

/// Parsec \[m\]
pub const PARSEC: f64 = 3.0856775814913673e16;

/// Solar mass \[kg\]
pub const M_SUN: f64 = 1.98847e30;

/// Solar luminosity \[W\]
pub const L_SUN: f64 = 3.828e26;

/// Solar radius \[m\]
pub const R_SUN: f64 = 6.957e8;

// Cosmological constants

/// Hubble constant (Planck 2018) [km/s/Mpc]
pub const H0_PLANCK: f64 = 67.66;

/// Critical density coefficient [kg/(m³·(km/s/Mpc)²)]
pub const RHO_CRIT_COEFF: f64 = 3.0 / (8.0 * PI * G);

/// CMB temperature today \[K\]
pub const T_CMB: f64 = 2.7255;

/// CMB photon number density today [m⁻³]
pub const N_GAMMA: f64 = 4.105e8;

/// Neutrino temperature today \[K\]
pub const T_NU: f64 = 1.95;

// Energy conversions

/// Electron volt to Joules
pub const EV_TO_J: f64 = 1.602176634e-19;

/// Joules to electron volt
pub const J_TO_EV: f64 = 6.241509074e18;

// Time conversions

/// Year in seconds
pub const YEAR: f64 = 3.15576e7;

/// Gigayear in seconds
pub const GYR: f64 = 3.15576e16;

/// Hubble time \[s\]
pub fn hubble_time(h0_km_s_mpc: f64) -> f64 {
    PARSEC * 1e6 / (h0_km_s_mpc * 1e3)
}

/// Critical density today [kg/m³]
pub fn critical_density(h0_km_s_mpc: f64) -> f64 {
    let h_si = h0_km_s_mpc * 1e3 / (PARSEC * 1e6); // Convert to SI
    3.0 * h_si * h_si / (8.0 * PI * G)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_speed_of_light() {
        assert_eq!(C, 2.99792458e8);
    }

    #[test]
    fn test_hubble_time() {
        let t_h = hubble_time(70.0);
        let expected = 4.4e17; // approximately 14 Gyr in seconds
        assert_relative_eq!(t_h, expected, epsilon = 1e16);
    }

    #[test]
    fn test_critical_density() {
        let rho_c = critical_density(70.0);
        assert!(rho_c > 0.0);
        assert!(rho_c < 1e-25); // Should be very small
    }
}
