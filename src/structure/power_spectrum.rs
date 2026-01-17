//! Matter power spectrum calculations

use crate::constants::T_CMB;
use std::f64::consts::PI;

/// Primordial power spectrum from inflation
///
/// P(k) = A_s (k/k_pivot)^(n_s - 1)
///
/// # Arguments
/// * `k` - Wavenumber [h/Mpc]
/// * `amplitude` - A_s (typically ~2.1e-9)
/// * `spectral_index` - n_s (typically ~0.96)
/// * `k_pivot` - Pivot scale [h/Mpc] (typically 0.05)
pub fn primordial_power_spectrum(
    k: f64,
    amplitude: f64,
    spectral_index: f64,
    k_pivot: f64,
) -> f64 {
    amplitude * (k / k_pivot).powf(spectral_index - 1.0)
}

/// Transfer function (simplified Eisenstein & Hu form)
///
/// T(k) suppresses power on small scales
pub fn transfer_function_eh(k: f64, omega_m: f64, omega_b: f64, h: f64) -> f64 {
    // Eisenstein & Hu (1998) fitting formula
    let omega_m_h2 = omega_m * h * h;
    let omega_b_h2 = omega_b * h * h;
    let theta_cmb = T_CMB / 2.7;

    // Sound horizon (not used in simplified version)
    let _s = 44.5 * (omega_m_h2 / theta_cmb.powi(4)).ln()
        / (1.0 + 10.0 * omega_b_h2.powf(0.75)).sqrt();

    // Silk damping scale (not used in simplified version)
    let _k_silk = 1.6 * omega_b_h2.powf(0.52) * omega_m_h2.powf(0.73)
        * (1.0 + (10.4 * omega_m_h2).powf(-0.95));

    let q = k / (13.41 * omega_m_h2.powf(0.5));
    let l = (2.0 * 2.718 + 1.8 * q).ln();
    let c = 14.2 + 731.0 / (1.0 + 62.5 * q);

    let t = l / (l + c * q * q);

    t
}

/// Linear matter power spectrum P(k)
///
/// P(k) = A_s (k/k_pivot)^(n_s-1) T²(k) D²(z) (2π²/k³)
pub fn matter_power_spectrum(
    k: f64,
    z: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
    amplitude: f64,
    spectral_index: f64,
) -> f64 {
    let p_prim = primordial_power_spectrum(k, amplitude, spectral_index, 0.05);
    let t = transfer_function_eh(k, omega_m, omega_b, h);

    // Growth factor (simplified for matter-dominated)
    let d = 1.0 / (1.0 + z);

    // Dimensional factor
    let factor = 2.0 * PI * PI / (k * k * k);

    p_prim * t * t * d * d * factor
}

/// Dimensionless power spectrum Δ²(k)
pub fn dimensionless_power(
    k: f64,
    z: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
    amplitude: f64,
    spectral_index: f64,
) -> f64 {
    let p_k = matter_power_spectrum(k, z, omega_m, omega_b, h, amplitude, spectral_index);
    k * k * k * p_k / (2.0 * PI * PI)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primordial_power_at_pivot() {
        let p = primordial_power_spectrum(0.05, 2.1e-9, 0.96, 0.05);
        assert!((p - 2.1e-9).abs() < 1e-15);
    }

    #[test]
    fn test_transfer_function_limiting_cases() {
        // Should be ~1 at large scales (small k)
        let t_large = transfer_function_eh(1e-4, 0.3, 0.05, 0.7);
        assert!(t_large > 0.9 && t_large <= 1.0);

        // Should be <1 at small scales (large k)
        let t_small = transfer_function_eh(10.0, 0.3, 0.05, 0.7);
        assert!(t_small < 0.5);
    }
}
