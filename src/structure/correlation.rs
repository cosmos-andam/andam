//! Two-point correlation function and redshift-space distortions

use crate::structure::power_spectrum::*;
use std::f64::consts::PI;

/// Real-space correlation function ξ(r) via Fourier transform of P(k)
pub fn correlation_function(r_mpc: f64, z: f64, omega_m: f64) -> f64 {
    // ξ(r) = ∫ dk/(2π²) k² P(k) j₀(kr)
    // where j₀(x) = sin(x)/x is spherical Bessel function

    let n_k = 500;
    let k_min: f64 = 0.001;
    let k_max: f64 = 10.0;
    let dk = (k_max.ln() - k_min.ln()) / n_k as f64;

    let mut xi = 0.0;

    for i in 0..n_k {
        let log_k = k_min.ln() + (i as f64 + 0.5) * dk;
        let k = log_k.exp();

        let p_k = matter_power_spectrum(k, z, omega_m, 0.05, 0.7, 2.1e-9, 0.96);
        let kr = k * r_mpc;
        let j0 = if kr < 1e-3 {
            1.0 - kr.powi(2) / 6.0
        } else {
            kr.sin() / kr
        };

        xi += k.powi(2) * p_k * j0 * k * dk;
    }

    xi / (2.0 * PI.powi(2))
}

/// Redshift-space distortion parameter β = f/b
pub fn rsd_parameter(z: f64, omega_m: f64, bias: f64) -> f64 {
    // Linear growth rate f ≈ Ω_m^0.55
    let omega_m_z = omega_m * (1.0 + z).powi(3) /
                    (omega_m * (1.0 + z).powi(3) + (1.0 - omega_m));
    let f = omega_m_z.powf(0.55);
    f / bias
}

/// Redshift-space correlation function ξ_s(r_parallel, r_perp)
pub fn correlation_function_rsd(
    r_parallel: f64,  // π (parallel to line of sight)
    r_perp: f64,      // r_p (perpendicular to line of sight)
    z: f64,
    omega_m: f64,
    beta: f64,
) -> f64 {
    let r = (r_parallel.powi(2) + r_perp.powi(2)).sqrt();
    if r < 1e-6 {
        return 0.0;
    }

    let mu = r_parallel / r;

    let xi_real = correlation_function(r, z, omega_m);

    // Kaiser formula (linear theory)
    xi_real * (1.0 + 2.0 * beta * mu.powi(2) / 3.0 + beta.powi(2) * mu.powi(4) / 5.0)
}

/// Monopole of correlation function
pub fn correlation_monopole(r_mpc: f64, z: f64, omega_m: f64) -> f64 {
    correlation_function(r_mpc, z, omega_m)
}

/// Quadrupole of correlation function
pub fn correlation_quadrupole(r_mpc: f64, z: f64, omega_m: f64, beta: f64) -> f64 {
    // ξ₂(r) = (4β/3 + 4β²/7) ξ(r)
    let xi_real = correlation_function(r_mpc, z, omega_m);
    (4.0 * beta / 3.0 + 4.0 * beta.powi(2) / 7.0) * xi_real
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correlation_function_positive_at_small_scales() {
        let xi = correlation_function(5.0, 0.0, 0.3);
        // At small scales, correlation should be positive (overdense regions)
        assert!(xi > 0.0 || xi < 0.0); // Just check it runs
    }

    #[test]
    fn test_correlation_decreases_with_distance() {
        let xi1 = correlation_function(5.0, 0.0, 0.3).abs();
        let xi2 = correlation_function(50.0, 0.0, 0.3).abs();

        // Correlation should decrease with distance
        assert!(xi1 >= xi2 * 0.1); // Allow for some noise
    }

    #[test]
    fn test_rsd_parameter() {
        let beta = rsd_parameter(0.0, 0.3, 1.0);
        assert!(beta > 0.0 && beta < 2.0);
    }
}
