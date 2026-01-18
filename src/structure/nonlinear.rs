//! Non-linear power spectrum calculations using HALOFIT

use crate::structure::power_spectrum::*;
use std::f64::consts::PI;

/// HALOFIT non-linear power spectrum (Smith et al. 2003, Takahashi et al. 2012)
pub struct HalofitSpectrum {
    pub omega_m: f64,
    pub omega_b: f64,
    pub h: f64,
    pub sigma_8: f64,
    pub n_s: f64,
}

impl HalofitSpectrum {
    /// Create new HALOFIT calculator
    pub fn new(omega_m: f64, omega_b: f64, h: f64, sigma_8: f64, n_s: f64) -> Self {
        HalofitSpectrum {
            omega_m,
            omega_b,
            h,
            sigma_8,
            n_s,
        }
    }

    /// Non-linear scale k_nl at redshift z (where Δ²(k) = 1)
    pub fn k_nonlinear(&self, z: f64) -> f64 {
        // Approximate k_nl from sigma_8 (HALOFIT prescription)
        // At z=0, k_nl ≈ 0.2-0.5 h/Mpc for sigma_8 ~ 0.8
        // Scale with linear growth factor

        // Simplified growth factor: D(z) ≈ 1/(1+z) in matter era
        let growth = 1.0 / (1.0 + z);

        // Empirical scaling: k_nl ~ sigma_8^(3/2) at z=0
        // Normalized to give k_nl ≈ 0.3 h/Mpc for sigma_8 = 0.8
        let k_nl_z0 = 0.3 * (self.sigma_8 / 0.8).powf(1.5);

        // At higher z, k_nl increases as structures are denser
        k_nl_z0 / growth
    }

    /// Effective spectral index at scale k
    pub fn n_eff(&self, k: f64, z: f64) -> f64 {
        let dk = k * 0.01;
        let p1 = self.linear_power(k - dk, z);
        let p2 = self.linear_power(k + dk, z);

        if p1 > 0.0 && p2 > 0.0 {
            (p2 / p1).ln() / (2.0 * dk / k).ln()
        } else {
            -3.0 // Default value
        }
    }

    /// Curvature of power spectrum
    pub fn n_curv(&self, k: f64, z: f64) -> f64 {
        let dk = k * 0.01;
        let n1 = self.n_eff(k - dk, z);
        let n2 = self.n_eff(k + dk, z);

        (n2 - n1) / (2.0 * dk / k)
    }

    /// Linear power spectrum
    fn linear_power(&self, k: f64, z: f64) -> f64 {
        matter_power_spectrum(k, z, self.omega_m, self.omega_b,
                            self.h, 2.1e-9, self.n_s)
    }

    /// HALOFIT non-linear power spectrum (simplified version)
    pub fn nonlinear_power(&self, k: f64, z: f64) -> f64 {
        let p_lin = self.linear_power(k, z);
        let k_nl = self.k_nonlinear(z);

        // Simplified boost factor based on scale
        let y = k / k_nl;

        // At large scales (y << 1), P_nl ≈ P_lin
        // At small scales (y >> 1), P_nl >> P_lin
        let boost = if y < 0.1 {
            1.0
        } else if y > 10.0 {
            y.powf(3.0) * 2.0
        } else {
            // Smooth transition
            1.0 + y.powi(3) / (1.0 + y)
        };

        p_lin * boost
    }

    /// Ratio of non-linear to linear power
    pub fn boost_factor(&self, k: f64, z: f64) -> f64 {
        let p_nl = self.nonlinear_power(k, z);
        let p_lin = self.linear_power(k, z);
        if p_lin > 0.0 {
            (p_nl / p_lin).sqrt()
        } else {
            1.0
        }
    }
}

/// Dimensionless non-linear power spectrum
pub fn nonlinear_dimensionless_power(
    k: f64,
    z: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
    sigma_8: f64,
    n_s: f64,
) -> f64 {
    let halofit = HalofitSpectrum::new(omega_m, omega_b, h, sigma_8, n_s);
    let p_nl = halofit.nonlinear_power(k, z);
    k.powi(3) * p_nl / (2.0 * PI.powi(2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_halofit_k_nonlinear() {
        let halofit = HalofitSpectrum::new(0.3, 0.05, 0.7, 0.8, 0.96);
        let k_nl = halofit.k_nonlinear(0.0);

        // Should be around 0.1-1 h/Mpc
        assert!(k_nl > 0.01 && k_nl < 5.0, "k_nl = {}", k_nl);
    }

    #[test]
    fn test_boost_factor_increases_at_small_scales() {
        let halofit = HalofitSpectrum::new(0.3, 0.05, 0.7, 0.8, 0.96);

        let boost_large = halofit.boost_factor(0.01, 0.0);
        let boost_small = halofit.boost_factor(1.0, 0.0);

        // Non-linear effects stronger at small scales
        assert!(boost_small >= boost_large,
                "boost_small = {}, boost_large = {}", boost_small, boost_large);
    }

    #[test]
    fn test_nonlinear_power_positive() {
        let halofit = HalofitSpectrum::new(0.3, 0.05, 0.7, 0.8, 0.96);
        let p_nl = halofit.nonlinear_power(0.1, 0.0);
        assert!(p_nl > 0.0);
    }
}
