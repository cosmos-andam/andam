//! Halo mass function and bias

use crate::constants::*;
use std::f64::consts::PI;

/// Halo mass function formalism
#[derive(Debug, Clone, Copy)]
pub enum MassFunctionType {
    PressSchechter,
    ShethTormen,
    Tinker,
}

/// Halo mass function calculator
pub struct HaloMassFunction {
    pub omega_m: f64,
    pub sigma_8: f64,
    pub mf_type: MassFunctionType,
}

impl HaloMassFunction {
    /// Create new mass function calculator
    pub fn new(omega_m: f64, sigma_8: f64, mf_type: MassFunctionType) -> Self {
        HaloMassFunction {
            omega_m,
            sigma_8,
            mf_type,
        }
    }

    /// RMS fluctuation σ(M) at mass scale M
    pub fn sigma_mass(&self, mass_msun: f64, z: f64) -> f64 {
        // σ(M) ∝ M^(-1/3) approximately
        // Normalized so σ(M*) = σ_8 at M* ≈ 10^13 M_sun
        let m_star = 1e13;
        let sigma_star = self.sigma_8;

        // Power law scaling
        let sigma = sigma_star * (mass_msun / m_star).powf(-1.0 / 3.0);

        // Redshift evolution: σ(z) = σ(0) / D(z)
        // Simplified: D(z) ≈ 1/(1+z) in matter era
        sigma / (1.0 + z)
    }

    /// Peak height ν = δ_c / σ(M)
    pub fn peak_height(&self, mass_msun: f64, z: f64) -> f64 {
        let delta_c = 1.686; // Critical overdensity for collapse
        let sigma = self.sigma_mass(mass_msun, z);
        delta_c / sigma
    }

    /// Multiplicity function f(ν)
    pub fn multiplicity_function(&self, nu: f64) -> f64 {
        match self.mf_type {
            MassFunctionType::PressSchechter => {
                // f(ν) = √(2/π) ν exp(-ν²/2)
                (2.0 / PI).sqrt() * nu * (-nu.powi(2) / 2.0).exp()
            }
            MassFunctionType::ShethTormen => {
                // Sheth-Tormen (1999)
                let a = 0.707;
                let p = 0.3;
                let a_nu2 = a * nu.powi(2);

                let a_factor = (2.0 * a / PI).sqrt();
                a_factor * (1.0 + a_nu2.powf(-p)) * (a_nu2).sqrt() * (-a_nu2 / 2.0).exp()
            }
            MassFunctionType::Tinker => {
                // Tinker et al. (2008) - simplified
                let alpha = 0.368;
                let beta = 0.589;
                let gamma = 0.864;
                let phi = -0.729;

                alpha
                    * (1.0 + (beta / nu).powf(2.0 * phi))
                    * nu.powf(2.0 * beta)
                    * (-gamma * nu.powi(2) / 2.0).exp()
            }
        }
    }

    /// Halo mass function dn/dM [h⁴ Mpc⁻³ M_sun⁻¹]
    pub fn dn_dm(&self, mass_msun: f64, z: f64) -> f64 {
        // dn/dM = (ρ̄/M) f(ν) |dlnσ/dlnM|
        let rho_crit = critical_density(70.0); // kg/m³
        let rho_m = self.omega_m * rho_crit;
        let rho_m_msun_mpc3 = rho_m * (PARSEC * 1e6).powi(3) / M_SUN;

        let nu = self.peak_height(mass_msun, z);
        let f_nu = self.multiplicity_function(nu);

        let sigma = self.sigma_mass(mass_msun, z);
        let dln_sigma_dln_m = 1.0 / 3.0; // For σ ∝ M^(-1/3)

        (rho_m_msun_mpc3 / mass_msun) * f_nu * (dln_sigma_dln_m / sigma)
    }

    /// Number density of halos in mass range [M_min, M_max]
    pub fn number_density(&self, m_min: f64, m_max: f64, z: f64, n_bins: usize) -> f64 {
        let log_m_min = m_min.ln();
        let log_m_max = m_max.ln();
        let dlogm = (log_m_max - log_m_min) / n_bins as f64;

        let mut n_total = 0.0;

        for i in 0..n_bins {
            let log_m = log_m_min + (i as f64 + 0.5) * dlogm;
            let m = log_m.exp();
            let dn_dm = self.dn_dm(m, z);
            n_total += dn_dm * m * dlogm;
        }

        n_total
    }

    /// Halo bias b(M)
    pub fn halo_bias(&self, mass_msun: f64, z: f64) -> f64 {
        let nu = self.peak_height(mass_msun, z);
        let delta_c = 1.686;

        match self.mf_type {
            MassFunctionType::PressSchechter => {
                // b = 1 + (ν² - 1)/δ_c
                1.0 + (nu.powi(2) - 1.0) / delta_c
            }
            MassFunctionType::ShethTormen => {
                // Sheth-Tormen bias
                let a = 0.707;
                let b_param = 0.5;
                let c = 0.6;

                let a_nu2 = a * nu.powi(2);
                1.0 + (a_nu2 - 1.0) / delta_c + 2.0 * b_param / (delta_c * (1.0 + a_nu2.powf(c)))
            }
            MassFunctionType::Tinker => {
                // Tinker et al. (2010) simplified bias
                let y = (nu.powi(2)).max(0.01);
                1.0 + 0.24 * y.ln() + 0.4 * nu.powi(2)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_halo_mass_function() {
        let hmf = HaloMassFunction::new(0.3, 0.8, MassFunctionType::ShethTormen);

        let m = 1e14; // M_sun
        let z = 0.0;

        let dn_dm = hmf.dn_dm(m, z);
        assert!(dn_dm > 0.0);
    }

    #[test]
    fn test_halo_bias_increases_with_mass() {
        let hmf = HaloMassFunction::new(0.3, 0.8, MassFunctionType::ShethTormen);

        let b1 = hmf.halo_bias(1e12, 0.0);
        let b2 = hmf.halo_bias(1e15, 0.0);

        // More massive halos should be more biased
        assert!(b2 > b1, "b1 = {}, b2 = {}", b1, b2);
    }

    #[test]
    fn test_sigma_mass_decreases_with_mass() {
        let hmf = HaloMassFunction::new(0.3, 0.8, MassFunctionType::ShethTormen);

        let s1 = hmf.sigma_mass(1e12, 0.0);
        let s2 = hmf.sigma_mass(1e14, 0.0);

        assert!(s2 < s1);
    }
}
