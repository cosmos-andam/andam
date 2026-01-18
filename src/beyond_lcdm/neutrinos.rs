//! Massive neutrino cosmology

/// Neutrino mass hierarchy
#[derive(Debug, Clone, Copy)]
pub enum MassHierarchy {
    Normal,
    Inverted,
    Degenerate,
}

/// Massive neutrino component
pub struct MassiveNeutrinos {
    pub total_mass: f64, // Sum of masses [eV]
    pub hierarchy: MassHierarchy,
    pub n_species: usize,
}

impl MassiveNeutrinos {
    /// Create neutrino model
    pub fn new(total_mass: f64, hierarchy: MassHierarchy) -> Self {
        MassiveNeutrinos {
            total_mass,
            hierarchy,
            n_species: 3,
        }
    }

    /// Neutrino density parameter
    pub fn omega_nu(&self, h: f64) -> f64 {
        // Ω_ν h² = Σm_ν / (93.14 eV)
        self.total_mass / 93.14 / (h * h)
    }

    /// Suppression of matter power spectrum
    pub fn power_suppression(&self, k: f64, _z: f64) -> f64 {
        // Simplified: exp(-(k/k_fs)^α) where k_fs depends on neutrino mass
        // For m_ν ~ 0.06 eV, k_fs ~ 0.05-0.1 h/Mpc
        let k_fs = 0.3 * (self.total_mass / 0.06); // Free-streaming scale [h/Mpc]
        let alpha = 1.5; // Power-law suppression
        (-(k / k_fs).powf(alpha)).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neutrino_omega() {
        let neutrinos = MassiveNeutrinos::new(0.06, MassHierarchy::Normal);
        let omega_nu = neutrinos.omega_nu(0.7);

        // Should be small
        assert!(omega_nu < 0.01);
        assert!(omega_nu > 0.0);
    }

    #[test]
    fn test_power_suppression() {
        let neutrinos = MassiveNeutrinos::new(0.06, MassHierarchy::Normal);

        // At small k, minimal suppression
        let s_small = neutrinos.power_suppression(0.01, 0.0);
        assert!(s_small > 0.9);

        // At large k, significant suppression
        let s_large = neutrinos.power_suppression(1.0, 0.0);
        assert!(s_large < 0.5);
    }
}
