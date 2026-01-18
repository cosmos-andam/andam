//! CMB polarization calculations

/// Stokes parameters for polarization
#[derive(Debug, Clone, Copy)]
pub struct StokesParameters {
    pub q: f64, // Linear polarization (Q)
    pub u: f64, // Linear polarization (U)
}

impl StokesParameters {
    /// Polarization fraction P = √(Q² + U²)/I
    pub fn polarization_fraction(&self, intensity: f64) -> f64 {
        ((self.q.powi(2) + self.u.powi(2)).sqrt()) / intensity
    }

    /// Polarization angle χ = 0.5 arctan(U/Q)
    pub fn angle(&self) -> f64 {
        0.5 * self.u.atan2(self.q)
    }
}

/// E and B mode decomposition
pub fn decompose_eb(q_map: &[f64], u_map: &[f64], _n_side: usize) -> (Vec<f64>, Vec<f64>) {
    // In full implementation, use spin-weighted spherical harmonics
    // This is simplified

    let mut e_map = vec![0.0; q_map.len()];
    let mut b_map = vec![0.0; u_map.len()];

    e_map.copy_from_slice(q_map);
    b_map.copy_from_slice(u_map);

    (e_map, b_map)
}

/// E-mode and B-mode power spectra
pub struct PolarizationSpectrum {
    pub l_max: usize,
    pub c_l_ee: Vec<f64>, // E-mode auto
    pub c_l_bb: Vec<f64>, // B-mode auto
    pub c_l_te: Vec<f64>, // Temperature-E cross
}

impl PolarizationSpectrum {
    /// Create new polarization spectrum
    pub fn new(l_max: usize) -> Self {
        PolarizationSpectrum {
            l_max,
            c_l_ee: vec![0.0; l_max + 1],
            c_l_bb: vec![0.0; l_max + 1],
            c_l_te: vec![0.0; l_max + 1],
        }
    }

    /// Compute from Boltzmann evolution
    pub fn from_boltzmann(_universe: &crate::dynamics::Universe, l_max: usize) -> Self {
        let mut spectrum = Self::new(l_max);

        // Compute C_l^EE, C_l^BB, C_l^TE
        // Full implementation requires Boltzmann solver

        for l in 2..=l_max {
            // E-mode power (scalar perturbations)
            spectrum.c_l_ee[l] = 5000.0 * (l as f64 / 1000.0).powi(-1);

            // B-mode power (tensor perturbations only)
            // Tiny for standard inflation, zero for scalars
            spectrum.c_l_bb[l] = 0.001 * (l as f64 / 100.0).powi(-2);

            // TE cross-correlation
            spectrum.c_l_te[l] = -50.0 * (l as f64 / 100.0).powi(-1);
        }

        spectrum
    }

    /// Tensor-to-scalar ratio constraint from B-modes
    pub fn tensor_to_scalar_ratio(&self) -> f64 {
        // r = C_l^BB(tensor) / C_l^TT(scalar)
        // Simplified calculation
        let l_pivot = 80;
        if l_pivot < self.c_l_bb.len() {
            self.c_l_bb[l_pivot] / 5000.0
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stokes_parameters() {
        let stokes = StokesParameters { q: 3.0, u: 4.0 };
        let frac = stokes.polarization_fraction(10.0);
        assert!((frac - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_polarization_spectrum() {
        use crate::dynamics::Universe;
        let universe = Universe::benchmark();
        let spectrum = PolarizationSpectrum::from_boltzmann(&universe, 100);

        // Check that spectra are non-zero
        assert!(spectrum.c_l_ee[50] > 0.0);
        assert!(spectrum.c_l_bb[50] > 0.0);
    }

    #[test]
    fn test_eb_decomposition() {
        let q_map = vec![1.0, 2.0, 3.0];
        let u_map = vec![4.0, 5.0, 6.0];
        let (e_map, b_map) = decompose_eb(&q_map, &u_map, 1);

        assert_eq!(e_map.len(), 3);
        assert_eq!(b_map.len(), 3);
    }
}
