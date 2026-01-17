//! Weak gravitational lensing calculations

use ndarray::Array2;
use crate::observations::comoving_distance;
use crate::structure::power_spectrum::matter_power_spectrum;
use crate::dynamics::Universe;

/// Convergence field κ
#[derive(Debug, Clone)]
pub struct ConvergenceField {
    pub kappa: Array2<f64>,
    pub size: usize,
}

impl ConvergenceField {
    /// Create new convergence field
    pub fn new(size: usize) -> Self {
        ConvergenceField {
            kappa: Array2::zeros((size, size)),
            size,
        }
    }

    /// Compute convergence from matter power spectrum (simplified)
    pub fn from_power_spectrum(
        universe: &Universe,
        source_z: f64,
        size: usize,
        field_size_deg: f64,
    ) -> Self {
        let mut field = Self::new(size);

        let _chi_s = comoving_distance(source_z, universe);

        // Simplified field generation
        for i in 0..size {
            for j in 0..size {
                let x = (i as f64 / size as f64 - 0.5) * field_size_deg;
                let y = (j as f64 / size as f64 - 0.5) * field_size_deg;

                // Simple sinusoidal pattern as placeholder
                field.kappa[[i, j]] = 0.01 * (x * 10.0).sin() * (y * 10.0).cos();
            }
        }

        field
    }
}

/// Shear components (γ1, γ2)
#[derive(Debug, Clone, Copy)]
pub struct Shear {
    pub gamma1: f64,
    pub gamma2: f64,
}

impl Shear {
    /// Magnitude of shear
    pub fn magnitude(&self) -> f64 {
        (self.gamma1 * self.gamma1 + self.gamma2 * self.gamma2).sqrt()
    }

    /// Position angle
    pub fn angle(&self) -> f64 {
        0.5 * self.gamma2.atan2(self.gamma1)
    }
}

/// Compute shear from convergence using Kaiser-Squires inversion (simplified)
pub fn convergence_to_shear(kappa: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
    let size = kappa.nrows();
    let mut gamma1 = Array2::zeros((size, size));
    let mut gamma2 = Array2::zeros((size, size));

    // Simplified version using finite differences
    for i in 1..size-1 {
        for j in 1..size-1 {
            // Second derivatives
            let d2_dx2 = kappa[[i+1, j]] - 2.0 * kappa[[i, j]] + kappa[[i-1, j]];
            let d2_dy2 = kappa[[i, j+1]] - 2.0 * kappa[[i, j]] + kappa[[i, j-1]];
            let d2_dxdy = (kappa[[i+1, j+1]] - kappa[[i+1, j-1]]
                         - kappa[[i-1, j+1]] + kappa[[i-1, j-1]]) / 4.0;

            gamma1[[i, j]] = d2_dx2 - d2_dy2;
            gamma2[[i, j]] = 2.0 * d2_dxdy;
        }
    }

    (gamma1, gamma2)
}

/// Lensing power spectrum C_ℓ^κκ (simplified)
pub fn lensing_power_spectrum(
    l: usize,
    source_z: f64,
    universe: &Universe,
) -> f64 {
    let chi_s = comoving_distance(source_z, universe);

    // Limber approximation
    let k = (l as f64) / chi_s;

    // Simplified weight function
    let omega_m = 0.3111;
    let h = universe.h0 / 100.0;
    let weight = 1.5 * omega_m * (100.0 * h * h) * chi_s * 0.5;

    let p_k = matter_power_spectrum(k, source_z, omega_m, 0.049, h, 2.1e-9, 0.9665);

    weight * weight * p_k / (chi_s * chi_s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convergence_field() {
        let field = ConvergenceField::new(64);
        assert_eq!(field.size, 64);
        assert_eq!(field.kappa.shape(), &[64, 64]);
    }

    #[test]
    fn test_shear_magnitude() {
        let shear = Shear { gamma1: 0.03, gamma2: 0.04 };
        let mag = shear.magnitude();
        assert!((mag - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_convergence_to_shear() {
        let kappa = Array2::zeros((10, 10));
        let (gamma1, gamma2) = convergence_to_shear(&kappa);
        assert_eq!(gamma1.shape(), &[10, 10]);
        assert_eq!(gamma2.shape(), &[10, 10]);
    }
}
