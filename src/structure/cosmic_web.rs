//! 3D cosmic web generation and visualization

use ndarray::Array3;
use rand::Rng;

/// 3D density field
pub struct DensityField {
    pub grid: Array3<f64>,
    pub box_size: f64, // Mpc/h
    pub n_grid: usize,
}

impl DensityField {
    /// Create new density field
    pub fn new(n_grid: usize, box_size: f64) -> Self {
        DensityField {
            grid: Array3::zeros((n_grid, n_grid, n_grid)),
            box_size,
            n_grid,
        }
    }

    /// Generate Gaussian random field (simplified)
    pub fn generate_gaussian(&mut self, _omega_m: f64, _z: f64) {
        let mut rng = rand::thread_rng();

        // Generate uncorrelated Gaussian noise
        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                for k in 0..self.n_grid {
                    self.grid[[i, j, k]] = rng.gen_range(-1.0..1.0);
                }
            }
        }

        // Smooth with Gaussian kernel to create correlations
        self.smooth(2.0);
    }

    /// Smooth field with Gaussian kernel
    fn smooth(&mut self, sigma: f64) {
        let kernel_size = (3.0 * sigma).max(1.0) as usize;
        let mut smoothed = self.grid.clone();

        for i in kernel_size..(self.n_grid.saturating_sub(kernel_size)) {
            for j in kernel_size..(self.n_grid.saturating_sub(kernel_size)) {
                for k in kernel_size..(self.n_grid.saturating_sub(kernel_size)) {
                    let mut sum = 0.0;
                    let mut weight_sum = 0.0;

                    for di in 0..kernel_size {
                        for dj in 0..kernel_size {
                            for dk in 0..kernel_size {
                                let dx = di as f64 - kernel_size as f64 / 2.0;
                                let dy = dj as f64 - kernel_size as f64 / 2.0;
                                let dz = dk as f64 - kernel_size as f64 / 2.0;
                                let r2 = dx * dx + dy * dy + dz * dz;
                                let weight = (-r2 / (2.0 * sigma * sigma)).exp();

                                let ii = i + di - kernel_size / 2;
                                let jj = j + dj - kernel_size / 2;
                                let kk = k + dk - kernel_size / 2;

                                if ii < self.n_grid && jj < self.n_grid && kk < self.n_grid {
                                    sum += self.grid[[ii, jj, kk]] * weight;
                                    weight_sum += weight;
                                }
                            }
                        }
                    }

                    if weight_sum > 0.0 {
                        smoothed[[i, j, k]] = sum / weight_sum;
                    }
                }
            }
        }

        self.grid = smoothed;
    }

    /// Extract particles above density threshold
    pub fn extract_particles(&self, threshold: f64) -> Vec<(f64, f64, f64)> {
        let mut particles = Vec::new();
        let dx = self.box_size / self.n_grid as f64;

        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                for k in 0..self.n_grid {
                    if self.grid[[i, j, k]] > threshold {
                        let x = i as f64 * dx;
                        let y = j as f64 * dx;
                        let z = k as f64 * dx;
                        particles.push((x, y, z));
                    }
                }
            }
        }

        particles
    }

    /// Extract 2D slice for visualization
    pub fn slice(&self, axis: usize, index: usize) -> Vec<Vec<f64>> {
        let mut slice_data = vec![vec![0.0; self.n_grid]; self.n_grid];

        #[allow(clippy::needless_range_loop)]
        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                slice_data[i][j] = match axis {
                    0 => self.grid[[index.min(self.n_grid - 1), i, j]],
                    1 => self.grid[[i, index.min(self.n_grid - 1), j]],
                    2 => self.grid[[i, j, index.min(self.n_grid - 1)]],
                    _ => 0.0,
                };
            }
        }

        slice_data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_field_creation() {
        let field = DensityField::new(32, 100.0);
        assert_eq!(field.n_grid, 32);
        assert_eq!(field.box_size, 100.0);
    }

    #[test]
    fn test_gaussian_generation() {
        let mut field = DensityField::new(16, 50.0);
        field.generate_gaussian(0.3, 0.0);

        // Check that field has been populated
        let mut has_nonzero = false;
        for i in 0..field.n_grid {
            for j in 0..field.n_grid {
                for k in 0..field.n_grid {
                    if field.grid[[i, j, k]].abs() > 1e-10 {
                        has_nonzero = true;
                        break;
                    }
                }
            }
        }
        assert!(has_nonzero);
    }

    #[test]
    fn test_extract_particles() {
        let mut field = DensityField::new(16, 50.0);
        field.generate_gaussian(0.3, 0.0);

        let particles = field.extract_particles(0.0);
        // Should have some particles
        assert!(!particles.is_empty());
    }
}
