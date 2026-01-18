//! Fisher matrix forecasts

use nalgebra::{DMatrix, DVector};

/// Fisher information matrix
pub struct FisherMatrix {
    pub matrix: DMatrix<f64>,
    pub parameter_names: Vec<String>,
}

impl FisherMatrix {
    /// Create Fisher matrix from derivatives
    pub fn from_derivatives<F>(
        params_fiducial: &[f64],
        param_names: Vec<String>,
        observables_fn: F,
        covariance: &DMatrix<f64>,
    ) -> Self
    where
        F: Fn(&[f64]) -> DVector<f64>,
    {
        let n_params = params_fiducial.len();
        let mut fisher = DMatrix::zeros(n_params, n_params);

        // Compute derivatives numerically
        let mut derivatives = Vec::new();
        for i in 0..n_params {
            let mut params_plus = params_fiducial.to_vec();
            let mut params_minus = params_fiducial.to_vec();
            let delta = params_fiducial[i] * 0.01;

            params_plus[i] += delta;
            params_minus[i] -= delta;

            let obs_plus = observables_fn(&params_plus);
            let obs_minus = observables_fn(&params_minus);

            let deriv = (obs_plus - obs_minus) / (2.0 * delta);
            derivatives.push(deriv);
        }

        // Fisher matrix: F_ij = Σ_α (∂O_α/∂θ_i) C^(-1)_αβ (∂O_β/∂θ_j)
        let cov_inv = covariance.clone().try_inverse().expect("Covariance not invertible");

        for i in 0..n_params {
            for j in 0..n_params {
                let mut sum = 0.0;
                for alpha in 0..derivatives[i].len() {
                    for beta in 0..derivatives[j].len() {
                        sum += derivatives[i][alpha] * cov_inv[(alpha, beta)] * derivatives[j][beta];
                    }
                }
                fisher[(i, j)] = sum;
            }
        }

        FisherMatrix {
            matrix: fisher,
            parameter_names: param_names,
        }
    }

    /// Marginalized 1σ error on parameter i
    pub fn marginalized_error(&self, param_idx: usize) -> f64 {
        let cov = self.covariance_matrix();
        cov[(param_idx, param_idx)].sqrt()
    }

    /// Covariance matrix (inverse of Fisher)
    pub fn covariance_matrix(&self) -> DMatrix<f64> {
        self.matrix.clone().try_inverse().expect("Fisher matrix not invertible")
    }

    /// Correlation matrix
    pub fn correlation_matrix(&self) -> DMatrix<f64> {
        let cov = self.covariance_matrix();
        let n = cov.nrows();
        let mut corr = DMatrix::zeros(n, n);

        for i in 0..n {
            for j in 0..n {
                corr[(i, j)] = cov[(i, j)] / (cov[(i, i)] * cov[(j, j)]).sqrt();
            }
        }

        corr
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fisher_matrix() {
        // Simple test: linear model y = mx + b
        let params_fiducial = vec![2.0, 1.0]; // m, b
        let param_names = vec!["slope".to_string(), "intercept".to_string()];

        // Generate mock observables
        let x_values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let observables_fn = |params: &[f64]| {
            let m = params[0];
            let b = params[1];
            DVector::from_vec(x_values.iter().map(|&x| m * x + b).collect())
        };

        // Covariance (diagonal with variance 1.0)
        let covariance = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0; 5]));

        let fisher = FisherMatrix::from_derivatives(
            &params_fiducial,
            param_names,
            observables_fn,
            &covariance,
        );

        // Check that Fisher matrix is positive definite
        assert!(fisher.matrix[(0, 0)] > 0.0);
        assert!(fisher.matrix[(1, 1)] > 0.0);

        // Check errors are finite
        let error_slope = fisher.marginalized_error(0);
        let error_intercept = fisher.marginalized_error(1);
        assert!(error_slope > 0.0 && error_slope.is_finite());
        assert!(error_intercept > 0.0 && error_intercept.is_finite());
    }
}
