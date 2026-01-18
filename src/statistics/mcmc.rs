//! Markov Chain Monte Carlo sampling

use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Parameter for MCMC
#[derive(Debug, Clone)]
pub struct Parameter {
    pub name: String,
    pub initial: f64,
    pub min: f64,
    pub max: f64,
    pub proposal_width: f64,
}

/// MCMC chain
pub struct Chain {
    pub samples: Vec<Vec<f64>>,
    pub log_probs: Vec<f64>,
    pub parameter_names: Vec<String>,
}

impl Chain {
    /// Get samples for parameter
    pub fn get_parameter(&self, name: &str) -> Option<Vec<f64>> {
        self.parameter_names
            .iter()
            .position(|n| n == name)
            .map(|idx| self.samples.iter().map(|s| s[idx]).collect())
    }

    /// Mean of parameter
    pub fn mean(&self, name: &str) -> Option<f64> {
        self.get_parameter(name)
            .map(|samples| samples.iter().sum::<f64>() / samples.len() as f64)
    }

    /// Standard deviation
    pub fn std(&self, name: &str) -> Option<f64> {
        self.get_parameter(name).map(|samples| {
            let mean = samples.iter().sum::<f64>() / samples.len() as f64;
            let var =
                samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
            var.sqrt()
        })
    }

    /// Percentile
    pub fn percentile(&self, name: &str, p: f64) -> Option<f64> {
        self.get_parameter(name).map(|mut samples| {
            samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let idx = ((samples.len() - 1) as f64 * p / 100.0) as usize;
            samples[idx]
        })
    }
}

/// MCMC sampler
pub struct MCMCSampler<F>
where
    F: Fn(&[f64]) -> f64,
{
    pub parameters: Vec<Parameter>,
    pub log_likelihood: F,
    pub n_walkers: usize,
    pub n_steps: usize,
}

impl<F> MCMCSampler<F>
where
    F: Fn(&[f64]) -> f64,
{
    /// Create new sampler
    pub fn new(
        parameters: Vec<Parameter>,
        log_likelihood: F,
        n_walkers: usize,
        n_steps: usize,
    ) -> Self {
        MCMCSampler {
            parameters,
            log_likelihood,
            n_walkers,
            n_steps,
        }
    }

    /// Run MCMC
    #[allow(deprecated)]
    pub fn run(&self, burn_in: usize) -> Chain {
        let mut rng = rand::thread_rng();

        // Initialize walkers
        let mut walkers: Vec<Vec<f64>> = (0..self.n_walkers)
            .map(|_| {
                self.parameters
                    .iter()
                    .map(|p| p.initial + rng.gen_range(-0.1..0.1) * p.proposal_width)
                    .collect()
            })
            .collect();

        let mut walker_log_probs: Vec<f64> =
            walkers.iter().map(|w| (self.log_likelihood)(w)).collect();

        let mut chain = Chain {
            samples: Vec::new(),
            log_probs: Vec::new(),
            parameter_names: self.parameters.iter().map(|p| p.name.clone()).collect(),
        };

        // MCMC iterations
        for step in 0..self.n_steps {
            for walker_idx in 0..self.n_walkers {
                // Propose new position
                let mut proposed = walkers[walker_idx].clone();

                for (i, param) in self.parameters.iter().enumerate() {
                    let normal = Normal::new(0.0, param.proposal_width).unwrap();
                    let delta = normal.sample(&mut rng);
                    proposed[i] += delta;

                    // Reflect at boundaries
                    if proposed[i] < param.min {
                        proposed[i] = 2.0 * param.min - proposed[i];
                    }
                    if proposed[i] > param.max {
                        proposed[i] = 2.0 * param.max - proposed[i];
                    }
                }

                let proposed_log_prob = (self.log_likelihood)(&proposed);

                // Metropolis-Hastings acceptance
                let log_ratio = proposed_log_prob - walker_log_probs[walker_idx];
                if log_ratio > 0.0 || rng.gen::<f64>().ln() < log_ratio {
                    walkers[walker_idx] = proposed;
                    walker_log_probs[walker_idx] = proposed_log_prob;
                }
            }

            // Store samples (after burn-in)
            if step >= burn_in {
                for (walker, &log_prob) in walkers.iter().zip(walker_log_probs.iter()) {
                    chain.samples.push(walker.clone());
                    chain.log_probs.push(log_prob);
                }
            }

            if step % 100 == 0 {
                let mean_log_prob = walker_log_probs.iter().sum::<f64>() / self.n_walkers as f64;
                println!("Step {}: <log L> = {:.2}", step, mean_log_prob);
            }
        }

        chain
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mcmc_gaussian() {
        // Test with simple Gaussian likelihood
        let params = vec![Parameter {
            name: "x".to_string(),
            initial: 0.0,
            min: -5.0,
            max: 5.0,
            proposal_width: 0.5,
        }];

        let log_likelihood = |theta: &[f64]| {
            let x = theta[0];
            -0.5 * x * x // Gaussian centered at 0
        };

        let sampler = MCMCSampler::new(params, log_likelihood, 10, 100);
        let chain = sampler.run(20);

        // Check we got samples
        assert!(!chain.samples.is_empty());

        // Mean should be close to 0
        let mean = chain.mean("x").unwrap();
        assert!(mean.abs() < 1.0);
    }

    #[test]
    fn test_chain_statistics() {
        let chain = Chain {
            samples: vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]],
            log_probs: vec![-1.0, -2.0, -3.0, -4.0, -5.0],
            parameter_names: vec!["x".to_string()],
        };

        let mean = chain.mean("x").unwrap();
        assert!((mean - 3.0).abs() < 1e-10);

        let std = chain.std("x").unwrap();
        assert!((std - (2.0_f64).sqrt()).abs() < 0.01);
    }
}
