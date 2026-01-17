//! Simplified Boltzmann equation framework for cosmological perturbations

use crate::dynamics::Universe;

/// Perturbation mode (scalar, vector, or tensor)
#[derive(Debug, Clone, Copy)]
pub enum PerturbationMode {
    Scalar,
    Vector,
    Tensor,
}

/// State of perturbation variables (simplified)
#[derive(Debug, Clone)]
pub struct PerturbationState {
    /// Scale factor
    pub a: f64,
    /// Photon monopole Θ_0
    pub theta_0: f64,
    /// Photon dipole Θ_1
    pub theta_1: f64,
    /// Photon quadrupole Θ_2
    pub theta_2: f64,
}

/// Boltzmann equation solver (simplified version)
pub struct BoltzmannSolver {
    /// Universe model
    pub universe: Universe,
    /// Wavenumber k [Mpc⁻¹]
    pub k: f64,
    /// Perturbation mode
    pub mode: PerturbationMode,
    /// Maximum multipole
    pub l_max: usize,
}

impl BoltzmannSolver {
    /// Create new Boltzmann solver
    pub fn new(universe: Universe, k: f64, l_max: usize) -> Self {
        BoltzmannSolver {
            universe,
            k,
            mode: PerturbationMode::Scalar,
            l_max,
        }
    }

    /// Initialize perturbations
    pub fn initial_conditions(&self) -> PerturbationState {
        PerturbationState {
            a: 1e-5,
            theta_0: -0.5,
            theta_1: 0.1,
            theta_2: 0.0,
        }
    }

    /// Evolve perturbations (simplified)
    pub fn evolve(&mut self, a_init: f64, a_final: f64) -> Vec<PerturbationState> {
        let mut results = Vec::new();
        let n_steps = 100;
        let dlna = (a_final / a_init).ln() / n_steps as f64;

        let mut state = self.initial_conditions();
        state.a = a_init;
        results.push(state.clone());

        for _ in 1..=n_steps {
            state.a *= dlna.exp();
            // Simplified evolution
            state.theta_0 += -0.001 * state.theta_1;
            state.theta_1 += 0.0005 * (state.theta_0 - 2.0 * state.theta_2);
            state.theta_2 += 0.0001 * state.theta_1;

            results.push(state.clone());
        }

        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boltzmann_initialization() {
        let universe = crate::dynamics::Universe::benchmark();
        let solver = BoltzmannSolver::new(universe, 0.1, 10);
        let state = solver.initial_conditions();

        assert!(state.theta_0.abs() > 0.0);
    }

    #[test]
    fn test_boltzmann_evolution() {
        let universe = crate::dynamics::Universe::benchmark();
        let mut solver = BoltzmannSolver::new(universe, 0.01, 10);
        let results = solver.evolve(1e-5, 1e-2);

        assert!(results.len() > 10);
    }
}
