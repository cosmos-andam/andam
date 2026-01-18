//! Nuclear reaction network solver

use super::reactions::*;
use crate::constants::*;
use std::collections::HashMap;

/// Abundance evolution state
#[derive(Debug, Clone)]
pub struct AbundanceState {
    pub time: f64,           // Time [s]
    pub temperature: f64,    // Temperature [K]
    pub abundances: HashMap<Nuclide, f64>, // Number densities [cm⁻³]
}

impl AbundanceState {
    /// Create initial state
    pub fn initial(temp_k: f64, baryon_density: f64) -> Self {
        let mut abundances = HashMap::new();

        // Initial conditions: all baryons in p and n
        // n/p ratio from equilibrium at high temperature
        let n_to_p = neutron_proton_equilibrium(temp_k);

        let n_p = baryon_density / (1.0 + n_to_p);
        let n_n = n_to_p * n_p;

        abundances.insert(Nuclide::Proton, n_p);
        abundances.insert(Nuclide::Neutron, n_n);
        abundances.insert(Nuclide::Deuterium, 0.0);
        abundances.insert(Nuclide::Tritium, 0.0);
        abundances.insert(Nuclide::Helium3, 0.0);
        abundances.insert(Nuclide::Helium4, 0.0);
        abundances.insert(Nuclide::Lithium7, 0.0);
        abundances.insert(Nuclide::Beryllium7, 0.0);

        AbundanceState {
            time: 0.0,
            temperature: temp_k,
            abundances,
        }
    }

    /// Total baryon number (should be conserved)
    pub fn total_baryon_number(&self) -> f64 {
        self.abundances.iter()
            .map(|(nuclide, &n)| n * nuclide.mass_number() as f64)
            .sum()
    }

    /// Mass fraction of a nuclide
    pub fn mass_fraction(&self, nuclide: Nuclide) -> f64 {
        let n = self.abundances.get(&nuclide).copied().unwrap_or(0.0);
        let a = nuclide.mass_number() as f64;

        let total_mass: f64 = self.abundances.iter()
            .map(|(nuc, &num)| num * nuc.mass_number() as f64)
            .sum();

        if total_mass > 0.0 {
            n * a / total_mass
        } else {
            0.0
        }
    }

    /// Helium-4 mass fraction Y_p
    pub fn helium_mass_fraction(&self) -> f64 {
        self.mass_fraction(Nuclide::Helium4)
    }
}

/// Neutron-to-proton ratio from weak equilibrium
fn neutron_proton_equilibrium(temp_k: f64) -> f64 {
    let delta_m = (M_N - M_P) * C * C; // Mass difference [J]
    let delta_m_mev = delta_m * J_TO_EV * 1e-6; // [MeV]
    let t_mev = K_B * temp_k * J_TO_EV * 1e-6; // Temperature [MeV]

    (-delta_m_mev / t_mev).exp()
}

/// Nuclear reaction network solver
pub struct NetworkSolver {
    reactions: Vec<Reaction>,
}

impl NetworkSolver {
    /// Create new network solver
    pub fn new(reactions: Vec<Reaction>) -> Self {
        NetworkSolver { reactions }
    }

    /// Standard BBN network
    pub fn standard_bbn() -> Self {
        NetworkSolver::new(standard_bbn_reactions())
    }

    /// Compute time derivatives of abundances
    pub fn derivatives(&self, state: &AbundanceState) -> HashMap<Nuclide, f64> {
        let mut dndt = HashMap::new();

        // Initialize all derivatives to zero
        for nuclide in [
            Nuclide::Neutron, Nuclide::Proton, Nuclide::Deuterium,
            Nuclide::Tritium, Nuclide::Helium3, Nuclide::Helium4,
            Nuclide::Lithium7, Nuclide::Beryllium7,
        ] {
            dndt.insert(nuclide, 0.0);
        }

        // Add contributions from each reaction
        for reaction in &self.reactions {
            let rate = reaction.rate_coefficient(state.temperature);

            // Get reactant abundances
            let mut n_product = rate;
            for reactant in &reaction.reactants {
                let n = state.abundances.get(reactant).copied().unwrap_or(0.0);
                n_product *= n;
            }

            // Decrease reactants
            for reactant in &reaction.reactants {
                *dndt.get_mut(reactant).unwrap() -= n_product;
            }

            // Increase products
            for product in &reaction.products {
                *dndt.get_mut(product).unwrap() += n_product;
            }

            // Add reverse reaction
            let reverse_rate = reaction.reverse_rate(state.temperature);
            let mut n_reverse = reverse_rate;
            for product in &reaction.products {
                let n = state.abundances.get(product).copied().unwrap_or(0.0);
                n_reverse *= n;
            }

            for product in &reaction.products {
                *dndt.get_mut(product).unwrap() -= n_reverse;
            }
            for reactant in &reaction.reactants {
                *dndt.get_mut(reactant).unwrap() += n_reverse;
            }
        }

        dndt
    }

    /// Evolve network from t_start to t_end with adaptive time stepping
    pub fn evolve(
        &self,
        initial_state: AbundanceState,
        t_start: f64,
        t_end: f64,
        n_steps: usize,
    ) -> Vec<AbundanceState> {
        let mut states = vec![initial_state.clone()];

        // Save initial temperature for temperature scaling
        let initial_temp = initial_state.temperature;

        let mut current = initial_state;
        let mut t = t_start;

        // Adaptive time stepping
        let dt_max = 10.0; // Maximum step size (10 seconds)
        let dt_min = 1e-6; // Minimum step (1 microsecond)
        let mut dt: f64; // Will be calculated adaptively each step

        let mut step_count = 0;
        const MAX_STEPS: usize = 1000000; // Safety limit (increased)

        while t < t_end && step_count < MAX_STEPS {
            // Calculate derivatives
            let derivs = self.derivatives(&current);

            // Estimate appropriate time step based on rate of change
            let mut dt_safe: f64 = dt_max;
            for (nuclide, &dndt) in &derivs {
                let n = current.abundances.get(nuclide).copied().unwrap_or(0.0);
                if n > 1e10 && dndt.abs() > 1e-10 {
                    // Don't let abundance change by more than 10% per step
                    let dt_suggested = 0.1 * n / dndt.abs();
                    dt_safe = dt_safe.min(dt_suggested);
                }
            }

            // Clamp time step
            dt = dt_safe.max(dt_min).min(dt_max);

            // Don't overshoot end time
            if t + dt > t_end {
                dt = t_end - t;
            }

            // RK2 (midpoint method) for better stability
            // Step 1: Half step with Euler
            let mut mid_abundances = current.abundances.clone();
            for (nuclide, &dndt) in &derivs {
                let old_n = current.abundances.get(nuclide).copied().unwrap_or(0.0);
                let mid_n = (old_n + dndt * dt / 2.0).max(0.0);
                mid_abundances.insert(*nuclide, mid_n);
            }

            let mid_time = t + dt / 2.0;
            let mid_temp = initial_temp * (t_start / mid_time).sqrt();
            let mid_state = AbundanceState {
                time: mid_time,
                temperature: mid_temp,
                abundances: mid_abundances,
            };

            // Step 2: Full step using midpoint derivative
            let mid_derivs = self.derivatives(&mid_state);
            let mut new_abundances = current.abundances.clone();
            for (nuclide, &dndt) in &mid_derivs {
                let old_n = current.abundances.get(nuclide).copied().unwrap_or(0.0);
                let new_n = (old_n + dndt * dt).max(0.0);
                new_abundances.insert(*nuclide, new_n);
            }

            // Update state
            t += dt;
            let new_temp = initial_temp * (t_start / t).sqrt();

            current = AbundanceState {
                time: t,
                temperature: new_temp,
                abundances: new_abundances,
            };

            // Save state periodically
            if step_count % (MAX_STEPS / n_steps).max(1) == 0 || t >= t_end {
                states.push(current.clone());
            }

            step_count += 1;
        }

        // Ensure we have the final state
        if states.last().map(|s| s.time).unwrap_or(0.0) < t_end {
            states.push(current);
        }

        states
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_network_evolution() {
        let solver = NetworkSolver::standard_bbn();
        let initial = AbundanceState::initial(1e9, 1e18);

        let evolution = solver.evolve(initial, 0.1, 100.0, 100);

        // Should have at least some evolution states
        assert!(evolution.len() >= 2, "Evolution should have at least initial and final states");
        assert!(evolution.len() <= 100000, "Evolution should not exceed safety limit");
    }

    #[test]
    fn test_initial_state() {
        let state = AbundanceState::initial(1e10, 1e18);
        let n_p = state.abundances.get(&Nuclide::Proton).unwrap();
        let n_n = state.abundances.get(&Nuclide::Neutron).unwrap();

        assert!(n_p > &0.0);
        assert!(n_n > &0.0);
    }
}
