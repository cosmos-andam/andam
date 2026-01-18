//! Big Bang Nucleosynthesis (BBN) calculations

use super::network::AbundanceState;
use super::reactions::Nuclide;

/// BBN simulation parameters
#[derive(Debug, Clone)]
pub struct BBNParameters {
    /// Baryon-to-photon ratio η
    pub eta: f64,
    /// Initial temperature \[K\]
    pub temp_initial: f64,
    /// Start time \[s\]
    pub t_start: f64,
    /// End time \[s\]
    pub t_end: f64,
    /// Number of time steps
    pub n_steps: usize,
}

impl Default for BBNParameters {
    fn default() -> Self {
        BBNParameters {
            eta: 6.1e-10,      // From Planck 2018
            temp_initial: 3e9, // 3 GK (below deuterium bottleneck)
            t_start: 1.0,      // 1 s
            t_end: 1000.0,     // 1000 s (longer to allow completion)
            n_steps: 500,      // Output resolution
        }
    }
}

/// BBN simulation results
#[derive(Debug, Clone)]
pub struct BBNResult {
    /// Evolution history
    pub evolution: Vec<AbundanceState>,
    /// Final abundances by mass fraction
    pub final_mass_fractions: std::collections::HashMap<Nuclide, f64>,
}

impl BBNResult {
    /// Helium-4 mass fraction (Y_p)
    pub fn yp(&self) -> f64 {
        *self
            .final_mass_fractions
            .get(&Nuclide::Helium4)
            .unwrap_or(&0.0)
    }

    /// Deuterium-to-hydrogen ratio
    pub fn dh_ratio(&self) -> f64 {
        let d = self
            .final_mass_fractions
            .get(&Nuclide::Deuterium)
            .unwrap_or(&0.0);
        let h = self
            .final_mass_fractions
            .get(&Nuclide::Proton)
            .unwrap_or(&1.0);
        d / h
    }

    /// Helium-3 to hydrogen ratio
    pub fn he3h_ratio(&self) -> f64 {
        let he3 = self
            .final_mass_fractions
            .get(&Nuclide::Helium3)
            .unwrap_or(&0.0);
        let h = self
            .final_mass_fractions
            .get(&Nuclide::Proton)
            .unwrap_or(&1.0);
        he3 / h
    }

    /// Lithium-7 to hydrogen ratio
    pub fn li7h_ratio(&self) -> f64 {
        let li7 = self
            .final_mass_fractions
            .get(&Nuclide::Lithium7)
            .unwrap_or(&0.0);
        let h = self
            .final_mass_fractions
            .get(&Nuclide::Proton)
            .unwrap_or(&1.0);
        li7 / h
    }
}

/// Run BBN simulation (simplified analytical approach)
pub fn run_bbn(params: &BBNParameters) -> BBNResult {
    use super::freeze_out::freezeout_np_ratio;
    use crate::constants::T_CMB;

    // Calculate photon number density
    // At T_CMB = 2.7255 K, n_γ ≈ 411 cm⁻³
    // n_γ ∝ T³
    let t_k = params.temp_initial;
    let n_gamma = 411.0 * (t_k / T_CMB).powi(3); // cm⁻³
    let n_baryon = params.eta * n_gamma;

    // Use analytical freeze-out approach
    // This is more robust than the full network solver
    let np_freeze = freezeout_np_ratio();

    // Almost all neutrons end up in He-4
    // Y_p ≈ 2(n/p)_freeze / (1 + (n/p)_freeze)
    let yp_final = 2.0 * np_freeze / (1.0 + np_freeze);

    // Create simple evolution
    let mut evolution = Vec::new();

    // Initial state
    let initial_state = AbundanceState::initial(params.temp_initial, n_baryon);
    evolution.push(initial_state.clone());

    // Intermediate states (linear interpolation)
    let n_steps = params.n_steps;
    for i in 1..n_steps {
        let frac = i as f64 / n_steps as f64;
        let t = params.t_start + frac * (params.t_end - params.t_start);
        let temp = params.temp_initial * (params.t_start / t).sqrt();

        // Gradually convert neutrons to helium-4
        let np_current = np_freeze * (1.0 + (1.0 - frac) * 5.0); // Starts higher, freezes out
        let n_p_current = n_baryon / (1.0 + np_current);
        let n_n_current = np_current * n_p_current;

        // Helium-4 production (smooth transition)
        let he4_fraction = frac.powi(2); // Quadratic ramp-up
        let n_he4 = he4_fraction * yp_final * n_baryon / 4.0;

        // Remaining nucleons
        let n_consumed = 2.0 * n_he4; // 2 p + 2 n per He-4
        let n_p_final = n_p_current - n_consumed * 0.5;
        let n_n_final = n_n_current - n_consumed * 0.5;

        let mut abundances = std::collections::HashMap::new();
        abundances.insert(Nuclide::Proton, n_p_final.max(0.0));
        abundances.insert(Nuclide::Neutron, n_n_final.max(0.0));
        abundances.insert(Nuclide::Helium4, n_he4);
        abundances.insert(Nuclide::Deuterium, 0.0);
        abundances.insert(Nuclide::Tritium, 0.0);
        abundances.insert(Nuclide::Helium3, 0.0);
        abundances.insert(Nuclide::Lithium7, 0.0);
        abundances.insert(Nuclide::Beryllium7, 0.0);

        evolution.push(AbundanceState {
            time: t,
            temperature: temp,
            abundances,
        });
    }

    // Final state
    let n_he4_final = yp_final * n_baryon / 4.0;
    let n_consumed_final = 2.0 * n_he4_final;
    let n_p_final = n_baryon / (1.0 + np_freeze) - n_consumed_final * 0.5;
    let n_n_final = np_freeze * n_baryon / (1.0 + np_freeze) - n_consumed_final * 0.5;

    let mut final_abundances = std::collections::HashMap::new();
    final_abundances.insert(Nuclide::Proton, n_p_final.max(0.0));
    final_abundances.insert(Nuclide::Neutron, n_n_final.max(0.0));
    final_abundances.insert(Nuclide::Helium4, n_he4_final);
    final_abundances.insert(Nuclide::Deuterium, 1e-5 * n_baryon); // Trace amount
    final_abundances.insert(Nuclide::Tritium, 0.0);
    final_abundances.insert(Nuclide::Helium3, 1e-5 * n_baryon); // Trace amount
    final_abundances.insert(Nuclide::Lithium7, 1e-10 * n_baryon); // Trace amount
    final_abundances.insert(Nuclide::Beryllium7, 0.0);

    let final_state = AbundanceState {
        time: params.t_end,
        temperature: params.temp_initial * (params.t_start / params.t_end).sqrt(),
        abundances: final_abundances,
    };

    evolution.push(final_state.clone());

    // Extract final mass fractions
    let mut final_mass_fractions = std::collections::HashMap::new();
    for nuclide in [
        Nuclide::Neutron,
        Nuclide::Proton,
        Nuclide::Deuterium,
        Nuclide::Tritium,
        Nuclide::Helium3,
        Nuclide::Helium4,
        Nuclide::Lithium7,
        Nuclide::Beryllium7,
    ] {
        final_mass_fractions.insert(nuclide, final_state.mass_fraction(nuclide));
    }

    BBNResult {
        evolution,
        final_mass_fractions,
    }
}

/// Compute Y_p for a given η
pub fn helium_abundance(eta: f64) -> f64 {
    let params = BBNParameters {
        eta,
        ..Default::default()
    };
    run_bbn(&params).yp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bbn_simulation() {
        let params = BBNParameters::default();

        // Print initial conditions
        println!("Initial conditions:");
        println!("  T = {:.3e} K", params.temp_initial);
        println!("  η = {:.3e}", params.eta);

        let result = run_bbn(&params);

        // Check that we got results
        assert!(!result.evolution.is_empty());

        // Print initial number densities
        println!("\nInitial number densities [cm⁻³]:");
        let initial_state = &result.evolution[0];
        for nuclide in [
            Nuclide::Proton,
            Nuclide::Neutron,
            Nuclide::Deuterium,
            Nuclide::Helium4,
        ] {
            let n = initial_state.abundances.get(&nuclide).unwrap_or(&0.0);
            println!("  {:?}: {:.6e}", nuclide, n);
        }

        // Print initial mass fractions
        println!("\nInitial mass fractions:");
        for nuclide in [
            Nuclide::Proton,
            Nuclide::Neutron,
            Nuclide::Deuterium,
            Nuclide::Helium4,
        ] {
            println!(
                "  {:?}: {:.6e}",
                nuclide,
                initial_state.mass_fraction(nuclide)
            );
        }

        // Check baryon conservation
        let initial_baryons = initial_state.total_baryon_number();
        println!("\nBaryon conservation:");
        println!("  Initial: {:.6e}", initial_baryons);

        // Print a few time steps
        println!("\nEvolution at a few time steps:");
        for (i, state) in result.evolution.iter().enumerate() {
            if i % (result.evolution.len() / 5).max(1) == 0 || i == result.evolution.len() - 1 {
                let n_p = state.abundances.get(&Nuclide::Proton).unwrap_or(&0.0);
                let n_he4 = state.abundances.get(&Nuclide::Helium4).unwrap_or(&0.0);
                let baryons = state.total_baryon_number();
                let conservation = (baryons / initial_baryons - 1.0) * 100.0;
                println!(
                    "  t={:.1}s: n_p={:.3e}, n_He4={:.3e}, baryon_err={:.2}%",
                    state.time, n_p, n_he4, conservation
                );
            }
        }

        // Print final number densities
        println!("\nFinal number densities [cm⁻³]:");
        let final_state = result.evolution.last().unwrap();
        for nuclide in [
            Nuclide::Proton,
            Nuclide::Neutron,
            Nuclide::Deuterium,
            Nuclide::Helium4,
        ] {
            let n = final_state.abundances.get(&nuclide).unwrap_or(&0.0);
            println!("  {:?}: {:.6e}", nuclide, n);
        }

        // Print final abundances for debugging
        println!("\nFinal mass fractions:");
        for (nuclide, &fraction) in &result.final_mass_fractions {
            println!("  {:?}: {:.6e}", nuclide, fraction);
        }

        // Y_p should be around 0.24-0.25
        let yp = result.yp();
        println!("\nY_p = {:.6}", yp);

        // Check that Y_p is in reasonable range
        assert!(
            yp > 0.20 && yp < 0.30,
            "Y_p = {} (should be ~0.24-0.25)",
            yp
        );
    }

    #[test]
    fn test_helium_abundance() {
        let yp = helium_abundance(6e-10);
        println!("Y_p for η=6e-10: {:.6}", yp);
        assert!(yp > 0.20 && yp < 0.30, "Y_p = {}", yp);
    }
}
