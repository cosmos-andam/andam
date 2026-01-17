//! Friedmann equation solver and universe evolution

use crate::constants::*;
use crate::dynamics::components::Component;

/// Represents a cosmological model with multiple components
#[derive(Debug, Clone)]
pub struct Universe {
    /// Hubble constant today [km/s/Mpc]
    pub h0: f64,
    /// Components (matter, radiation, dark energy, curvature)
    pub components: Vec<Component>,
}

impl Universe {
    /// Create a new universe with given Hubble constant
    pub fn new(h0: f64) -> Self {
        Universe {
            h0,
            components: Vec::new(),
        }
    }

    /// Add a component to the universe
    pub fn add_component(&mut self, component: Component) {
        self.components.push(component);
    }

    /// Create a benchmark ΛCDM model (Planck 2018)
    pub fn benchmark() -> Self {
        let mut universe = Universe::new(67.66);
        universe.add_component(Component::matter(0.3111));
        universe.add_component(Component::radiation(9.3e-5));
        universe.add_component(Component::dark_energy(0.6889));
        universe
    }

    /// Create an Einstein-de Sitter universe (Ω_m = 1, flat)
    pub fn einstein_de_sitter(h0: f64) -> Self {
        let mut universe = Universe::new(h0);
        universe.add_component(Component::matter(1.0));
        universe
    }

    /// Total Omega (should be 1 for flat universe)
    pub fn omega_total(&self) -> f64 {
        self.components.iter().map(|c| c.omega_0).sum()
    }

    /// Hubble parameter at scale factor a: H(a) / H_0
    pub fn hubble_normalized(&self, a: f64) -> f64 {
        let sum: f64 = self.components
            .iter()
            .map(|c| c.density_evolution(a))
            .sum();
        sum.sqrt()
    }

    /// Hubble parameter at scale factor a [km/s/Mpc]
    pub fn hubble(&self, a: f64) -> f64 {
        self.h0 * self.hubble_normalized(a)
    }

    /// Hubble parameter at redshift z [km/s/Mpc]
    pub fn hubble_z(&self, z: f64) -> f64 {
        let a = 1.0 / (1.0 + z);
        self.hubble(a)
    }

    /// Deceleration parameter q(a) = -ä·a/ȧ²
    pub fn deceleration(&self, a: f64) -> f64 {
        let mut sum_rho_1_w = 0.0;
        let mut sum_rho = 0.0;

        for component in &self.components {
            let rho = component.density_evolution(a);
            sum_rho += rho;
            sum_rho_1_w += rho * (1.0 + component.w);
        }

        0.5 * sum_rho_1_w / sum_rho
    }

    /// Age of universe at scale factor a [Gyr]
    pub fn age(&self, a: f64) -> f64 {
        // Integrate t = ∫ da / (a·H(a))
        let n_steps = 1000;
        let da = a / n_steps as f64;

        let mut age = 0.0;
        for i in 1..=n_steps {
            let a_i = i as f64 * da;
            let h = self.hubble(a_i); // km/s/Mpc
            let h_si = h * 1e3 / (PARSEC * 1e6); // Convert to SI
            let dt = da / (a_i * h_si); // dt = da/(a*H) in seconds
            age += dt;
        }

        age / GYR // Convert to Gyr
    }

    /// Age of universe today [Gyr]
    pub fn age_today(&self) -> f64 {
        self.age(1.0)
    }

    /// Critical density today [kg/m³]
    pub fn critical_density(&self) -> f64 {
        critical_density(self.h0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_benchmark_universe() {
        let universe = Universe::benchmark();
        let omega_tot = universe.omega_total();
        assert_relative_eq!(omega_tot, 1.0, epsilon = 1e-3);
    }

    #[test]
    fn test_einstein_de_sitter() {
        let universe = Universe::einstein_de_sitter(70.0);
        assert_relative_eq!(universe.omega_total(), 1.0, epsilon = 1e-10);

        // H(a) = H_0 * a^(-3/2) for matter-only
        assert_relative_eq!(universe.hubble_normalized(0.5), 2.828427, epsilon = 1e-5);
    }

    #[test]
    fn test_age_calculation() {
        let universe = Universe::benchmark();
        let age = universe.age_today();
        // Should be around 13.8 Gyr
        assert!(age > 13.0 && age < 14.5);
    }
}
