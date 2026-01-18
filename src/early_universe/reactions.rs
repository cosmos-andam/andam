//! Nuclear reactions in the early universe

use crate::constants::*;

/// Nuclear species
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuclide {
    Neutron,
    Proton,
    Deuterium,
    Tritium,
    Helium3,
    Helium4,
    Lithium7,
    Beryllium7,
}

impl Nuclide {
    /// Atomic mass number
    pub fn mass_number(&self) -> u8 {
        match self {
            Nuclide::Neutron => 1,
            Nuclide::Proton => 1,
            Nuclide::Deuterium => 2,
            Nuclide::Tritium => 3,
            Nuclide::Helium3 => 3,
            Nuclide::Helium4 => 4,
            Nuclide::Lithium7 => 7,
            Nuclide::Beryllium7 => 7,
        }
    }

    /// Binding energy \[MeV\]
    pub fn binding_energy(&self) -> f64 {
        match self {
            Nuclide::Neutron => 0.0,
            Nuclide::Proton => 0.0,
            Nuclide::Deuterium => 2.224,
            Nuclide::Tritium => 8.482,
            Nuclide::Helium3 => 7.718,
            Nuclide::Helium4 => 28.296,
            Nuclide::Lithium7 => 39.245,
            Nuclide::Beryllium7 => 37.600,
        }
    }

    /// Symbol for plotting
    pub fn symbol(&self) -> &'static str {
        match self {
            Nuclide::Neutron => "n",
            Nuclide::Proton => "p",
            Nuclide::Deuterium => "²H",
            Nuclide::Tritium => "³H",
            Nuclide::Helium3 => "³He",
            Nuclide::Helium4 => "⁴He",
            Nuclide::Lithium7 => "⁷Li",
            Nuclide::Beryllium7 => "⁷Be",
        }
    }

    /// LaTeX notation
    pub fn latex(&self) -> &'static str {
        match self {
            Nuclide::Neutron => "n",
            Nuclide::Proton => "p",
            Nuclide::Deuterium => r"^{2}\mathrm{H}",
            Nuclide::Tritium => r"^{3}\mathrm{H}",
            Nuclide::Helium3 => r"^{3}\mathrm{He}",
            Nuclide::Helium4 => r"^{4}\mathrm{He}",
            Nuclide::Lithium7 => r"^{7}\mathrm{Li}",
            Nuclide::Beryllium7 => r"^{7}\mathrm{Be}",
        }
    }
}

/// Nuclear reaction
#[derive(Debug, Clone)]
pub struct Reaction {
    pub reactants: Vec<Nuclide>,
    pub products: Vec<Nuclide>,
    pub q_value: f64, // Energy release \[MeV\]
    pub label: String,
}

impl Reaction {
    /// Create new reaction
    pub fn new(reactants: Vec<Nuclide>, products: Vec<Nuclide>, label: &str) -> Self {
        // Calculate Q-value from binding energies
        let be_products: f64 = products.iter().map(|n| n.binding_energy()).sum();
        let be_reactants: f64 = reactants.iter().map(|n| n.binding_energy()).sum();
        let q_value = be_products - be_reactants;

        Reaction {
            reactants,
            products,
            q_value,
            label: label.to_string(),
        }
    }

    /// Reaction rate coefficient [cm³/s] at temperature T \[K\]
    /// These rates are from Coc et al. and similar BBN codes
    /// Units: cm³/s (for use with number densities in cm⁻³)
    pub fn rate_coefficient(&self, temp_k: f64) -> f64 {
        let t9 = temp_k / 1e9; // Temperature in GK

        // Scaling factor to convert from typical parameterizations
        // Literature rates are often in units that need scaling
        let scale = 1e-20; // Scale down by factor of 1e20 for stability

        // Use fitting formulas from literature (Coc et al., Cyburt et al.)
        let base_rate = match self.label.as_str() {
            "p(n,γ)d" => {
                // n + p → d + γ (dominant at early times)
                // This is the key bottleneck reaction
                4.742e-4 * (1.0 / t9)
            }
            "d(p,γ)³He" => {
                // d + p → ³He + γ
                2.24e3
                    * t9.powf(-2.0 / 3.0)
                    * (-3.720 / t9.powf(1.0 / 3.0)).exp()
                    * (1.0 + 0.112 * t9 + 3.38 * t9.powi(2))
            }
            "d(d,n)³He" => {
                // d + d → ³He + n
                3.95e8
                    * t9.powf(-2.0 / 3.0)
                    * (-4.259 / t9.powf(1.0 / 3.0)).exp()
                    * (1.0 + 0.0979 * t9)
            }
            "d(d,p)t" => {
                // d + d → t + p
                4.17e8
                    * t9.powf(-2.0 / 3.0)
                    * (-4.258 / t9.powf(1.0 / 3.0)).exp()
                    * (1.0 + 0.0978 * t9)
            }
            "t(d,n)⁴He" => {
                // t + d → ⁴He + n (produces helium-4)
                1.63e11 * t9.powf(-2.0 / 3.0) * (-4.871 / t9.powf(1.0 / 3.0)).exp()
            }
            "³He(d,p)⁴He" => {
                // ³He + d → ⁴He + p (produces helium-4)
                5.59e10
                    * t9.powf(-2.0 / 3.0)
                    * (-7.183 / t9.powf(1.0 / 3.0)).exp()
                    * (1.0 + 0.0488 * t9)
            }
            "³He(n,p)t" => {
                // ³He + n → t + p
                7.21e8 * (1.0 - 0.86 * t9 + 0.429 * t9.powi(2))
            }
            "⁴He(t,γ)⁷Li" => {
                // ⁴He + t → ⁷Li + γ
                2.39e4 * t9.powf(-2.0 / 3.0) * (-8.09 / t9.powf(1.0 / 3.0)).exp()
            }
            _ => 1e-40, // Default very small value
        };

        base_rate * scale
    }

    /// Reverse reaction rate using detailed balance
    pub fn reverse_rate(&self, temp_k: f64) -> f64 {
        // At these temperatures, reverse rates are negligible for bound nuclei
        // Simplified approach: only include reverse for deuterium photodissociation
        if self.label == "p(n,γ)d" {
            // d → p + n (photodissociation)
            // Suppressed by Boltzmann factor at low T
            let forward_rate = self.rate_coefficient(temp_k);
            let t_mev = K_B * temp_k * J_TO_EV * 1e-6;
            forward_rate * (-self.q_value / t_mev).exp()
        } else {
            // Other reverse reactions negligible
            0.0
        }
    }

    /// LaTeX representation of reaction
    pub fn latex_equation(&self) -> String {
        let reactants_str: Vec<_> = self.reactants.iter().map(|n| n.latex()).collect();
        let products_str: Vec<_> = self.products.iter().map(|n| n.latex()).collect();

        format!(
            "${}\\rightarrow {}$",
            reactants_str.join(" + "),
            products_str.join(" + ")
        )
    }
}

/// Standard BBN reaction network
pub fn standard_bbn_reactions() -> Vec<Reaction> {
    vec![
        Reaction::new(
            vec![Nuclide::Proton, Nuclide::Neutron],
            vec![Nuclide::Deuterium],
            "p(n,γ)d",
        ),
        Reaction::new(
            vec![Nuclide::Deuterium, Nuclide::Proton],
            vec![Nuclide::Helium3],
            "d(p,γ)³He",
        ),
        Reaction::new(
            vec![Nuclide::Deuterium, Nuclide::Deuterium],
            vec![Nuclide::Helium3, Nuclide::Neutron],
            "d(d,n)³He",
        ),
        Reaction::new(
            vec![Nuclide::Deuterium, Nuclide::Deuterium],
            vec![Nuclide::Tritium, Nuclide::Proton],
            "d(d,p)t",
        ),
        Reaction::new(
            vec![Nuclide::Tritium, Nuclide::Deuterium],
            vec![Nuclide::Helium4, Nuclide::Neutron],
            "t(d,n)⁴He",
        ),
        Reaction::new(
            vec![Nuclide::Helium3, Nuclide::Deuterium],
            vec![Nuclide::Helium4, Nuclide::Proton],
            "³He(d,p)⁴He",
        ),
        Reaction::new(
            vec![Nuclide::Helium3, Nuclide::Neutron],
            vec![Nuclide::Tritium, Nuclide::Proton],
            "³He(n,p)t",
        ),
        Reaction::new(
            vec![Nuclide::Helium4, Nuclide::Tritium],
            vec![Nuclide::Lithium7],
            "⁴He(t,γ)⁷Li",
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binding_energies() {
        assert!((Nuclide::Deuterium.binding_energy() - 2.224).abs() < 0.001);
        assert!((Nuclide::Helium4.binding_energy() - 28.296).abs() < 0.001);
    }

    #[test]
    fn test_reaction_rates() {
        let reaction = Reaction::new(
            vec![Nuclide::Proton, Nuclide::Neutron],
            vec![Nuclide::Deuterium],
            "p(n,γ)d",
        );

        let rate = reaction.rate_coefficient(1e9);
        assert!(rate > 0.0);
    }
}
