//! Universe components and their properties

/// Types of components in the universe
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ComponentType {
    Matter,
    Radiation,
    DarkEnergy,
    Curvature,
}

/// Represents a component of the universe with equation of state
#[derive(Debug, Clone)]
pub struct Component {
    pub name: String,
    pub component_type: ComponentType,
    pub omega_0: f64, // Density parameter today
    pub w: f64,       // Equation of state parameter
}

impl Component {
    /// Create a matter component (w = 0)
    pub fn matter(omega_0: f64) -> Self {
        Component {
            name: "Matter".to_string(),
            component_type: ComponentType::Matter,
            omega_0,
            w: 0.0,
        }
    }

    /// Create a radiation component (w = 1/3)
    pub fn radiation(omega_0: f64) -> Self {
        Component {
            name: "Radiation".to_string(),
            component_type: ComponentType::Radiation,
            omega_0,
            w: 1.0 / 3.0,
        }
    }

    /// Create a cosmological constant (w = -1)
    pub fn dark_energy(omega_0: f64) -> Self {
        Component {
            name: "Dark Energy".to_string(),
            component_type: ComponentType::DarkEnergy,
            omega_0,
            w: -1.0,
        }
    }

    /// Create a curvature component (w = -1/3)
    pub fn curvature(omega_0: f64) -> Self {
        Component {
            name: "Curvature".to_string(),
            component_type: ComponentType::Curvature,
            omega_0,
            w: -1.0 / 3.0,
        }
    }

    /// Energy density at scale factor a, normalized to today
    /// ρ(a) / ρ_0 = a^(-3(1+w))
    pub fn density_evolution(&self, a: f64) -> f64 {
        self.omega_0 * a.powf(-3.0 * (1.0 + self.w))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_matter_scaling() {
        let matter = Component::matter(0.3);
        assert_relative_eq!(matter.density_evolution(0.5), 0.3 * 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_radiation_scaling() {
        let radiation = Component::radiation(1e-4);
        assert_relative_eq!(
            radiation.density_evolution(0.5),
            1e-4 * 16.0,
            epsilon = 1e-14
        );
    }
}
