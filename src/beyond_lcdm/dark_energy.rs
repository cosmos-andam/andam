//! Dark energy models beyond cosmological constant

use crate::dynamics::Universe;

/// Dark energy equation of state
#[derive(Debug, Clone)]
pub enum DarkEnergyModel {
    /// Cosmological constant: w = -1
    Lambda,
    /// Constant w
    ConstantW(f64),
    /// CPL parametrization: w(a) = w_0 + w_a(1-a)
    CPL { w_0: f64, w_a: f64 },
    /// Early dark energy
    EarlyDE { w_0: f64, omega_e: f64 },
}

impl DarkEnergyModel {
    /// Equation of state at scale factor a
    pub fn w(&self, a: f64) -> f64 {
        match self {
            DarkEnergyModel::Lambda => -1.0,
            DarkEnergyModel::ConstantW(w) => *w,
            DarkEnergyModel::CPL { w_0, w_a } => {
                w_0 + w_a * (1.0 - a)
            },
            DarkEnergyModel::EarlyDE { w_0, .. } => *w_0,
        }
    }

    /// Dark energy density evolution
    pub fn rho_de(&self, a: f64, omega_de_0: f64) -> f64 {
        match self {
            DarkEnergyModel::Lambda => omega_de_0,
            DarkEnergyModel::ConstantW(w) => {
                omega_de_0 * a.powf(-3.0 * (1.0 + w))
            },
            DarkEnergyModel::CPL { w_0, w_a } => {
                // ρ_DE ∝ exp(-3∫(1+w(a'))da'/a')
                let integral = (1.0 + w_0) * a.ln() + w_a * (a - 1.0);
                omega_de_0 * (-3.0 * integral).exp()
            },
            DarkEnergyModel::EarlyDE { w_0, .. } => {
                // Simplified early DE
                omega_de_0 * a.powf(-3.0 * (1.0 + w_0))
            },
        }
    }

    /// Modified Hubble parameter
    pub fn hubble_modified(&self, a: f64, _universe: &Universe) -> f64 {
        let omega_m = 0.3; // Get from universe
        let omega_de_0 = 0.7;

        let rho_m = omega_m * a.powf(-3.0);
        let rho_de = self.rho_de(a, omega_de_0);

        (rho_m + rho_de).sqrt()
    }
}

/// Compare models
pub fn compare_dark_energy_models(
    models: Vec<(&str, DarkEnergyModel)>,
    universe: &Universe,
) -> Vec<(String, Vec<(f64, f64)>)> {
    let mut results = Vec::new();

    for (name, model) in models {
        let mut data = Vec::new();

        for i in 0..100 {
            let a = 0.1 + (i as f64) * 0.009;
            let h = model.hubble_modified(a, universe);
            data.push((a, h));
        }

        results.push((name.to_string(), data));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambda_cdm() {
        let model = DarkEnergyModel::Lambda;
        assert_eq!(model.w(0.5), -1.0);
        assert_eq!(model.rho_de(0.5, 0.7), 0.7);
    }

    #[test]
    fn test_constant_w() {
        let model = DarkEnergyModel::ConstantW(-0.9);
        assert_eq!(model.w(0.5), -0.9);

        // Check density evolution
        let rho = model.rho_de(0.5, 1.0);
        assert!(rho > 0.0);
    }

    #[test]
    fn test_cpl() {
        let model = DarkEnergyModel::CPL { w_0: -0.9, w_a: -0.1 };
        let w_at_a = model.w(0.5);
        assert!((w_at_a - (-0.9 - 0.1 * 0.5)).abs() < 1e-10);
    }
}
