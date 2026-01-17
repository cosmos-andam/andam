//! Neutron-proton freeze-out calculations

use crate::constants::*;

/// Weak interaction rate [s^-1] at temperature T [K]
pub fn weak_interaction_rate(temp_k: f64) -> f64 {
    // Simplified formula for n <-> p + e + ν_e weak interaction rate
    let t_mev = K_B * temp_k * J_TO_EV * 1e-6; // Temperature in MeV

    // Rate proportional to T^5
    let rate = 1.13 * t_mev.powi(5);

    rate
}

/// Hubble parameter [s^-1] at temperature T [K] in radiation-dominated era
pub fn hubble_rate_radiation(temp_k: f64) -> f64 {
    // H = sqrt(8πGρ/3) where ρ is radiation density
    // ρ = (π^2/30) g_* k_B^4 T^4 / (ℏ^3 c^5)
    let g_star = 10.75; // Effective degrees of freedom (photons + e± + 3 neutrinos)

    let rho = (std::f64::consts::PI.powi(2) / 30.0) * g_star
        * K_B.powi(4) * temp_k.powi(4)
        / (HBAR.powi(3) * C.powi(5));

    let hubble = ((8.0 * std::f64::consts::PI * G * rho) / 3.0).sqrt();

    hubble
}

/// Freeze-out temperature [K] where weak rate equals Hubble rate
pub fn freezeout_temperature() -> f64 {
    // Solve Γ_weak(T) = H(T)
    // This occurs at T ~ 0.7 MeV
    let t_freeze_mev = 0.7; // MeV
    let t_freeze_k = t_freeze_mev * 1e6 / (K_B * J_TO_EV);

    t_freeze_k
}

/// Neutron-to-proton ratio at freeze-out
pub fn freezeout_np_ratio() -> f64 {
    let t_freeze = freezeout_temperature();

    // n/p = exp(-Δm/kT)
    let delta_m = (M_N - M_P) * C * C; // Mass difference [J]
    let delta_m_mev = delta_m * J_TO_EV * 1e-6; // [MeV]
    let t_mev = K_B * t_freeze * J_TO_EV * 1e-6;

    (-delta_m_mev / t_mev).exp()
}

/// Estimate final helium-4 mass fraction from freeze-out
pub fn estimate_yp_from_freezeout() -> f64 {
    let np_freeze = freezeout_np_ratio();

    // Almost all neutrons end up in He-4
    // Y_p ≈ 2(n/p) / (1 + n/p)
    2.0 * np_freeze / (1.0 + np_freeze)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_freezeout_temperature() {
        let t_freeze = freezeout_temperature();

        // Should be around 10^10 K
        assert!(t_freeze > 5e9 && t_freeze < 2e10);
    }

    #[test]
    fn test_freezeout_np_ratio() {
        let np = freezeout_np_ratio();

        // Should be around 1/6 to 1/7
        assert!(np > 0.1 && np < 0.2);
    }

    #[test]
    fn test_yp_estimate() {
        let yp = estimate_yp_from_freezeout();

        // Should be around 0.25
        assert!(yp > 0.20 && yp < 0.30);
    }
}
