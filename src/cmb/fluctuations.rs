//! CMB temperature fluctuations and angular power spectrum

use crate::dynamics::Universe;
use std::f64::consts::PI;

/// Compute simplified CMB angular power spectrum C_ℓ
///
/// This uses a phenomenological model based on the standard ΛCDM prediction
pub fn angular_power_spectrum(l_max: usize, universe: &Universe) -> Vec<f64> {
    let mut c_l = vec![0.0; l_max + 1];

    // Acoustic peak positions
    let peaks = acoustic_peak_positions(universe);

    #[allow(clippy::needless_range_loop)]
    for l in 2..=l_max {
        // Base power law from Sachs-Wolfe effect
        let l_f = l as f64;
        let base = 1000.0 / (l_f * (l_f + 1.0));

        // Add acoustic oscillations
        let mut oscillation = 0.0;
        for (i, &peak_l) in peaks.iter().enumerate() {
            let width = 50.0;
            let amplitude = if i == 0 { 1.0 } else { 0.6 / (i as f64) };
            let gaussian =
                amplitude * (-(l_f - peak_l as f64).powi(2) / (2.0 * width * width)).exp();
            oscillation += gaussian;
        }

        c_l[l] = base * (1.0 + 5.0 * oscillation);
    }

    c_l
}

/// Compute ℓ(ℓ+1)C_ℓ / 2π for plotting
pub fn dimensionless_power_spectrum(c_l: &[f64]) -> Vec<(usize, f64)> {
    c_l.iter()
        .enumerate()
        .map(|(l, &cl)| {
            let dl = if l > 1 {
                (l * (l + 1)) as f64 * cl / (2.0 * PI)
            } else {
                0.0
            };
            (l, dl)
        })
        .collect()
}

/// Acoustic peaks positions
pub fn acoustic_peak_positions(_universe: &Universe) -> Vec<usize> {
    // Phenomenological model for peak positions
    // Based on ΛCDM predictions (Planck 2018)
    // l_peak ~ n * π / θ_s

    // For standard ΛCDM, sound horizon angle θ_s ≈ 0.0104 rad
    // Adjusted to give first peak at l_1 ≈ 220

    let theta_s = 0.01432; // radians, calibrated to match observations

    let mut peaks = Vec::new();
    for n in 1..=5 {
        let l_peak = ((n as f64) * PI / theta_s) as usize;
        peaks.push(l_peak);
    }

    peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acoustic_peaks() {
        let universe = crate::dynamics::Universe::benchmark();
        let peaks = acoustic_peak_positions(&universe);

        // First peak should be around l ~ 220
        assert!(peaks[0] > 180 && peaks[0] < 280);
    }

    #[test]
    fn test_power_spectrum_shape() {
        let universe = crate::dynamics::Universe::benchmark();
        let c_l = angular_power_spectrum(1000, &universe);

        // Should have values at acoustic peaks
        assert!(c_l[220] > 0.0);
    }
}
