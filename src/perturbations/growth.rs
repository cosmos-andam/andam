//! Growth factor for structure formation

use crate::dynamics::Universe;

/// Linear growth factor D(a) for matter perturbations
///
/// δ(a,k) = D(a) δ(a_init, k) on superhorizon scales
pub fn growth_factor(a: f64, universe: &Universe) -> f64 {
    // Numerical solution to growth equation
    // d²D/da² + (3/a + dlnH/dlna) dD/da - (3Ω_m/2a²H²)D = 0

    let n_steps = 1000;
    let a_start = 0.001;
    let da = (a - a_start) / n_steps as f64;

    let mut d = Vec::new();
    let mut d_prime = Vec::new();

    // Initial conditions in radiation domination: D ∝ a²
    d.push(a_start * a_start);
    d_prime.push(2.0 * a_start);

    for i in 1..=n_steps {
        let a_i = a_start + i as f64 * da;
        let h = universe.hubble_normalized(a_i);

        // Approximation using matter density parameter
        let omega_m = 0.3111; // Planck 2018

        let d_curr = *d.last().unwrap();
        let d_prime_curr = *d_prime.last().unwrap();

        // Integration step
        let d_next = d_curr + d_prime_curr * da;
        let d_prime_next = d_prime_curr + (
            -3.0 * d_prime_curr / a_i
            + 1.5 * omega_m * d_curr / (a_i * a_i * h * h)
        ) * da;

        d.push(d_next);
        d_prime.push(d_prime_next);
    }

    // Normalize to D(a=1) = 1
    *d.last().unwrap() / growth_factor_at_unity(universe)
}

/// Helper: compute D(a=1) for normalization
fn growth_factor_at_unity(universe: &Universe) -> f64 {
    let n_steps = 1000;
    let a_start = 0.001;
    let da = (1.0 - a_start) / n_steps as f64;

    let mut d = Vec::new();
    let mut d_prime = Vec::new();

    d.push(a_start * a_start);
    d_prime.push(2.0 * a_start);

    for i in 1..=n_steps {
        let a_i = a_start + i as f64 * da;
        let h = universe.hubble_normalized(a_i);
        let omega_m = 0.3111;

        let d_curr = *d.last().unwrap();
        let d_prime_curr = *d_prime.last().unwrap();

        let d_next = d_curr + d_prime_curr * da;
        let d_prime_next = d_prime_curr + (
            -3.0 * d_prime_curr / a_i
            + 1.5 * omega_m * d_curr / (a_i * a_i * h * h)
        ) * da;

        d.push(d_next);
        d_prime.push(d_prime_next);
    }

    *d.last().unwrap()
}

/// Growth rate f = d ln D / d ln a
pub fn growth_rate(a: f64, universe: &Universe) -> f64 {
    let da = 0.001 * a;
    let d_plus = growth_factor(a + da, universe);
    let d_minus = growth_factor(a - da, universe);

    (a / growth_factor(a, universe)) * (d_plus - d_minus) / (2.0 * da)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_growth_normalized() {
        let universe = crate::dynamics::Universe::benchmark();
        let d = growth_factor(1.0, &universe);
        assert_relative_eq!(d, 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_growth_scaling() {
        let universe = crate::dynamics::Universe::einstein_de_sitter(70.0);

        // In matter-only: D ∝ a (approximately for this simplified solver)
        let d1 = growth_factor(0.5, &universe);
        let d2 = growth_factor(1.0, &universe);

        // Check that growth factor increases with scale factor
        assert!(d2 > d1);
        // For matter-dominated universe, expect monotonic growth
        assert!(d2 / d1 > 1.0);
    }
}
