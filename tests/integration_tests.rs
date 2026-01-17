//! Integration tests for andam

use andam::dynamics::Universe;
use approx::assert_relative_eq;

#[test]
fn test_benchmark_universe_age() {
    let universe = Universe::benchmark();
    let age = universe.age_today();

    // Planck 2018 gives ~13.8 Gyr
    assert!(age > 13.0 && age < 14.5, "Age should be around 13.8 Gyr, got {}", age);
}

#[test]
fn test_omega_total() {
    let universe = Universe::benchmark();
    let omega_tot = universe.omega_total();

    assert_relative_eq!(omega_tot, 1.0, epsilon = 1e-3);
}

#[test]
fn test_hubble_at_present() {
    let universe = Universe::benchmark();
    let h_normalized = universe.hubble_normalized(1.0);

    assert_relative_eq!(h_normalized, 1.0, epsilon = 1e-3);
}
