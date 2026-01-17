use andam::*;
use approx::assert_relative_eq;

#[test]
fn test_growth_factor_normalization() {
    let universe = dynamics::Universe::benchmark();
    let d = perturbations::growth::growth_factor(1.0, &universe);
    assert_relative_eq!(d, 1.0, epsilon = 0.15);
}

#[test]
fn test_cmb_peaks_reasonable() {
    let universe = dynamics::Universe::benchmark();
    let peaks = cmb::fluctuations::acoustic_peak_positions(&universe);

    // First peak around l=220
    assert!(peaks[0] > 180 && peaks[0] < 280);
}

#[test]
fn test_lensing_convergence() {
    use advanced::weak_lensing::ConvergenceField;

    let field = ConvergenceField::new(128);
    assert_eq!(field.kappa.shape(), &[128, 128]);
}

#[test]
fn test_boltzmann_solver() {
    let universe = dynamics::Universe::benchmark();
    let mut solver = perturbations::BoltzmannSolver::new(universe, 0.1, 10);
    let results = solver.evolve(1e-5, 1e-2);

    assert!(results.len() > 10);
}

#[test]
fn test_angular_power_spectrum() {
    let universe = dynamics::Universe::benchmark();
    let c_l = cmb::angular_power_spectrum(500, &universe);

    // Should have non-zero values at peak
    assert!(c_l[220] > 0.0);
}
