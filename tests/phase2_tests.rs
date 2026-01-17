use andam::*;
use approx::assert_relative_eq;

#[test]
fn test_distance_consistency() {
    let universe = Universe::benchmark();
    let z = 1.0;

    let d_l = luminosity_distance(z, &universe);
    let d_a = angular_diameter_distance(z, &universe);

    // Etherington reciprocity relation
    assert_relative_eq!(d_l / d_a, (1.0 + z).powi(2), epsilon = 1e-6);
}

#[test]
fn test_recombination_physics() {
    let universe = Universe::benchmark();
    let z_rec = recombination_redshift(&universe);

    let x_e_rec = ionization_fraction(z_rec, &universe);

    // Should be around 0.5 at recombination
    assert!((x_e_rec - 0.5).abs() < 0.1);
}

#[test]
fn test_power_spectrum_shape() {
    use structure::power_spectrum::dimensionless_power;

    let omega_m = 0.3;
    let omega_b = 0.05;
    let h = 0.7;
    let a_s = 2.1e-9;
    let n_s = 0.96;

    // Power should have reasonable values
    let p_small = dimensionless_power(0.001, 0.0, omega_m, omega_b, h, a_s, n_s);
    let p_mid = dimensionless_power(0.01, 0.0, omega_m, omega_b, h, a_s, n_s);
    let p_large = dimensionless_power(0.1, 0.0, omega_m, omega_b, h, a_s, n_s);

    // All powers should be positive
    assert!(p_small > 0.0);
    assert!(p_mid > 0.0);
    assert!(p_large > 0.0);

    // Power spectrum should suppress at small scales due to transfer function
    assert!(p_large < p_mid);
}
