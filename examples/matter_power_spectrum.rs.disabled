//! Example: Plot matter power spectrum

use andam::structure::power_spectrum::dimensionless_power;
use andam::visualization::plotly_plots::create_loglog_interactive;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Planck 2018 parameters
    let omega_m = 0.3111;
    let omega_b = 0.0490;
    let h = 0.6766;
    let a_s = 2.1e-9;
    let n_s = 0.9665;

    let mut data_z0 = Vec::new();
    let mut data_z1 = Vec::new();
    let mut data_z10 = Vec::new();

    // k from 10^-4 to 1 h/Mpc
    for i in 0..200 {
        let log_k = -4.0 + (i as f64) * 0.025;
        let k = 10_f64.powf(log_k);

        let delta2_z0 = dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s);
        let delta2_z1 = dimensionless_power(k, 1.0, omega_m, omega_b, h, a_s, n_s);
        let delta2_z10 = dimensionless_power(k, 10.0, omega_m, omega_b, h, a_s, n_s);

        data_z0.push((k, delta2_z0));
        data_z1.push((k, delta2_z1));
        data_z10.push((k, delta2_z10));
    }

    let datasets = vec![
        (data_z0.as_slice(), "z = 0"),
        (data_z1.as_slice(), "z = 1"),
        (data_z10.as_slice(), "z = 10"),
    ];

    create_loglog_interactive(
        "matter_power_spectrum.html",
        datasets,
        "Matter Power Spectrum",
        "Wavenumber k [h/Mpc]",
        "Δ²(k)",
    )?;

    println!("Created matter_power_spectrum.html");
    Ok(())
}
