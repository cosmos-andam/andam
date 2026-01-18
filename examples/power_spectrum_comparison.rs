//! Compare linear and non-linear matter power spectra using HALOFIT
//!
//! This example demonstrates:
//! - Linear matter power spectrum P_lin(k)
//! - Non-linear power spectrum P_nl(k) from HALOFIT
//! - Boost factor B(k) = P_nl(k) / P_lin(k)
//! - Evolution with redshift

use andam::structure::{
    matter_power_spectrum,
    HalofitSpectrum,
};
use plotters::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Matter Power Spectrum: Linear vs Non-linear ===\n");

    // Planck 2018 cosmological parameters
    let omega_m = 0.3;
    let omega_b = 0.05;
    let h = 0.7;
    let sigma_8 = 0.8;
    let n_s = 0.96;

    // Create HALOFIT calculator
    let halofit = HalofitSpectrum::new(omega_m, omega_b, h, sigma_8, n_s);

    // Redshifts to plot
    let redshifts = [0.0, 0.5, 1.0, 2.0];

    println!("Cosmological parameters:");
    println!("  Omega_m = {}", omega_m);
    println!("  Omega_b = {}", omega_b);
    println!("  h = {}", h);
    println!("  sigma_8 = {}", sigma_8);
    println!("  n_s = {}\n", n_s);

    // Calculate k_nl at different redshifts
    println!("Non-linear scale k_nl:");
    for &z in &redshifts {
        let k_nl = halofit.k_nonlinear(z);
        println!("  z = {:.1}: k_nl = {:.4} h/Mpc", z, k_nl);
    }
    println!();

    // Generate k values (log-spaced from 0.001 to 10 h/Mpc)
    let n_k = 100;
    let k_values: Vec<f64> = (0..n_k)
        .map(|i| {
            let log_k = -3.0 + 4.0 * i as f64 / (n_k - 1) as f64;
            10_f64.powf(log_k)
        })
        .collect();

    // Create plots
    create_power_spectrum_plot(&k_values, &halofit, &redshifts, omega_m, omega_b, h, n_s)?;
    create_boost_factor_plot(&k_values, &halofit, &redshifts, omega_m, omega_b, h, n_s)?;

    println!("\nGenerated plots:");
    println!("  - power_spectrum_comparison.png");
    println!("  - boost_factor.png");
    println!("\nDone!");

    Ok(())
}

fn create_power_spectrum_plot(
    k_values: &[f64],
    halofit: &HalofitSpectrum,
    redshifts: &[f64],
    omega_m: f64,
    omega_b: f64,
    h: f64,
    n_s: f64,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("power_spectrum_comparison.png", (1200, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Matter Power Spectrum: Linear vs Non-linear (HALOFIT)",
                 ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(
            (0.001_f64..10.0_f64).log_scale(),
            (1e-2_f64..1e5_f64).log_scale(),
        )?;

    chart.configure_mesh()
        .x_desc("Wavenumber k [h/Mpc]")
        .y_desc("Power Spectrum P(k) [(Mpc/h)^3]")
        .draw()?;

    let colors = [&RED, &BLUE, &GREEN, &MAGENTA];

    // Plot for each redshift
    for (idx, &z) in redshifts.iter().enumerate() {
        let color = colors[idx % colors.len()];

        // Linear power spectrum
        let p_lin: Vec<(f64, f64)> = k_values
            .iter()
            .map(|&k| {
                let p = matter_power_spectrum(k, z, omega_m, omega_b, h, 2.1e-9, n_s);
                (k, p)
            })
            .collect();

        // Non-linear power spectrum
        let p_nl: Vec<(f64, f64)> = k_values
            .iter()
            .map(|&k| {
                let p = halofit.nonlinear_power(k, z);
                (k, p)
            })
            .collect();

        // Plot linear (dashed line)
        chart.draw_series(LineSeries::new(p_lin, color.stroke_width(2).stroke_width(2)))?
            .label(format!("Linear z={:.1}", z))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            });

        // Plot non-linear (solid line)
        chart.draw_series(LineSeries::new(p_nl, color.stroke_width(3)))?
            .label(format!("Non-linear z={:.1}", z))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(3))
            });
    }

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}

fn create_boost_factor_plot(
    k_values: &[f64],
    halofit: &HalofitSpectrum,
    redshifts: &[f64],
    omega_m: f64,
    omega_b: f64,
    h: f64,
    n_s: f64,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("boost_factor.png", (1200, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Non-linear Boost Factor B(k) = P_nl(k) / P_lin(k)",
                 ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(
            (0.001_f64..10.0_f64).log_scale(),
            0.8_f64..10.0_f64,
        )?;

    chart.configure_mesh()
        .x_desc("Wavenumber k [h/Mpc]")
        .y_desc("Boost Factor B(k)")
        .draw()?;

    // Draw horizontal line at B=1
    chart.draw_series(LineSeries::new(
        vec![(0.001, 1.0), (10.0, 1.0)],
        &BLACK.mix(0.3),
    ))?;

    let colors = [&RED, &BLUE, &GREEN, &MAGENTA];

    // Plot boost factor for each redshift
    for (idx, &z) in redshifts.iter().enumerate() {
        let color = colors[idx % colors.len()];

        let boost: Vec<(f64, f64)> = k_values
            .iter()
            .map(|&k| {
                let p_lin = matter_power_spectrum(k, z, omega_m, omega_b, h, 2.1e-9, n_s);
                let p_nl = halofit.nonlinear_power(k, z);
                let b = if p_lin > 0.0 { p_nl / p_lin } else { 1.0 };
                (k, b)
            })
            .collect();

        chart.draw_series(LineSeries::new(boost, color.stroke_width(3)))?
            .label(format!("z = {:.1}", z))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(3))
            });

        // Mark k_nl
        let k_nl = halofit.k_nonlinear(z);
        chart.draw_series(std::iter::once(Circle::new(
            (k_nl, 1.5),
            5,
            color.filled(),
        )))?;
    }

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}
