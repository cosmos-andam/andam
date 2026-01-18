//! Halo mass function and bias calculations
//!
//! This example demonstrates:
//! - Different halo mass function models (Press-Schechter, Sheth-Tormen, Tinker)
//! - Halo bias as a function of mass
//! - Number density of halos dn/dM
//! - Mass variance sigma(M)

use andam::structure::{HaloMassFunction, MassFunctionType};
use plotters::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Halo Mass Function and Bias ===\n");

    // Planck 2018 cosmological parameters
    let omega_m = 0.3;
    let omega_b = 0.05;
    let h = 0.7;
    let sigma_8 = 0.8;
    let n_s = 0.96;

    println!("Cosmological parameters:");
    println!("  Omega_m = {}", omega_m);
    println!("  Omega_b = {}", omega_b);
    println!("  h = {}", h);
    println!("  sigma_8 = {}", sigma_8);
    println!("  n_s = {}\n", n_s);

    // Create mass functions for different models
    let ps_mf = HaloMassFunction::new(
        omega_m,
        sigma_8,
        MassFunctionType::PressSchechter,
    );

    let st_mf = HaloMassFunction::new(
        omega_m,
        sigma_8,
        MassFunctionType::ShethTormen,
    );

    let tinker_mf = HaloMassFunction::new(
        omega_m,
        sigma_8,
        MassFunctionType::Tinker,
    );

    // Test at specific masses
    let z = 0.0;
    let test_masses = [1e11, 1e12, 1e13, 1e14, 1e15]; // M_sun

    println!("Halo abundance dn/dM at z = {}:", z);
    println!("  Mass [M_sun]        Press-Schechter    Sheth-Tormen      Tinker");
    for &mass in &test_masses {
        let dn_ps = ps_mf.dn_dm(mass, z);
        let dn_st = st_mf.dn_dm(mass, z);
        let dn_tinker = tinker_mf.dn_dm(mass, z);
        println!("  {:.2e}        {:.4e}      {:.4e}     {:.4e}",
                 mass, dn_ps, dn_st, dn_tinker);
    }
    println!();

    println!("Halo bias b(M) at z = {}:", z);
    println!("  Mass [M_sun]        Press-Schechter    Sheth-Tormen      Tinker");
    for &mass in &test_masses {
        let bias_ps = ps_mf.halo_bias(mass, z);
        let bias_st = st_mf.halo_bias(mass, z);
        let bias_tinker = tinker_mf.halo_bias(mass, z);
        println!("  {:.2e}        {:.4}             {:.4}            {:.4}",
                 mass, bias_ps, bias_st, bias_tinker);
    }
    println!();

    // Generate mass values (log-spaced)
    let n_mass = 100;
    let mass_values: Vec<f64> = (0..n_mass)
        .map(|i| {
            let log_m = 10.0 + 6.0 * i as f64 / (n_mass - 1) as f64;
            10_f64.powf(log_m)
        })
        .collect();

    // Create plots
    create_mass_function_plot(&mass_values, &ps_mf, &st_mf, &tinker_mf)?;
    create_bias_plot(&mass_values, &ps_mf, &st_mf, &tinker_mf)?;
    create_sigma_mass_plot(&mass_values, &ps_mf)?;

    println!("Generated plots:");
    println!("  - halo_mass_function.png");
    println!("  - halo_bias.png");
    println!("  - sigma_mass.png");
    println!("\nDone!");

    Ok(())
}

fn create_mass_function_plot(
    mass_values: &[f64],
    ps_mf: &HaloMassFunction,
    st_mf: &HaloMassFunction,
    tinker_mf: &HaloMassFunction,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("halo_mass_function.png", (1200, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Halo Mass Function dn/dM at z=0",
                 ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(80)
        .build_cartesian_2d(
            (1e10_f64..1e16_f64).log_scale(),
            (1e-30_f64..1e-5_f64).log_scale(),
        )?;

    chart.configure_mesh()
        .x_desc("Halo Mass M [M_sun]")
        .y_desc("dn/dM [(M_sun)^-1 (Mpc/h)^-3]")
        .draw()?;

    let z = 0.0;

    // Press-Schechter
    let dn_ps: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, ps_mf.dn_dm(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(dn_ps, &RED))?
        .label("Press-Schechter (1974)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Sheth-Tormen
    let dn_st: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, st_mf.dn_dm(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(dn_st, &BLUE))?
        .label("Sheth-Tormen (1999)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Tinker
    let dn_tinker: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, tinker_mf.dn_dm(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(dn_tinker, &GREEN))?
        .label("Tinker et al. (2008)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}

fn create_bias_plot(
    mass_values: &[f64],
    ps_mf: &HaloMassFunction,
    st_mf: &HaloMassFunction,
    tinker_mf: &HaloMassFunction,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("halo_bias.png", (1200, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Halo Bias b(M) at z=0",
                 ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(
            (1e10_f64..1e16_f64).log_scale(),
            0.5_f64..20.0_f64,
        )?;

    chart.configure_mesh()
        .x_desc("Halo Mass M [M_sun]")
        .y_desc("Halo Bias b(M)")
        .draw()?;

    // Draw horizontal line at b=1
    chart.draw_series(LineSeries::new(
        vec![(1e10, 1.0), (1e16, 1.0)],
        &BLACK.mix(0.3),
    ))?;

    let z = 0.0;

    // Press-Schechter
    let bias_ps: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, ps_mf.halo_bias(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(bias_ps, &RED))?
        .label("Press-Schechter")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Sheth-Tormen
    let bias_st: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, st_mf.halo_bias(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(bias_st, &BLUE))?
        .label("Sheth-Tormen")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Tinker
    let bias_tinker: Vec<(f64, f64)> = mass_values
        .iter()
        .map(|&m| (m, tinker_mf.halo_bias(m, z)))
        .collect();

    chart.draw_series(LineSeries::new(bias_tinker, &GREEN))?
        .label("Tinker et al.")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}

fn create_sigma_mass_plot(
    mass_values: &[f64],
    ps_mf: &HaloMassFunction,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("sigma_mass.png", (1200, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("RMS Mass Fluctuation sigma(M)",
                 ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(
            (1e10_f64..1e16_f64).log_scale(),
            (0.01_f64..10.0_f64).log_scale(),
        )?;

    chart.configure_mesh()
        .x_desc("Halo Mass M [M_sun]")
        .y_desc("sigma(M)")
        .draw()?;

    let redshifts = [0.0, 0.5, 1.0, 2.0];
    let colors = [&RED, &BLUE, &GREEN, &MAGENTA];

    for (idx, &z) in redshifts.iter().enumerate() {
        let color = colors[idx % colors.len()];

        let sigma: Vec<(f64, f64)> = mass_values
            .iter()
            .map(|&m| (m, ps_mf.sigma_mass(m, z)))
            .collect();

        chart.draw_series(LineSeries::new(sigma, color.stroke_width(3)))?
            .label(format!("z = {:.1}", z))
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(3))
            });
    }

    // Draw horizontal line at sigma = delta_c = 1.686
    chart.draw_series(LineSeries::new(
        vec![(1e10, 1.686), (1e16, 1.686)],
        &BLACK.mix(0.5),
    ))?
    .label("Critical density (delta_c = 1.686)")
    .legend(|(x, y)| {
        PathElement::new(vec![(x, y), (x + 20, y)], &BLACK.mix(0.5))
    });

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}
