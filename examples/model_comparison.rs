//! Compare different cosmological models
//!
//! This example demonstrates:
//! - ΛCDM (cosmological constant)
//! - Constant w dark energy
//! - CPL parametrization w(a) = w_0 + w_a(1-a)
//! - Hubble parameter evolution H(a)
//! - Equation of state w(a)

use andam::beyond_lcdm::dark_energy::*;
use andam::dynamics::Universe;
use plotters::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Beyond ΛCDM: Dark Energy Models ===\n");

    let universe = Universe::benchmark();

    // Define models
    let models = vec![
        ("ΛCDM (w = -1)", DarkEnergyModel::Lambda),
        ("Constant w = -0.9", DarkEnergyModel::ConstantW(-0.9)),
        (
            "CPL (w_0=-0.9, w_a=-0.1)",
            DarkEnergyModel::CPL {
                w_0: -0.9,
                w_a: -0.1,
            },
        ),
        (
            "Early DE",
            DarkEnergyModel::EarlyDE {
                w_0: -0.8,
                omega_e: 0.05,
            },
        ),
    ];

    println!("Comparing {} models:", models.len());
    for (name, _) in &models {
        println!("  - {}", name);
    }
    println!();

    let results = compare_dark_energy_models(models.clone(), &universe);

    let filename = "model_comparison.png";
    let root = BitMapBackend::new(filename, (2400, 1800)).into_drawing_area();
    root.fill(&WHITE)?;

    // Split into 2 panels
    let areas = root.split_evenly((2, 1));

    // Panel 1: Hubble parameter H(a)/H_0
    {
        let mut chart = ChartBuilder::on(&areas[0])
            .caption("Hubble Parameter Evolution", ("sans-serif", 48))
            .margin(20)
            .x_label_area_size(80)
            .y_label_area_size(100)
            .build_cartesian_2d(0.1..1.0, 0.5..1.8)?;

        chart
            .configure_mesh()
            .x_desc("Scale Factor a")
            .y_desc("H(a)/H_0")
            .x_label_style(("sans-serif", 36))
            .y_label_style(("sans-serif", 36))
            .draw()?;

        let colors = [BLUE, RED, GREEN, MAGENTA];

        for ((name, data), &color) in results.iter().zip(colors.iter()) {
            chart
                .draw_series(LineSeries::new(
                    data.iter().map(|&(a, h)| (a, h)),
                    color.stroke_width(4),
                ))?
                .label(name)
                .legend(move |(x, y)| {
                    PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(4))
                });
        }

        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.9))
            .border_style(BLACK)
            .label_font(("sans-serif", 32))
            .draw()?;

        // Add reference line at a=1 (today)
        chart.draw_series(LineSeries::new(
            vec![(1.0, 0.5), (1.0, 1.8)],
            BLACK.stroke_width(2),
        ))?;

        // Add annotation
        chart.draw_series(std::iter::once(Text::new(
            "Today (a=1)",
            (1.02, 1.5),
            ("sans-serif", 28).into_font(),
        )))?;
    }

    // Panel 2: Equation of state w(a)
    {
        let mut chart = ChartBuilder::on(&areas[1])
            .caption("Dark Energy Equation of State", ("sans-serif", 48))
            .margin(20)
            .x_label_area_size(80)
            .y_label_area_size(100)
            .build_cartesian_2d(0.1..1.0, -1.1..-0.6)?;

        chart
            .configure_mesh()
            .x_desc("Scale Factor a")
            .y_desc("w(a)")
            .x_label_style(("sans-serif", 36))
            .y_label_style(("sans-serif", 36))
            .draw()?;

        let colors = [BLUE, RED, GREEN, MAGENTA];

        // Generate w(a) data for each model
        for ((name, model), &color) in models.iter().zip(colors.iter()) {
            let w_data: Vec<_> = (0..100)
                .map(|i| {
                    let a = 0.1 + (i as f64) * 0.009;
                    let w = model.w(a);
                    (a, w)
                })
                .collect();

            chart
                .draw_series(LineSeries::new(w_data, color.stroke_width(4)))?
                .label(*name)
                .legend(move |(x, y)| {
                    PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(4))
                });
        }

        // Reference line at w = -1 (cosmological constant)
        chart.draw_series(LineSeries::new(
            vec![(0.1, -1.0), (1.0, -1.0)],
            BLACK.stroke_width(2),
        ))?;

        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.9))
            .border_style(BLACK)
            .label_font(("sans-serif", 32))
            .draw()?;

        // Add annotations
        chart.draw_series(std::iter::once(Text::new(
            "w = -1 (ΛCDM)",
            (0.15, -0.98),
            ("sans-serif", 28).into_font(),
        )))?;

        chart.draw_series(std::iter::once(Text::new(
            "CPL: w(a) = w_0 + w_a(1-a)",
            (0.15, -0.65),
            ("sans-serif", 28).into_font(),
        )))?;
    }

    root.present()?;
    println!("Created: {}", filename);

    // Print model details
    println!("\n=== Model Details ===");
    for (name, model) in &models {
        println!("\n{}", name);
        println!("  w(a=0.5) = {:.3}", model.w(0.5));
        println!("  w(a=1.0) = {:.3}", model.w(1.0));
        println!(
            "  ρ_DE(a=0.5)/ρ_DE(a=1.0) = {:.3}",
            model.rho_de(0.5, 0.7) / model.rho_de(1.0, 0.7)
        );
    }

    println!("\n=== Key Insights ===");
    println!("- ΛCDM: w = -1 always (constant vacuum energy)");
    println!("- Constant w: Similar to ΛCDM but w ≠ -1");
    println!("- CPL: Time-varying w(a), allows for dynamical dark energy");
    println!("- Early DE: Non-negligible dark energy at high redshift");
    println!("\nCurrent observations constrain w_0 ≈ -1 ± 0.05");

    Ok(())
}
