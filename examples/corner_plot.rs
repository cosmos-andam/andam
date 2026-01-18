//! Corner plot for MCMC chains
//!
//! This example demonstrates:
//! - MCMC sampling with multiple parameters
//! - Parameter posterior distributions
//! - 1D histograms and 2D scatter plots
//! - Confidence contours

use andam::statistics::mcmc::*;
use plotters::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== MCMC Parameter Estimation ===\n");

    // Example: fit parameters omega_m and sigma_8
    let params = vec![
        Parameter {
            name: "Omega_m".to_string(),
            initial: 0.3,
            min: 0.2,
            max: 0.4,
            proposal_width: 0.01,
        },
        Parameter {
            name: "sigma_8".to_string(),
            initial: 0.8,
            min: 0.6,
            max: 1.0,
            proposal_width: 0.02,
        },
    ];

    // Mock likelihood (Gaussian around truth)
    let log_likelihood = |theta: &[f64]| {
        let omega_m = theta[0];
        let sigma_8 = theta[1];

        let truth = [0.315, 0.81];
        let sigma = [0.01, 0.02];

        let chi2 = ((omega_m - truth[0]) / sigma[0]).powi(2)
                 + ((sigma_8 - truth[1]) / sigma[1]).powi(2);

        -0.5 * chi2
    };

    println!("Running MCMC sampler...");
    println!("  Walkers: 50");
    println!("  Steps: 1000");
    println!("  Burn-in: 200\n");

    let sampler = MCMCSampler::new(params, log_likelihood, 50, 1000);
    let chain = sampler.run(200);

    println!("\nMCMC complete!");
    println!("  Total samples: {}\n", chain.samples.len());

    // Create corner plot
    let filename = "corner_plot.png";
    let root = BitMapBackend::new(filename, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((2, 2));

    // 1D histogram for Omega_m (top-left)
    {
        let name = &chain.parameter_names[0];
        let samples = chain.get_parameter(name).unwrap();

        let min_val = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let max_val = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        let mut chart = ChartBuilder::on(&areas[0])
            .caption(name, ("sans-serif", 40))
            .margin(15)
            .x_label_area_size(50)
            .y_label_area_size(50)
            .build_cartesian_2d(min_val..max_val, 0u32..200u32)?;

        chart.configure_mesh()
            .x_desc(name)
            .y_desc("Count")
            .draw()?;

        // Create histogram
        let n_bins = 50;
        let bin_width = (max_val - min_val) / n_bins as f64;
        let mut bins = vec![0u32; n_bins];

        for &sample in &samples {
            let bin_idx = ((sample - min_val) / bin_width) as usize;
            if bin_idx < n_bins {
                bins[bin_idx] += 1;
            }
        }

        let histogram_data: Vec<_> = bins.iter().enumerate()
            .map(|(i, &count)| {
                let x = min_val + (i as f64 + 0.5) * bin_width;
                (x, count)
            })
            .collect();

        chart.draw_series(
            histogram_data.iter().map(|&(x, y)| {
                Rectangle::new(
                    [(x - bin_width / 2.0, 0), (x + bin_width / 2.0, y)],
                    BLUE.mix(0.7).filled(),
                )
            })
        )?;

        // Add mean line
        let mean = chain.mean(name).unwrap();
        chart.draw_series(LineSeries::new(
            vec![(mean, 0), (mean, 200)],
            &RED,
        ))?;
    }

    // Empty space (top-right)
    {
        let area = &areas[1];
        area.fill(&WHITE)?;
    }

    // 2D scatter plot (bottom-left)
    {
        let samples_x = chain.get_parameter(&chain.parameter_names[0]).unwrap();
        let samples_y = chain.get_parameter(&chain.parameter_names[1]).unwrap();

        let x_min = samples_x.iter().copied().fold(f64::INFINITY, f64::min);
        let x_max = samples_x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let y_min = samples_y.iter().copied().fold(f64::INFINITY, f64::min);
        let y_max = samples_y.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        let mut chart = ChartBuilder::on(&areas[2])
            .margin(15)
            .x_label_area_size(50)
            .y_label_area_size(50)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

        chart.configure_mesh()
            .x_desc(&chain.parameter_names[0])
            .y_desc(&chain.parameter_names[1])
            .draw()?;

        // Draw scatter points
        chart.draw_series(
            samples_x.iter().zip(samples_y.iter())
                .map(|(&x, &y)| Circle::new((x, y), 2, BLUE.mix(0.3).filled()))
        )?;

        // Add confidence ellipses
        let mean_x = chain.mean(&chain.parameter_names[0]).unwrap();
        let mean_y = chain.mean(&chain.parameter_names[1]).unwrap();
        let std_x = chain.std(&chain.parameter_names[0]).unwrap();
        let std_y = chain.std(&chain.parameter_names[1]).unwrap();

        for &n_sigma in &[1.0, 2.0, 3.0] {
            let mut ellipse = Vec::new();
            for i in 0..100 {
                let theta = 2.0 * std::f64::consts::PI * (i as f64) / 100.0;
                let x = mean_x + n_sigma * std_x * theta.cos();
                let y = mean_y + n_sigma * std_y * theta.sin();
                ellipse.push((x, y));
            }

            chart.draw_series(LineSeries::new(
                ellipse,
                &RED,
            ))?;
        }

        // Add mean point
        chart.draw_series(std::iter::once(Circle::new(
            (mean_x, mean_y),
            5,
            RED.filled(),
        )))?;
    }

    // 1D histogram for sigma_8 (bottom-right)
    {
        let name = &chain.parameter_names[1];
        let samples = chain.get_parameter(name).unwrap();

        let min_val = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let max_val = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        let mut chart = ChartBuilder::on(&areas[3])
            .caption(name, ("sans-serif", 40))
            .margin(15)
            .x_label_area_size(50)
            .y_label_area_size(50)
            .build_cartesian_2d(min_val..max_val, 0u32..200u32)?;

        chart.configure_mesh()
            .x_desc(name)
            .y_desc("Count")
            .draw()?;

        // Create histogram
        let n_bins = 50;
        let bin_width = (max_val - min_val) / n_bins as f64;
        let mut bins = vec![0u32; n_bins];

        for &sample in &samples {
            let bin_idx = ((sample - min_val) / bin_width) as usize;
            if bin_idx < n_bins {
                bins[bin_idx] += 1;
            }
        }

        let histogram_data: Vec<_> = bins.iter().enumerate()
            .map(|(i, &count)| {
                let x = min_val + (i as f64 + 0.5) * bin_width;
                (x, count)
            })
            .collect();

        chart.draw_series(
            histogram_data.iter().map(|&(x, y)| {
                Rectangle::new(
                    [(x - bin_width / 2.0, 0), (x + bin_width / 2.0, y)],
                    BLUE.mix(0.7).filled(),
                )
            })
        )?;

        // Add mean line
        let mean = chain.mean(name).unwrap();
        chart.draw_series(LineSeries::new(
            vec![(mean, 0), (mean, 200)],
            &RED,
        ))?;
    }

    root.present()?;
    println!("Created: {}", filename);

    // Print results
    println!("\n=== MCMC Results ===");
    for name in &chain.parameter_names {
        let mean = chain.mean(name).unwrap();
        let std = chain.std(name).unwrap();
        let p16 = chain.percentile(name, 16.0).unwrap();
        let p84 = chain.percentile(name, 84.0).unwrap();

        println!("{}: {:.4} ± {:.4} [{:.4}, {:.4}]",
                 name, mean, std, p16, p84);
    }

    Ok(())
}
