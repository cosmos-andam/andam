//! Example: Plot scale factor evolution for different universe models

use andam::dynamics::{Universe, Component};
use andam::visualization::plots_2d::{create_multiline_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create different universe models
    let benchmark = Universe::benchmark();
    let einstein_de_sitter = Universe::einstein_de_sitter(70.0);

    let mut lambda_only = Universe::new(70.0);
    lambda_only.add_component(Component::dark_energy(1.0));

    // Generate data
    let n_points = 200;
    let mut data_benchmark = Vec::new();
    let mut data_eds = Vec::new();
    let mut data_lambda = Vec::new();

    for i in 0..n_points {
        let a = 0.01 + (i as f64 / n_points as f64) * 0.99; // a from 0.01 to 1.0
        let t_bench = benchmark.age(a);
        let t_eds = einstein_de_sitter.age(a);
        let t_lambda = lambda_only.age(a);

        data_benchmark.push((t_bench, a));
        data_eds.push((t_eds, a));
        data_lambda.push((t_lambda, a));
    }

    let datasets = [
        (data_benchmark.as_slice(), "Benchmark ΛCDM"),
        (data_eds.as_slice(), "Einstein-de Sitter"),
        (data_lambda.as_slice(), "Λ-only"),
    ];

    let config = PlotConfig {
        title: "Universe Evolution: Scale Factor vs Time".to_string(),
        x_label: "Time [Gyr]".to_string(),
        y_label: "Scale Factor a(t)".to_string(),
        ..Default::default()
    };

    create_multiline_plot("universe_evolution.png", &datasets, &config)?;
    println!("Created universe_evolution.png");

    Ok(())
}
