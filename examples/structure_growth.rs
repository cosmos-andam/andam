//! Example: Plot growth factor evolution

use andam::dynamics::Universe;
use andam::perturbations::growth::growth_factor;
use andam::visualization::plots_2d::{create_multiline_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let benchmark = Universe::benchmark();
    let eds = Universe::einstein_de_sitter(70.0);

    println!("Computing growth factor evolution...");

    let mut data_benchmark = Vec::new();
    let mut data_eds = Vec::new();

    for i in 1..200 {
        let a = 0.01 + (i as f64) * 0.99 / 200.0;

        let d_bench = growth_factor(a, &benchmark);
        let d_eds = growth_factor(a, &eds);

        data_benchmark.push((a, d_bench));
        data_eds.push((a, d_eds));
    }

    let datasets = [
        (data_benchmark.as_slice(), "Benchmark ΛCDM"),
        (data_eds.as_slice(), "Einstein-de Sitter"),
    ];

    let config = PlotConfig {
        title: "Growth Factor Evolution".to_string(),
        x_label: "Scale Factor a".to_string(),
        y_label: "Growth Factor D(a)".to_string(),
        ..Default::default()
    };

    create_multiline_plot("growth_factor.png", &datasets, &config)?;
    println!("Created growth_factor.png");

    Ok(())
}
