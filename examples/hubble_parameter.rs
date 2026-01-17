//! Example: Plot Hubble parameter evolution

use andam::dynamics::Universe;
use andam::visualization::plots_2d::{create_loglog_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let universe = Universe::benchmark();

    // Generate H(a) data
    let mut data = Vec::new();
    let n_points = 200;

    for i in 1..=n_points {
        let a = (i as f64) / (n_points as f64);
        let h = universe.hubble_normalized(a);
        data.push((a, h));
    }

    let config = PlotConfig {
        title: "Hubble Parameter Evolution".to_string(),
        x_label: "Scale Factor a".to_string(),
        y_label: "H(a) / H₀".to_string(),
        ..Default::default()
    };

    create_loglog_plot("hubble_evolution.png", &data, &config)?;
    println!("Created hubble_evolution.png");

    Ok(())
}
