//! Example: Create a Hubble diagram showing velocity vs distance

use andam::dynamics::Universe;
use andam::visualization::plots_2d::{create_line_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let universe = Universe::benchmark();

    // Generate distance and velocity data
    let mut data = Vec::new();
    for i in 0..100 {
        let distance_mpc = i as f64 * 10.0; // 0 to 1000 Mpc
        let velocity = universe.h0 * distance_mpc; // v = H_0 * d
        data.push((distance_mpc, velocity));
    }

    let config = PlotConfig {
        title: "Hubble Diagram".to_string(),
        x_label: "Distance [Mpc]".to_string(),
        y_label: "Velocity [km/s]".to_string(),
        ..Default::default()
    };

    create_line_plot("hubble_diagram.png", &data, &config)?;
    println!("Created hubble_diagram.png");

    Ok(())
}
