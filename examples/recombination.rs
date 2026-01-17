//! Example: Plot ionization fraction during recombination

use andam::dynamics::Universe;
use andam::cmb::recombination::{ionization_fraction, recombination_redshift};
use andam::visualization::plots_2d::{create_line_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let universe = Universe::benchmark();

    let mut data = Vec::new();

    // From z=2000 to z=500
    for i in 0..200 {
        let z = 2000.0 - (i as f64) * 7.5;
        let x_e = ionization_fraction(z, &universe);
        data.push((z, x_e));
    }

    let z_rec = recombination_redshift(&universe);
    println!("Recombination redshift: z = {:.1}", z_rec);

    let config = PlotConfig {
        title: "Ionization Fraction during Recombination".to_string(),
        x_label: "Redshift z".to_string(),
        y_label: "Ionization Fraction X_e".to_string(),
        ..Default::default()
    };

    create_line_plot("recombination.png", &data, &config)?;
    println!("Created recombination.png");

    Ok(())
}
