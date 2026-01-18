//! Example: Compare different distance measures

use andam::dynamics::Universe;
use andam::observations::{angular_diameter_distance, comoving_distance, luminosity_distance};
use andam::visualization::plotly_plots::create_interactive_plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let universe = Universe::benchmark();

    let mut data_lum = Vec::new();
    let mut data_ang = Vec::new();
    let mut data_com = Vec::new();

    for i in 1..100 {
        let z = (i as f64) * 0.05; // z from 0.05 to 5.0
        let d_l = luminosity_distance(z, &universe);
        let d_a = angular_diameter_distance(z, &universe);
        let d_c = comoving_distance(z, &universe);

        data_lum.push((z, d_l));
        data_ang.push((z, d_a));
        data_com.push((z, d_c));
    }

    let datasets = vec![
        (data_lum.as_slice(), "Luminosity Distance"),
        (data_ang.as_slice(), "Angular Diameter Distance"),
        (data_com.as_slice(), "Comoving Distance"),
    ];

    create_interactive_plot(
        "distance_measures.html",
        datasets,
        "Cosmological Distance Measures",
        "Redshift z",
        "Distance [Mpc]",
    )?;

    println!("Created distance_measures.html");
    Ok(())
}
