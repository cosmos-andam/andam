//! Example: Compute and plot CMB angular power spectrum

use andam::cmb::fluctuations::{angular_power_spectrum, dimensionless_power_spectrum};
use andam::dynamics::Universe;
use andam::visualization::plotly_plots::create_interactive_plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let universe = Universe::benchmark();

    println!("Computing CMB angular power spectrum...");

    let l_max = 1000;
    let c_l = angular_power_spectrum(l_max, &universe);
    let d_l = dimensionless_power_spectrum(&c_l);

    // Convert to plottable data
    let data: Vec<(f64, f64)> = d_l
        .iter()
        .filter(|(l, _)| *l >= 2)
        .map(|(l, dl)| (*l as f64, *dl))
        .collect();

    let datasets = vec![(data.as_slice(), "CMB Power Spectrum")];

    create_interactive_plot(
        "cmb_power_spectrum.html",
        datasets,
        "CMB Angular Power Spectrum",
        "Multipole ℓ",
        "ℓ(ℓ+1)C_ℓ/2π [μK²]",
    )?;

    println!("Created cmb_power_spectrum.html");
    Ok(())
}
