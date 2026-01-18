//! High-resolution BBN evolution with equations

use andam::early_universe::*;
use andam::visualization::equation_plots::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Computing BBN evolution...");
    println!("Using analytical freeze-out approach\n");

    // BBN parameters
    let params = BBNParameters::default();

    // Run BBN simulation
    let result = run_bbn(&params);

    // Print results
    println!("=== BBN Results ===");
    println!("Initial temperature: {:.2e} K", params.temp_initial);
    println!("Baryon-to-photon ratio η: {:.2e}", params.eta);
    println!("\nFinal Abundances:");
    println!("  Y_p (⁴He): {:.6}", result.yp());
    println!("  D/H: {:.6e}", result.dh_ratio());
    println!("  ³He/H: {:.6e}", result.he3h_ratio());
    println!("  ⁷Li/H: {:.6e}", result.li7h_ratio());

    // Create high-resolution plot
    let config = PublicationConfig {
        width: 2400,
        height: 1600,
        dpi: 300,
        ..Default::default()
    };

    plot_abundance_evolution("bbn_evolution_hires.png", &result.evolution, &config)?;

    println!("\nCreated: bbn_evolution_hires.png");
    Ok(())
}
