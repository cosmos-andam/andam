//! Cosmic web density field visualization
//!
//! This example demonstrates:
//! - Generating 3D Gaussian random density fields
//! - Extracting 2D slices for visualization
//! - Identifying overdense regions (halos)
//! - Creating density maps

use andam::structure::DensityField;
use plotters::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Cosmic Web Density Field ===\n");

    // Field parameters
    let n_grid = 64; // Grid resolution (64^3 cells)
    let box_size = 100.0; // Mpc/h
    let omega_m = 0.3;
    let z = 0.0;

    println!("Field parameters:");
    println!("  Grid size: {}^3 cells", n_grid);
    println!("  Box size: {} Mpc/h", box_size);
    println!("  Resolution: {:.2} Mpc/h per cell", box_size / n_grid as f64);
    println!("  Redshift z = {}\n", z);

    // Generate density field
    println!("Generating Gaussian random field...");
    let mut field = DensityField::new(n_grid, box_size);
    field.generate_gaussian(omega_m, z);
    println!("Done!\n");

    // Extract 2D slices
    let slice_index = n_grid / 2; // Middle slice
    println!("Extracting 2D slices at index {}...", slice_index);

    let xy_slice = field.slice(2, slice_index); // xy plane (z = const)
    let xz_slice = field.slice(1, slice_index); // xz plane (y = const)
    let yz_slice = field.slice(0, slice_index); // yz plane (x = const)

    println!("Done!\n");

    // Find overdense regions
    let threshold = 0.3; // Density contrast threshold
    println!("Extracting particles above threshold δ > {}...", threshold);
    let particles = field.extract_particles(threshold);
    println!("Found {} overdense cells ({:.1}% of volume)\n",
             particles.len(),
             100.0 * particles.len() as f64 / (n_grid * n_grid * n_grid) as f64);

    // Create visualizations
    println!("Creating density maps...");
    create_density_map(&xy_slice, "density_xy.png", "XY Plane")?;
    create_density_map(&xz_slice, "density_xz.png", "XZ Plane")?;
    create_density_map(&yz_slice, "density_yz.png", "YZ Plane")?;

    println!("\nGenerated plots:");
    println!("  - density_xy.png (z = {:.1} Mpc/h)", slice_index as f64 * box_size / n_grid as f64);
    println!("  - density_xz.png (y = {:.1} Mpc/h)", slice_index as f64 * box_size / n_grid as f64);
    println!("  - density_yz.png (x = {:.1} Mpc/h)", slice_index as f64 * box_size / n_grid as f64);
    println!("\nDone!");

    Ok(())
}

fn create_density_map(
    slice_data: &Vec<Vec<f64>>,
    filename: &str,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    let n_grid = slice_data.len();

    // Find min/max for color scale
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;

    for row in slice_data {
        for &val in row {
            min_val = min_val.min(val);
            max_val = max_val.max(val);
        }
    }

    let root = BitMapBackend::new(filename, (800, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Density Field: {}", title), ("sans-serif", 30).into_font())
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0..n_grid, 0..n_grid)?;

    chart.configure_mesh()
        .disable_mesh()
        .draw()?;

    // Draw density field as colored rectangles
    for i in 0..n_grid {
        for j in 0..n_grid {
            let val = slice_data[i][j];

            // Normalize to [0, 1]
            let normalized = if max_val > min_val {
                (val - min_val) / (max_val - min_val)
            } else {
                0.5
            };

            // Color mapping: blue (underdense) -> white (mean) -> red (overdense)
            let color = if normalized < 0.5 {
                // Blue to white
                let t = normalized * 2.0;
                RGBColor(
                    (t * 255.0) as u8,
                    (t * 255.0) as u8,
                    255,
                )
            } else {
                // White to red
                let t = (normalized - 0.5) * 2.0;
                RGBColor(
                    255,
                    ((1.0 - t) * 255.0) as u8,
                    ((1.0 - t) * 255.0) as u8,
                )
            };

            chart.draw_series(std::iter::once(Rectangle::new(
                [(i, j), (i + 1, j + 1)],
                color.filled(),
            )))?;
        }
    }

    // Add colorbar annotation (as a text label in the drawing area)
    let colorbar_text = format!("Delta: [{:.2}, {:.2}]", min_val, max_val);
    root.titled(&colorbar_text, ("sans-serif", 15))?;

    root.present()?;
    Ok(())
}
