//! High-resolution plots with embedded equations

use plotters::prelude::*;

/// Configuration for publication-quality plots
pub struct PublicationConfig {
    pub width: u32,
    pub height: u32,
    pub dpi: u32,
    pub font_size_title: u32,
    pub font_size_label: u32,
    pub font_size_equation: u32,
    pub font_family: String,
}

impl Default for PublicationConfig {
    fn default() -> Self {
        PublicationConfig {
            width: 1920,   // High resolution
            height: 1080,
            dpi: 300,      // Publication quality
            font_size_title: 48,
            font_size_label: 36,
            font_size_equation: 32,
            font_family: "DejaVu Sans".to_string(),
        }
    }
}

/// Create plot with equation overlay
pub fn create_plot_with_equation(
    filename: &str,
    data_series: Vec<(&[(f64, f64)], &str, RGBColor)>,
    title: &str,
    x_label: &str,
    y_label: &str,
    equations: Vec<(String, (f64, f64))>, // (equation string, (x, y) position in data coords)
    config: &PublicationConfig,
) -> Result<(), Box<dyn std::error::Error>> {

    let root = BitMapBackend::new(filename, (config.width, config.height))
        .into_drawing_area();
    root.fill(&WHITE)?;

    // Find data ranges
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for (data, _, _) in &data_series {
        for &(x, y) in *data {
            x_min = x_min.min(x);
            x_max = x_max.max(x);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
    }

    // Add 5% padding
    let x_range = x_max - x_min;
    let y_range = y_max - y_min;
    x_min -= 0.05 * x_range;
    x_max += 0.05 * x_range;
    y_min -= 0.05 * y_range;
    y_max += 0.05 * y_range;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, (config.font_family.as_str(), config.font_size_title))
        .margin(15)
        .x_label_area_size(80)
        .y_label_area_size(100)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .x_label_style((config.font_family.as_str(), config.font_size_label))
        .y_label_style((config.font_family.as_str(), config.font_size_label))
        .axis_desc_style((config.font_family.as_str(), config.font_size_label))
        .draw()?;

    // Draw data series
    for (data, label, color) in data_series {
        chart.draw_series(LineSeries::new(
            data.iter().map(|&(x, y)| (x, y)),
            color.stroke_width(3),
        ))?
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(3)));
    }

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font((config.font_family.as_str(), config.font_size_label))
        .draw()?;

    // Add equation annotations
    for (equation_text, (eq_x, eq_y)) in &equations {
        // Draw background box for equation
        chart.draw_series(std::iter::once(Rectangle::new(
            [(eq_x - 0.02 * (x_max - x_min), eq_y - 0.02 * (y_max - y_min)),
             (eq_x + 0.15 * (x_max - x_min), eq_y + 0.04 * (y_max - y_min))],
            WHITE.mix(0.95).filled(),
        )))?;

        // Draw equation text (clone to satisfy 'static requirement)
        let text = equation_text.clone();
        chart.draw_series(std::iter::once(Text::new(
            text,
            (*eq_x, *eq_y),
            (config.font_family.as_str(), config.font_size_equation).into_font(),
        )))?;
    }

    root.present()?;
    Ok(())
}

/// Create abundance evolution plot with equations
pub fn plot_abundance_evolution(
    filename: &str,
    evolution: &[crate::early_universe::network::AbundanceState],
    config: &PublicationConfig,
) -> Result<(), Box<dyn std::error::Error>> {

    use crate::early_universe::reactions::Nuclide;

    // Extract data for each nuclide
    let nuclides = [
        (Nuclide::Neutron, "n", RGBColor(0, 0, 255)),
        (Nuclide::Proton, "p", RGBColor(255, 0, 0)),
        (Nuclide::Deuterium, "²H", RGBColor(0, 150, 0)),
        (Nuclide::Helium4, "⁴He", RGBColor(200, 0, 200)),
        (Nuclide::Lithium7, "⁷Li", RGBColor(255, 165, 0)),
    ];

    let mut data_series = Vec::new();

    for (nuclide, label, color) in nuclides {
        let data: Vec<(f64, f64)> = evolution.iter()
            .map(|state| {
                let t = state.time;
                let x = state.mass_fraction(nuclide);
                (t, x)
            })
            .filter(|&(_, x)| x > 1e-15) // Filter very small values
            .collect();

        if !data.is_empty() {
            data_series.push((data, label, color));
        }
    }

    // Convert to references
    let data_refs: Vec<_> = data_series.iter()
        .map(|(data, label, color)| (data.as_slice(), *label, *color))
        .collect();

    // Add key equations
    let equations = vec![
        ("n/p = exp(-Δm/kT)".to_string(), (50.0, 0.3)),
        ("Y_p ≈ 2(n/p)/(1+n/p)".to_string(), (200.0, 0.15)),
    ];

    create_plot_with_equation(
        filename,
        data_refs,
        "Big Bang Nucleosynthesis: Abundance Evolution",
        "Time [s]",
        "Mass Fraction",
        equations,
        config,
    )
}
