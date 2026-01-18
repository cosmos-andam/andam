//! 2D plotting functionality using plotters

use plotters::prelude::*;
use std::error::Error;

/// Configuration for 2D plots
#[derive(Debug, Clone)]
pub struct PlotConfig {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub width: u32,
    pub height: u32,
    pub x_log: bool,
    pub y_log: bool,
}

impl Default for PlotConfig {
    fn default() -> Self {
        PlotConfig {
            title: "Plot".to_string(),
            x_label: "x".to_string(),
            y_label: "y".to_string(),
            width: 1024,
            height: 768,
            x_log: false,
            y_log: false,
        }
    }
}

/// Create a simple line plot
pub fn create_line_plot(
    filename: &str,
    data: &[(f64, f64)],
    config: &PlotConfig,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(filename, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find data ranges
    let x_min = data.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
    let x_max = data
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = data.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
    let y_max = data
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(&config.title, ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc(&config.x_label)
        .y_desc(&config.y_label)
        .draw()?;

    chart.draw_series(LineSeries::new(data.iter().map(|(x, y)| (*x, *y)), &BLUE))?;

    root.present()?;
    Ok(())
}

/// Create a log-log plot
pub fn create_loglog_plot(
    filename: &str,
    data: &[(f64, f64)],
    config: &PlotConfig,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(filename, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find data ranges (positive values only for log)
    let x_min = data
        .iter()
        .map(|(x, _)| *x)
        .filter(|x| *x > 0.0)
        .fold(f64::INFINITY, f64::min);
    let x_max = data
        .iter()
        .map(|(x, _)| *x)
        .filter(|x| *x > 0.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = data
        .iter()
        .map(|(_, y)| *y)
        .filter(|y| *y > 0.0)
        .fold(f64::INFINITY, f64::min);
    let y_max = data
        .iter()
        .map(|(_, y)| *y)
        .filter(|y| *y > 0.0)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(&config.title, ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())?;

    chart
        .configure_mesh()
        .x_desc(&config.x_label)
        .y_desc(&config.y_label)
        .draw()?;

    chart.draw_series(LineSeries::new(
        data.iter()
            .filter(|(x, y)| *x > 0.0 && *y > 0.0)
            .map(|(x, y)| (*x, *y)),
        &RED,
    ))?;

    root.present()?;
    Ok(())
}

/// Create a multi-line plot with legend
pub fn create_multiline_plot(
    filename: &str,
    datasets: &[(&[(f64, f64)], &str)],
    config: &PlotConfig,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(filename, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find global data ranges
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for (data, _) in datasets {
        for (x, y) in *data {
            x_min = x_min.min(*x);
            x_max = x_max.max(*x);
            y_min = y_min.min(*y);
            y_max = y_max.max(*y);
        }
    }

    let mut chart = ChartBuilder::on(&root)
        .caption(&config.title, ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc(&config.x_label)
        .y_desc(&config.y_label)
        .draw()?;

    let colors = [&BLUE, &RED, &GREEN, &MAGENTA, &CYAN];

    for (i, (data, label)) in datasets.iter().enumerate() {
        let color = colors[i % colors.len()];
        chart
            .draw_series(LineSeries::new(data.iter().map(|(x, y)| (*x, *y)), color))?
            .label(*label)
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
    }

    chart.configure_series_labels().border_style(BLACK).draw()?;

    root.present()?;
    Ok(())
}
