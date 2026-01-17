//! Interactive plotting using Plotly

use plotly::{
    common::{Mode, Title},
    layout::{Axis, Layout},
    Plot, Scatter,
};
use std::error::Error;

/// Create an interactive line plot
pub fn create_interactive_plot(
    filename: &str,
    datasets: Vec<(&[(f64, f64)], &str)>,
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn Error>> {
    let mut plot = Plot::new();

    for (data, name) in datasets {
        let x: Vec<f64> = data.iter().map(|(x, _)| *x).collect();
        let y: Vec<f64> = data.iter().map(|(_, y)| *y).collect();

        let trace = Scatter::new(x, y)
            .mode(Mode::Lines)
            .name(name);

        plot.add_trace(trace);
    }

    let layout = Layout::new()
        .title(Title::new(title))
        .x_axis(Axis::new().title(Title::new(x_label)))
        .y_axis(Axis::new().title(Title::new(y_label)));

    plot.set_layout(layout);
    plot.write_html(filename);

    Ok(())
}

/// Create a log-log plot
pub fn create_loglog_interactive(
    filename: &str,
    datasets: Vec<(&[(f64, f64)], &str)>,
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn Error>> {
    let mut plot = Plot::new();

    for (data, name) in datasets {
        let x: Vec<f64> = data.iter().map(|(x, _)| *x).collect();
        let y: Vec<f64> = data.iter().map(|(_, y)| *y).collect();

        let trace = Scatter::new(x, y)
            .mode(Mode::Lines)
            .name(name);

        plot.add_trace(trace);
    }

    let layout = Layout::new()
        .title(Title::new(title))
        .x_axis(Axis::new()
            .title(Title::new(x_label))
            .type_(plotly::layout::AxisType::Log))
        .y_axis(Axis::new()
            .title(Title::new(y_label))
            .type_(plotly::layout::AxisType::Log));

    plot.set_layout(layout);
    plot.write_html(filename);

    Ok(())
}
