//! Visualization utilities for cosmological data

pub mod plots_2d;
pub mod plotly_plots;
pub mod colors;

pub use plots_2d::PlotConfig;
pub use plotly_plots::{create_interactive_plot, create_loglog_interactive};
