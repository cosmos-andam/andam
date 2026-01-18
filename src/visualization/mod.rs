//! Visualization utilities for cosmological data

pub mod colors;
pub mod equation_plots;
pub mod plotly_plots;
pub mod plots_2d;

pub use equation_plots::{plot_abundance_evolution, PublicationConfig};
pub use plotly_plots::{create_interactive_plot, create_loglog_interactive};
pub use plots_2d::PlotConfig;
