//! Visualization utilities for cosmological data

pub mod colors;
pub mod equation_plots;
pub mod plots_2d;

pub use equation_plots::{plot_abundance_evolution, PublicationConfig};
pub use plots_2d::PlotConfig;
