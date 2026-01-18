# Phase 1: Foundation (Weeks 1-3)

## Overview
Establish the foundational infrastructure for the cosmology Rust crate, including project structure, constants, unit system, basic Friedmann equation solver, and simple plotting capabilities.

---

## Objectives
- [x] Set up complete project structure
- [x] Implement physical constants module
- [x] Implement unit conversion system
- [x] Create basic Friedmann equation solver
- [x] Implement simple 2D plotting functionality
- [x] Set up testing framework
- [x] Create initial examples

---

## Task 1: Project Initialization

### Step 1.1: Create Project Structure
```bash
cargo new andam --lib
cd andam
```

### Step 1.2: Create Directory Structure
```bash
mkdir -p src/{fundamental,dynamics,observations,cmb,early_universe,structure,perturbations,advanced,visualization}
mkdir -p examples
mkdir -p tests
mkdir -p benches
mkdir -p docs
```

### Step 1.3: Set Up Cargo.toml

Create the main `Cargo.toml` with initial dependencies:

```toml
[package]
name = "andam"
version = "0.1.0"
edition = "2021"
authors = ["Andam Contributors"]
license = "MIT OR Apache-2.0"
description = "A comprehensive Rust library for cosmological calculations and visualizations"
repository = "https://github.com/cosmos-andam/andam"
keywords = ["cosmology", "astrophysics", "physics", "simulation"]
categories = ["science", "simulation"]

[dependencies]
# Numerical computation
ndarray = "0.15"
nalgebra = "0.32"
num-complex = "0.4"
num-traits = "0.2"

# ODE solvers
ode_solvers = "0.4"

# Plotting
plotly = { version = "0.8", features = ["kaleido"] }
plotters = "0.3"
plotters-backend = "0.3"

# Utilities
thiserror = "1.0"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"

[dev-dependencies]
criterion = "0.5"
approx = "0.5"

[[bench]]
name = "friedmann_bench"
harness = false

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
```

---

## Task 2: Constants Module

### Step 2.1: Create `src/constants.rs`

Implement comprehensive physical and cosmological constants:

```rust
//! Physical and cosmological constants
//! 
//! All constants are in SI units unless otherwise specified

use std::f64::consts::PI;

/// Speed of light in vacuum [m/s]
pub const C: f64 = 2.99792458e8;

/// Gravitational constant [m³/(kg·s²)]
pub const G: f64 = 6.67430e-11;

/// Planck constant [J·s]
pub const H_PLANCK: f64 = 6.62607015e-34;

/// Reduced Planck constant (ℏ) [J·s]
pub const HBAR: f64 = 1.054571817e-34;

/// Boltzmann constant [J/K]
pub const K_B: f64 = 1.380649e-23;

/// Stefan-Boltzmann constant [W/(m²·K⁴)]
pub const SIGMA_SB: f64 = 5.670374419e-8;

/// Electron mass [kg]
pub const M_E: f64 = 9.1093837015e-31;

/// Proton mass [kg]
pub const M_P: f64 = 1.67262192369e-27;

/// Neutron mass [kg]
pub const M_N: f64 = 1.67492749804e-27;

/// Atomic mass unit [kg]
pub const AMU: f64 = 1.66053906660e-27;

/// Elementary charge [C]
pub const E_CHARGE: f64 = 1.602176634e-19;

/// Thomson scattering cross-section [m²]
pub const SIGMA_T: f64 = 6.6524587321e-29;

/// Fine structure constant (dimensionless)
pub const ALPHA: f64 = 7.2973525693e-3;

// Astronomical constants

/// Astronomical unit [m]
pub const AU: f64 = 1.495978707e11;

/// Parsec [m]
pub const PARSEC: f64 = 3.0856775814913673e16;

/// Solar mass [kg]
pub const M_SUN: f64 = 1.98847e30;

/// Solar luminosity [W]
pub const L_SUN: f64 = 3.828e26;

/// Solar radius [m]
pub const R_SUN: f64 = 6.957e8;

// Cosmological constants

/// Hubble constant (Planck 2018) [km/s/Mpc]
pub const H0_PLANCK: f64 = 67.66;

/// Critical density coefficient [kg/(m³·(km/s/Mpc)²)]
pub const RHO_CRIT_COEFF: f64 = 3.0 / (8.0 * PI * G);

/// CMB temperature today [K]
pub const T_CMB: f64 = 2.7255;

/// CMB photon number density today [m⁻³]
pub const N_GAMMA: f64 = 4.105e8;

/// Neutrino temperature today [K]
pub const T_NU: f64 = 1.95;

// Energy conversions

/// Electron volt to Joules
pub const EV_TO_J: f64 = 1.602176634e-19;

/// Joules to electron volt
pub const J_TO_EV: f64 = 6.241509074e18;

// Time conversions

/// Year in seconds
pub const YEAR: f64 = 3.15576e7;

/// Gigayear in seconds
pub const GYR: f64 = 3.15576e16;

/// Hubble time [s]
pub fn hubble_time(h0_km_s_mpc: f64) -> f64 {
 PARSEC * 1e6 / (h0_km_s_mpc * 1e3)
}

/// Critical density today [kg/m³]
pub fn critical_density(h0_km_s_mpc: f64) -> f64 {
 let h_si = h0_km_s_mpc * 1e3 / (PARSEC * 1e6); // Convert to SI
 3.0 * h_si * h_si / (8.0 * PI * G)
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_speed_of_light() {
 assert_eq!(C, 2.99792458e8);
 }
 
 #[test]
 fn test_hubble_time() {
 let t_h = hubble_time(70.0);
 let expected = 4.4e17; // approximately 14 Gyr in seconds
 assert_relative_eq!(t_h, expected, epsilon = 1e16);
 }
 
 #[test]
 fn test_critical_density() {
 let rho_c = critical_density(70.0);
 assert!(rho_c > 0.0);
 assert!(rho_c < 1e-25); // Should be very small
 }
}
```

---

## Task 3: Units Module

### Step 3.1: Create `src/units.rs`

Implement unit conversion system:

```rust
//! Unit conversion utilities for cosmology
//! 
//! Provides conversions between common astronomical and cosmological units

use crate::constants::*;

/// Length units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Length {
 Meter(f64),
 Kilometer(f64),
 AU(f64),
 Parsec(f64),
 Kiloparsec(f64),
 Megaparsec(f64),
 Gigaparsec(f64),
}

impl Length {
 /// Convert to meters
 pub fn to_meters(&self) -> f64 {
 match self {
 Length::Meter(x) => *x,
 Length::Kilometer(x) => x * 1e3,
 Length::AU(x) => x * AU,
 Length::Parsec(x) => x * PARSEC,
 Length::Kiloparsec(x) => x * PARSEC * 1e3,
 Length::Megaparsec(x) => x * PARSEC * 1e6,
 Length::Gigaparsec(x) => x * PARSEC * 1e9,
 }
 }
 
 /// Convert to megaparsecs
 pub fn to_mpc(&self) -> f64 {
 self.to_meters() / (PARSEC * 1e6)
 }
 
 /// Convert to parsecs
 pub fn to_pc(&self) -> f64 {
 self.to_meters() / PARSEC
 }
}

/// Mass units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mass {
 Kilogram(f64),
 SolarMass(f64),
 EarthMass(f64),
}

impl Mass {
 /// Convert to kilograms
 pub fn to_kg(&self) -> f64 {
 match self {
 Mass::Kilogram(x) => *x,
 Mass::SolarMass(x) => x * M_SUN,
 Mass::EarthMass(x) => x * 5.972e24,
 }
 }
 
 /// Convert to solar masses
 pub fn to_solar_masses(&self) -> f64 {
 self.to_kg() / M_SUN
 }
}

/// Time units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Time {
 Second(f64),
 Year(f64),
 Kiloyear(f64),
 Megayear(f64),
 Gigayear(f64),
}

impl Time {
 /// Convert to seconds
 pub fn to_seconds(&self) -> f64 {
 match self {
 Time::Second(x) => *x,
 Time::Year(x) => x * YEAR,
 Time::Kiloyear(x) => x * YEAR * 1e3,
 Time::Megayear(x) => x * YEAR * 1e6,
 Time::Gigayear(x) => x * YEAR * 1e9,
 }
 }
 
 /// Convert to gigayears
 pub fn to_gyr(&self) -> f64 {
 self.to_seconds() / GYR
 }
}

/// Energy units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Energy {
 Joule(f64),
 ElectronVolt(f64),
 KiloElectronVolt(f64),
 MegaElectronVolt(f64),
 GeV(f64),
}

impl Energy {
 /// Convert to joules
 pub fn to_joules(&self) -> f64 {
 match self {
 Energy::Joule(x) => *x,
 Energy::ElectronVolt(x) => x * EV_TO_J,
 Energy::KiloElectronVolt(x) => x * EV_TO_J * 1e3,
 Energy::MegaElectronVolt(x) => x * EV_TO_J * 1e6,
 Energy::GeV(x) => x * EV_TO_J * 1e9,
 }
 }
 
 /// Convert to eV
 pub fn to_ev(&self) -> f64 {
 self.to_joules() * J_TO_EV
 }
}

/// Temperature conversion utilities
pub mod temperature {
 use super::*;
 
 /// Convert temperature to energy (E = kT)
 pub fn temp_to_energy_ev(temp_kelvin: f64) -> f64 {
 (K_B * temp_kelvin) * J_TO_EV
 }
 
 /// Convert energy to temperature (T = E/k)
 pub fn energy_to_temp(energy_ev: f64) -> f64 {
 (energy_ev * EV_TO_J) / K_B
 }
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_length_conversions() {
 let d = Length::Megaparsec(1.0);
 assert_relative_eq!(d.to_mpc(), 1.0, epsilon = 1e-10);
 assert_relative_eq!(d.to_pc(), 1e6, epsilon = 1e-4);
 }
 
 #[test]
 fn test_time_conversions() {
 let t = Time::Gigayear(1.0);
 assert_relative_eq!(t.to_gyr(), 1.0, epsilon = 1e-10);
 }
 
 #[test]
 fn test_temperature_energy() {
 let temp = 1e4; // 10,000 K
 let energy = temperature::temp_to_energy_ev(temp);
 let back = temperature::energy_to_temp(energy);
 assert_relative_eq!(temp, back, epsilon = 1e-6);
 }
}
```

---

## Task 4: Basic Universe Model

### Step 4.1: Create `src/dynamics/mod.rs`

```rust
//! Cosmic dynamics and evolution
//! 
//! Implements Friedmann equations and universe evolution

pub mod friedmann;
pub mod components;

pub use friedmann::Universe;
pub use components::{Component, ComponentType};
```

### Step 4.2: Create `src/dynamics/components.rs`

```rust
//! Universe components and their properties

/// Types of components in the universe
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ComponentType {
 Matter,
 Radiation,
 DarkEnergy,
 Curvature,
}

/// Represents a component of the universe with equation of state
#[derive(Debug, Clone)]
pub struct Component {
 pub name: String,
 pub component_type: ComponentType,
 pub omega_0: f64, // Density parameter today
 pub w: f64, // Equation of state parameter
}

impl Component {
 /// Create a matter component (w = 0)
 pub fn matter(omega_0: f64) -> Self {
 Component {
 name: "Matter".to_string(),
 component_type: ComponentType::Matter,
 omega_0,
 w: 0.0,
 }
 }
 
 /// Create a radiation component (w = 1/3)
 pub fn radiation(omega_0: f64) -> Self {
 Component {
 name: "Radiation".to_string(),
 component_type: ComponentType::Radiation,
 omega_0,
 w: 1.0 / 3.0,
 }
 }
 
 /// Create a cosmological constant (w = -1)
 pub fn dark_energy(omega_0: f64) -> Self {
 Component {
 name: "Dark Energy".to_string(),
 component_type: ComponentType::DarkEnergy,
 omega_0,
 w: -1.0,
 }
 }
 
 /// Create a curvature component (w = -1/3)
 pub fn curvature(omega_0: f64) -> Self {
 Component {
 name: "Curvature".to_string(),
 component_type: ComponentType::Curvature,
 omega_0,
 w: -1.0 / 3.0,
 }
 }
 
 /// Energy density at scale factor a, normalized to today
 /// ρ(a) / ρ_0 = a^(-3(1+w))
 pub fn density_evolution(&self, a: f64) -> f64 {
 self.omega_0 * a.powf(-3.0 * (1.0 + self.w))
 }
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_matter_scaling() {
 let matter = Component::matter(0.3);
 assert_relative_eq!(matter.density_evolution(0.5), 0.3 * 8.0, epsilon = 1e-10);
 }
 
 #[test]
 fn test_radiation_scaling() {
 let radiation = Component::radiation(1e-4);
 assert_relative_eq!(radiation.density_evolution(0.5), 1e-4 * 16.0, epsilon = 1e-14);
 }
}
```

### Step 4.3: Create `src/dynamics/friedmann.rs`

```rust
//! Friedmann equation solver and universe evolution

use crate::constants::*;
use crate::dynamics::components::{Component, ComponentType};
use std::f64::consts::PI;

/// Represents a cosmological model with multiple components
#[derive(Debug, Clone)]
pub struct Universe {
 /// Hubble constant today [km/s/Mpc]
 pub h0: f64,
 /// Components (matter, radiation, dark energy, curvature)
 pub components: Vec<Component>,
}

impl Universe {
 /// Create a new universe with given Hubble constant
 pub fn new(h0: f64) -> Self {
 Universe {
 h0,
 components: Vec::new(),
 }
 }
 
 /// Add a component to the universe
 pub fn add_component(&mut self, component: Component) {
 self.components.push(component);
 }
 
 /// Create a benchmark ΛCDM model (Planck 2018)
 pub fn benchmark() -> Self {
 let mut universe = Universe::new(67.66);
 universe.add_component(Component::matter(0.3111));
 universe.add_component(Component::radiation(9.3e-5));
 universe.add_component(Component::dark_energy(0.6889));
 universe
 }
 
 /// Create an Einstein-de Sitter universe (Ω_m = 1, flat)
 pub fn einstein_de_sitter(h0: f64) -> Self {
 let mut universe = Universe::new(h0);
 universe.add_component(Component::matter(1.0));
 universe
 }
 
 /// Total Omega (should be 1 for flat universe)
 pub fn omega_total(&self) -> f64 {
 self.components.iter().map(|c| c.omega_0).sum()
 }
 
 /// Hubble parameter at scale factor a: H(a) / H_0
 pub fn hubble_normalized(&self, a: f64) -> f64 {
 let sum: f64 = self.components
 .iter()
 .map(|c| c.density_evolution(a))
 .sum();
 sum.sqrt()
 }
 
 /// Hubble parameter at scale factor a [km/s/Mpc]
 pub fn hubble(&self, a: f64) -> f64 {
 self.h0 * self.hubble_normalized(a)
 }
 
 /// Hubble parameter at redshift z [km/s/Mpc]
 pub fn hubble_z(&self, z: f64) -> f64 {
 let a = 1.0 / (1.0 + z);
 self.hubble(a)
 }
 
 /// Deceleration parameter q(a) = -ä·a/ȧ²
 pub fn deceleration(&self, a: f64) -> f64 {
 let mut sum_rho_1_w = 0.0;
 let mut sum_rho = 0.0;
 
 for component in &self.components {
 let rho = component.density_evolution(a);
 sum_rho += rho;
 sum_rho_1_w += rho * (1.0 + component.w);
 }
 
 0.5 * sum_rho_1_w / sum_rho
 }
 
 /// Age of universe at scale factor a [Gyr]
 pub fn age(&self, a: f64) -> f64 {
 // Integrate dt/da = 1/(a·H(a))
 let n_steps = 1000;
 let da = a / n_steps as f64;
 
 let mut age = 0.0;
 for i in 1..=n_steps {
 let a_i = i as f64 * da;
 let h = self.hubble(a_i); // km/s/Mpc
 let h_si = h * 1e3 / (PARSEC * 1e6); // Convert to SI
 let dt = 1.0 / (a_i * h_si); // seconds
 age += dt;
 }
 
 age / GYR // Convert to Gyr
 }
 
 /// Age of universe today [Gyr]
 pub fn age_today(&self) -> f64 {
 self.age(1.0)
 }
 
 /// Critical density today [kg/m³]
 pub fn critical_density(&self) -> f64 {
 critical_density(self.h0)
 }
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_benchmark_universe() {
 let universe = Universe::benchmark();
 let omega_tot = universe.omega_total();
 assert_relative_eq!(omega_tot, 1.0, epsilon = 1e-3);
 }
 
 #[test]
 fn test_einstein_de_sitter() {
 let universe = Universe::einstein_de_sitter(70.0);
 assert_relative_eq!(universe.omega_total(), 1.0, epsilon = 1e-10);
 
 // H(a) = H_0 * a^(-3/2) for matter-only
 assert_relative_eq!(universe.hubble_normalized(0.5), 2.828427, epsilon = 1e-5);
 }
 
 #[test]
 fn test_age_calculation() {
 let universe = Universe::benchmark();
 let age = universe.age_today();
 // Should be around 13.8 Gyr
 assert!(age > 13.0 && age < 14.5);
 }
}
```

---

## Task 5: Simple Plotting

### Step 5.1: Create `src/visualization/mod.rs`

```rust
//! Visualization utilities for cosmological data

pub mod plots_2d;
pub mod colors;

pub use plots_2d::PlotConfig;
```

### Step 5.2: Create `src/visualization/plots_2d.rs`

```rust
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
 let root = BitMapBackend::new(filename, (config.width, config.height))
 .into_drawing_area();
 root.fill(&WHITE)?;
 
 // Find data ranges
 let x_min = data.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
 let x_max = data.iter().map(|(x, _)| *x).fold(f64::NEG_INFINITY, f64::max);
 let y_min = data.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
 let y_max = data.iter().map(|(_, y)| *y).fold(f64::NEG_INFINITY, f64::max);
 
 let mut chart = ChartBuilder::on(&root)
 .caption(&config.title, ("sans-serif", 40).into_font())
 .margin(10)
 .x_label_area_size(40)
 .y_label_area_size(50)
 .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
 
 chart.configure_mesh()
 .x_desc(&config.x_label)
 .y_desc(&config.y_label)
 .draw()?;
 
 chart.draw_series(LineSeries::new(
 data.iter().map(|(x, y)| (*x, *y)),
 &BLUE,
 ))?;
 
 root.present()?;
 Ok(())
}

/// Create a log-log plot
pub fn create_loglog_plot(
 filename: &str,
 data: &[(f64, f64)],
 config: &PlotConfig,
) -> Result<(), Box<dyn Error>> {
 let root = BitMapBackend::new(filename, (config.width, config.height))
 .into_drawing_area();
 root.fill(&WHITE)?;
 
 // Find data ranges (positive values only for log)
 let x_min = data.iter()
 .map(|(x, _)| *x)
 .filter(|x| *x > 0.0)
 .fold(f64::INFINITY, f64::min);
 let x_max = data.iter()
 .map(|(x, _)| *x)
 .filter(|x| *x > 0.0)
 .fold(f64::NEG_INFINITY, f64::max);
 let y_min = data.iter()
 .map(|(_, y)| *y)
 .filter(|y| *y > 0.0)
 .fold(f64::INFINITY, f64::min);
 let y_max = data.iter()
 .map(|(_, y)| *y)
 .filter(|y| *y > 0.0)
 .fold(f64::NEG_INFINITY, f64::max);
 
 let mut chart = ChartBuilder::on(&root)
 .caption(&config.title, ("sans-serif", 40).into_font())
 .margin(10)
 .x_label_area_size(40)
 .y_label_area_size(50)
 .build_cartesian_2d(
 (x_min..x_max).log_scale(),
 (y_min..y_max).log_scale()
 )?;
 
 chart.configure_mesh()
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
 let root = BitMapBackend::new(filename, (config.width, config.height))
 .into_drawing_area();
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
 
 chart.configure_mesh()
 .x_desc(&config.x_label)
 .y_desc(&config.y_label)
 .draw()?;
 
 let colors = [&BLUE, &RED, &GREEN, &MAGENTA, &CYAN];
 
 for (i, (data, label)) in datasets.iter().enumerate() {
 let color = colors[i % colors.len()];
 chart.draw_series(LineSeries::new(
 data.iter().map(|(x, y)| (*x, *y)),
 color,
 ))?
 .label(*label)
 .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
 }
 
 chart.configure_series_labels()
 .border_style(&BLACK)
 .draw()?;
 
 root.present()?;
 Ok(())
}
```

### Step 5.3: Create `src/visualization/colors.rs`

```rust
//! Color schemes for cosmological visualizations

use plotters::style::RGBColor;

/// CMB color scheme (blue to red)
pub fn cmb_colormap(value: f64) -> RGBColor {
 // value should be in [0, 1]
 let v = value.max(0.0).min(1.0);
 
 if v < 0.5 {
 let t = v * 2.0;
 RGBColor(
 (255.0 * (1.0 - t)) as u8,
 0,
 (255.0 * t) as u8,
 )
 } else {
 let t = (v - 0.5) * 2.0;
 RGBColor(
 (255.0 * t) as u8,
 0,
 (255.0 * (1.0 - t)) as u8,
 )
 }
}

/// Viridis-like colormap
pub fn viridis(value: f64) -> RGBColor {
 let v = value.max(0.0).min(1.0);
 
 let r = (253.0 * v.powf(0.8)) as u8;
 let g = (231.0 * (1.0 - (1.0 - v).powf(2.0))) as u8;
 let b = (37.0 + 218.0 * v) as u8;
 
 RGBColor(r, g, b)
}
```

---

## Task 6: Library Entry Point

### Step 6.1: Create `src/lib.rs`

```rust
//! # Andam
//! 
//! A comprehensive Rust library for cosmological calculations and visualizations.
//! 
//! ## Features
//! 
//! - Friedmann equation solver
//! - Multiple universe components (matter, radiation, dark energy)
//! - Unit conversions
//! - 2D plotting capabilities
//! 
//! ## Example
//! 
//! ```rust
//! use andam::dynamics::Universe;
//! 
//! let universe = Universe::benchmark();
//! let age = universe.age_today();
//! println!("Age of universe: {:.2} Gyr", age);
//! ```

pub mod constants;
pub mod units;
pub mod dynamics;
pub mod visualization;

// Re-export commonly used items
pub use constants::*;
pub use units::{Length, Mass, Time, Energy};
pub use dynamics::Universe;
```

---

## Task 7: Create Examples

### Step 7.1: Create `examples/hubble_diagram.rs`

```rust
//! Example: Create a Hubble diagram showing velocity vs distance

use andam::dynamics::Universe;
use andam::visualization::plots_2d::{create_line_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 // Generate distance and velocity data
 let mut data = Vec::new();
 for i in 0..100 {
 let distance_mpc = i as f64 * 10.0; // 0 to 1000 Mpc
 let velocity = universe.h0 * distance_mpc; // v = H_0 * d
 data.push((distance_mpc, velocity));
 }
 
 let config = PlotConfig {
 title: "Hubble Diagram".to_string(),
 x_label: "Distance [Mpc]".to_string(),
 y_label: "Velocity [km/s]".to_string(),
 ..Default::default()
 };
 
 create_line_plot("hubble_diagram.png", &data, &config)?;
 println!("Created hubble_diagram.png");
 
 Ok(())
}
```

### Step 7.2: Create `examples/universe_evolution.rs`

```rust
//! Example: Plot scale factor evolution for different universe models

use andam::dynamics::{Universe, Component};
use andam::visualization::plots_2d::{create_multiline_plot, PlotConfig};
use andam::constants::GYR;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 // Create different universe models
 let benchmark = Universe::benchmark();
 let einstein_de_sitter = Universe::einstein_de_sitter(70.0);
 
 let mut lambda_only = Universe::new(70.0);
 lambda_only.add_component(Component::dark_energy(1.0));
 
 // Generate data
 let n_points = 200;
 let mut data_benchmark = Vec::new();
 let mut data_eds = Vec::new();
 let mut data_lambda = Vec::new();
 
 for i in 0..n_points {
 let a = 0.01 + (i as f64 / n_points as f64) * 0.99; // a from 0.01 to 1.0
 let t_bench = benchmark.age(a);
 let t_eds = einstein_de_sitter.age(a);
 let t_lambda = lambda_only.age(a);
 
 data_benchmark.push((t_bench, a));
 data_eds.push((t_eds, a));
 data_lambda.push((t_lambda, a));
 }
 
 let datasets = [
 (data_benchmark.as_slice(), "Benchmark ΛCDM"),
 (data_eds.as_slice(), "Einstein-de Sitter"),
 (data_lambda.as_slice(), "Λ-only"),
 ];
 
 let config = PlotConfig {
 title: "Universe Evolution: Scale Factor vs Time".to_string(),
 x_label: "Time [Gyr]".to_string(),
 y_label: "Scale Factor a(t)".to_string(),
 ..Default::default()
 };
 
 create_multiline_plot("universe_evolution.png", &datasets, &config)?;
 println!("Created universe_evolution.png");
 
 Ok(())
}
```

### Step 7.3: Create `examples/hubble_parameter.rs`

```rust
//! Example: Plot Hubble parameter evolution

use andam::dynamics::Universe;
use andam::visualization::plots_2d::{create_loglog_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 // Generate H(a) data
 let mut data = Vec::new();
 let n_points = 200;
 
 for i in 1..=n_points {
 let a = (i as f64) / (n_points as f64);
 let h = universe.hubble_normalized(a);
 data.push((a, h));
 }
 
 let config = PlotConfig {
 title: "Hubble Parameter Evolution".to_string(),
 x_label: "Scale Factor a".to_string(),
 y_label: "H(a) / H₀".to_string(),
 ..Default::default()
 };
 
 create_loglog_plot("hubble_evolution.png", &data, &config)?;
 println!("Created hubble_evolution.png");
 
 Ok(())
}
```

---

## Task 8: Testing

### Step 8.1: Create `tests/integration_tests.rs`

```rust
//! Integration tests for andam

use andam::dynamics::Universe;
use approx::assert_relative_eq;

#[test]
fn test_benchmark_universe_age() {
 let universe = Universe::benchmark();
 let age = universe.age_today();
 
 // Planck 2018 gives ~13.8 Gyr
 assert!(age > 13.0 && age < 14.5, "Age should be around 13.8 Gyr, got {}", age);
}

#[test]
fn test_omega_total() {
 let universe = Universe::benchmark();
 let omega_tot = universe.omega_total();
 
 assert_relative_eq!(omega_tot, 1.0, epsilon = 1e-3);
}

#[test]
fn test_hubble_at_present() {
 let universe = Universe::benchmark();
 let h_normalized = universe.hubble_normalized(1.0);
 
 assert_relative_eq!(h_normalized, 1.0, epsilon = 1e-10);
}
```

---

## Task 9: Documentation

### Step 9.1: Add README.md

Create a comprehensive README with installation and usage instructions.

### Step 9.2: Generate Documentation

```bash
cargo doc --no-deps --open
```

---

## Phase 1 Completion Checklist

- [ ] All modules compile without errors
- [ ] All tests pass (`cargo test`)
- [ ] Examples run successfully
- [ ] Documentation builds (`cargo doc`)
- [ ] Code follows Rust conventions (`cargo clippy`)
- [ ] Code is formatted (`cargo fmt`)

---

## Expected Outputs

After completing Phase 1, you should have:

1. [DONE] Working project structure
2. [DONE] Constants and units modules
3. [DONE] Basic Friedmann equation solver
4. [DONE] Simple 2D plotting
5. [DONE] Three working examples generating plots:
 - `hubble_diagram.png`
 - `universe_evolution.png`
 - `hubble_evolution.png`

---

## Next Steps

Proceed to **Phase 2: Core Cosmology** for:
- Full universe evolution with ODE solver
- Distance measures (luminosity, angular diameter)
- CMB recombination physics
- Enhanced visualization capabilities
