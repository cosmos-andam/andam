# Phase 2: Core Cosmology (Weeks 4-8)

## Overview
Implement core cosmological calculations including sophisticated distance measures, CMB recombination physics, matter power spectrum, and enhanced visualization capabilities with interactive plots and color-coded diagrams.

---

## Prerequisites
[DONE] Phase 1 completed successfully
[DONE] All Phase 1 tests passing
[DONE] Basic plotting functionality working

---

## Objectives
- [x] Implement sophisticated ODE solver for universe evolution
- [x] Distance measures (luminosity, angular diameter, comoving)
- [x] CMB recombination physics (Saha equation, ionization fraction)
- [x] Matter power spectrum basics
- [x] Interactive plots with Plotly
- [x] Enhanced color schemes
- [x] Comprehensive examples

---

## Task 1: Advanced Universe Evolution

### Step 1.1: Add ODE Solver Dependencies

Update `Cargo.toml`:
```toml
[dependencies]
# ... existing dependencies ...
ode_solvers = "0.4"
```

### Step 1.2: Create `src/dynamics/solver.rs`

```rust
//! ODE solver for universe evolution equations

use ode_solvers::*;
use nalgebra::base::Vector3;
use crate::dynamics::Universe;
use crate::constants::*;

type State = Vector3<f64>;
type StateDerivative = Vector3<f64>;

/// State vector: [a, ȧ, t]
pub struct UniverseODE {
 universe: Universe,
}

impl UniverseODE {
 pub fn new(universe: Universe) -> Self {
 UniverseODE { universe }
 }
}

impl ode_solvers::System<State> for UniverseODE {
 fn system(&self, _t: f64, y: &State, dy: &mut State) {
 let a = y[0];
 let a_dot = y[1];
 // let t = y[2]; // Not used in RHS
 
 // Compute H² from Friedmann equation
 let h_normalized = self.universe.hubble_normalized(a);
 let h0_si = self.universe.h0 * 1e3 / (PARSEC * 1e6);
 let h = h0_si * h_normalized;
 
 // Compute acceleration from acceleration equation
 let q = self.universe.deceleration(a);
 let a_ddot = -q * h * h * a;
 
 dy[0] = a_dot;
 dy[1] = a_ddot;
 dy[2] = 1.0; // dt/dt = 1
 }
}

/// Solve universe evolution from a_start to a_end
pub fn solve_evolution(
 universe: &Universe,
 a_start: f64,
 a_end: f64,
 n_steps: usize,
) -> Vec<(f64, f64, f64)> {
 let mut system = UniverseODE::new(universe.clone());
 
 // Initial conditions: [a, ȧ, t]
 let h0_si = universe.h0 * 1e3 / (PARSEC * 1e6);
 let h_initial = h0_si * universe.hubble_normalized(a_start);
 let a_dot_initial = a_start * h_initial;
 
 let y0 = State::new(a_start, a_dot_initial, 0.0);
 
 let mut stepper = Rk4::new(system, 0.0, y0, 1.0, 1e-6);
 
 // Estimate total time span
 let t_span = universe.age(a_end) - universe.age(a_start);
 let t_span_si = t_span * GYR;
 let dt = t_span_si / n_steps as f64;
 
 let mut results = Vec::new();
 results.push((a_start, a_dot_initial, 0.0));
 
 let mut t = 0.0;
 for _ in 0..n_steps {
 if let Ok(_) = stepper.integrate() {
 let state = stepper.y_out();
 results.push((state[0], state[1], state[2]));
 t += dt;
 }
 }
 
 results
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_evolution_solver() {
 let universe = Universe::einstein_de_sitter(70.0);
 let results = solve_evolution(&universe, 0.1, 1.0, 100);
 
 assert!(results.len() > 50);
 assert_relative_eq!(results.last().unwrap().0, 1.0, epsilon = 0.1);
 }
}
```

---

## Task 2: Distance Measures

### Step 2.1: Create `src/observations/mod.rs`

```rust
//! Observational cosmology tools

pub mod distances;
pub mod dark_matter;

pub use distances::{
 luminosity_distance,
 angular_diameter_distance,
 comoving_distance,
 distance_modulus,
};
```

### Step 2.2: Create `src/observations/distances.rs`

```rust
//! Cosmological distance measures
//! 
//! Implements various distance definitions used in cosmology:
//! - Comoving distance
//! - Luminosity distance
//! - Angular diameter distance
//! - Distance modulus

use crate::dynamics::Universe;
use crate::constants::*;
use std::f64::consts::PI;

/// Integrand for comoving distance: 1/H(z)
fn comoving_integrand(z: f64, universe: &Universe) -> f64 {
 1.0 / universe.hubble_z(z)
}

/// Comoving distance to redshift z [Mpc]
/// 
/// d_c = (c/H_0) ∫₀^z dz'/H(z')
pub fn comoving_distance(z: f64, universe: &Universe) -> f64 {
 let c_over_h0 = C / (universe.h0 * 1e3) * PARSEC * 1e-6; // Mpc
 
 // Simpson's rule integration
 let n = 1000;
 let dz = z / n as f64;
 let mut sum = comoving_integrand(0.0, universe) + comoving_integrand(z, universe);
 
 for i in 1..n {
 let zi = i as f64 * dz;
 let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
 sum += weight * comoving_integrand(zi, universe);
 }
 
 c_over_h0 * sum * dz / 3.0
}

/// Transverse comoving distance [Mpc]
/// 
/// Accounts for curvature:
/// - Flat: d_M = d_c
/// - Open: d_M = (c/H_0)/√|Ω_k| sinh(√|Ω_k| d_c H_0/c)
/// - Closed: d_M = (c/H_0)/√Ω_k sin(√Ω_k d_c H_0/c)
pub fn transverse_comoving_distance(z: f64, universe: &Universe) -> f64 {
 let d_c = comoving_distance(z, universe);
 let omega_k = 1.0 - universe.omega_total();
 let c_over_h0 = C / (universe.h0 * 1e3) * PARSEC * 1e-6;
 
 if omega_k.abs() < 1e-5 {
 // Flat universe
 d_c
 } else if omega_k > 0.0 {
 // Open universe
 let sqrt_ok = omega_k.sqrt();
 (c_over_h0 / sqrt_ok) * (sqrt_ok * d_c / c_over_h0).sinh()
 } else {
 // Closed universe
 let sqrt_ok = (-omega_k).sqrt();
 (c_over_h0 / sqrt_ok) * (sqrt_ok * d_c / c_over_h0).sin()
 }
}

/// Luminosity distance [Mpc]
/// 
/// d_L = (1+z) d_M
pub fn luminosity_distance(z: f64, universe: &Universe) -> f64 {
 (1.0 + z) * transverse_comoving_distance(z, universe)
}

/// Angular diameter distance [Mpc]
/// 
/// d_A = d_M / (1+z)
pub fn angular_diameter_distance(z: f64, universe: &Universe) -> f64 {
 transverse_comoving_distance(z, universe) / (1.0 + z)
}

/// Distance modulus (m - M)
/// 
/// μ = 5 log₁₀(d_L / 10 pc)
pub fn distance_modulus(z: f64, universe: &Universe) -> f64 {
 let d_l = luminosity_distance(z, universe);
 let d_l_pc = d_l * 1e6; // Convert Mpc to pc
 5.0 * (d_l_pc / 10.0).log10()
}

/// Comoving volume element dV_c/dz/dΩ [Mpc³/sr]
pub fn comoving_volume_element(z: f64, universe: &Universe) -> f64 {
 let d_h = C / (universe.h0 * 1e3) * PARSEC * 1e-6; // Hubble distance [Mpc]
 let d_m = transverse_comoving_distance(z, universe);
 let e_z = universe.hubble_z(z) / universe.h0;
 
 d_h * d_m * d_m / e_z
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_comoving_distance_z0() {
 let universe = Universe::benchmark();
 let d = comoving_distance(0.0, &universe);
 assert_relative_eq!(d, 0.0, epsilon = 1e-6);
 }
 
 #[test]
 fn test_distance_relations() {
 let universe = Universe::benchmark();
 let z = 1.0;
 
 let d_l = luminosity_distance(z, &universe);
 let d_a = angular_diameter_distance(z, &universe);
 
 // d_L = (1+z)² d_A
 assert_relative_eq!(d_l, (1.0 + z).powi(2) * d_a, epsilon = 1e-6);
 }
 
 #[test]
 fn test_flat_universe_distances() {
 let universe = Universe::benchmark();
 let z = 0.5;
 
 let d_c = comoving_distance(z, &universe);
 let d_m = transverse_comoving_distance(z, &universe);
 
 // Should be equal for flat universe
 assert_relative_eq!(d_c, d_m, epsilon = 1e-6);
 }
}
```

---

## Task 3: CMB Recombination Physics

### Step 3.1: Create `src/cmb/mod.rs`

```rust
//! Cosmic Microwave Background physics

pub mod recombination;
pub mod fluctuations;

pub use recombination::{
 saha_equation,
 ionization_fraction,
 recombination_redshift,
 optical_depth,
};
```

### Step 3.2: Create `src/cmb/recombination.rs`

```rust
//! CMB recombination physics
//! 
//! Implements the Saha equation and ionization fraction calculation

use crate::constants::*;
use std::f64::consts::PI;

/// Saha equation for ionization equilibrium
/// 
/// n_e n_p / n_H = (m_e k T / 2π ℏ²)^(3/2) exp(-B/kT)
/// 
/// Returns ionization fraction X_e
pub fn saha_equation(temp_k: f64, n_h: f64) -> f64 {
 // Hydrogen binding energy [J]
 let binding_energy = 13.6 * EV_TO_J;
 
 // Prefactor
 let prefactor = ((M_E * K_B * temp_k) / (2.0 * PI * HBAR * HBAR)).powf(1.5);
 
 // Boltzmann factor
 let boltzmann = (-binding_energy / (K_B * temp_k)).exp();
 
 // Solve quadratic equation for X_e
 // X_e² / (1 - X_e) = (1/n_H) * prefactor * exp(-B/kT)
 let a = 1.0;
 let b = prefactor * boltzmann / n_h;
 let c = -b;
 
 // Quadratic formula
 let discriminant = b * b - 4.0 * a * c;
 if discriminant < 0.0 {
 return 1.0; // Fully ionized
 }
 
 let x_e = (-b + discriminant.sqrt()) / (2.0 * a);
 x_e.min(1.0).max(0.0)
}

/// Ionization fraction as a function of redshift
/// 
/// Uses Saha equation with corrections
pub fn ionization_fraction(z: f64, universe: &crate::dynamics::Universe) -> f64 {
 let temp = T_CMB * (1.0 + z);
 
 // Baryon number density today [m⁻³]
 let omega_b = 0.049; // Approximate baryon density parameter
 let rho_b = omega_b * universe.critical_density();
 let n_b = rho_b / M_P;
 
 // Scale with redshift
 let n_h = n_b * (1.0 + z).powi(3);
 
 saha_equation(temp, n_h)
}

/// Recombination redshift (when X_e ~ 0.5)
pub fn recombination_redshift(universe: &crate::dynamics::Universe) -> f64 {
 // Iterate to find z where X_e = 0.5
 let mut z_low = 500.0;
 let mut z_high = 2000.0;
 
 for _ in 0..50 {
 let z_mid = (z_low + z_high) / 2.0;
 let x_e = ionization_fraction(z_mid, universe);
 
 if x_e > 0.5 {
 z_low = z_mid;
 } else {
 z_high = z_mid;
 }
 }
 
 (z_low + z_high) / 2.0
}

/// Optical depth to Thomson scattering
/// 
/// τ = ∫ n_e σ_T c dt
pub fn optical_depth(z: f64, universe: &crate::dynamics::Universe) -> f64 {
 let n_steps = 1000;
 let dz = z / n_steps as f64;
 
 let mut tau = 0.0;
 
 for i in 0..n_steps {
 let zi = i as f64 * dz;
 let x_e = ionization_fraction(zi, universe);
 
 // Baryon number density
 let omega_b = 0.049;
 let rho_b = omega_b * universe.critical_density();
 let n_b = rho_b / M_P;
 let n_e = x_e * n_b * (1.0 + zi).powi(3);
 
 // dt/dz
 let h = universe.hubble_z(zi) * 1e3 / (PARSEC * 1e6); // SI units
 let dt_dz = 1.0 / ((1.0 + zi) * h);
 
 tau += n_e * SIGMA_T * C * dt_dz * dz;
 }
 
 tau
}

/// Visibility function g(z) = -dτ/dz exp(-τ)
/// 
/// Peaks at last scattering surface
pub fn visibility_function(z: f64, universe: &crate::dynamics::Universe) -> f64 {
 let dz = 1.0;
 let tau = optical_depth(z, universe);
 let tau_plus = optical_depth(z + dz, universe);
 let dtau_dz = (tau_plus - tau) / dz;
 
 -dtau_dz * (-tau).exp()
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_saha_high_temp() {
 let temp = 1e5; // High temperature
 let n_h = 1e6;
 let x_e = saha_equation(temp, n_h);
 
 // Should be nearly fully ionized
 assert!(x_e > 0.99);
 }
 
 #[test]
 fn test_saha_low_temp() {
 let temp = 1000.0; // Low temperature
 let n_h = 1e6;
 let x_e = saha_equation(temp, n_h);
 
 // Should be mostly neutral
 assert!(x_e < 0.1);
 }
 
 #[test]
 fn test_recombination_redshift() {
 let universe = crate::dynamics::Universe::benchmark();
 let z_rec = recombination_redshift(&universe);
 
 // Should be around z ~ 1100
 assert!(z_rec > 900.0 && z_rec < 1400.0);
 }
}
```

---

## Task 4: Matter Power Spectrum Basics

### Step 4.1: Create `src/structure/mod.rs`

```rust
//! Structure formation and power spectra

pub mod power_spectrum;
pub mod transfer_function;

pub use power_spectrum::{
 primordial_power_spectrum,
 matter_power_spectrum,
};
```

### Step 4.2: Create `src/structure/power_spectrum.rs`

```rust
//! Matter power spectrum calculations

use std::f64::consts::PI;

/// Primordial power spectrum from inflation
/// 
/// P(k) = A_s (k/k_pivot)^(n_s - 1)
/// 
/// # Arguments
/// * `k` - Wavenumber [h/Mpc]
/// * `amplitude` - A_s (typically ~2.1e-9)
/// * `spectral_index` - n_s (typically ~0.96)
/// * `k_pivot` - Pivot scale [h/Mpc] (typically 0.05)
pub fn primordial_power_spectrum(
 k: f64,
 amplitude: f64,
 spectral_index: f64,
 k_pivot: f64,
) -> f64 {
 amplitude * (k / k_pivot).powf(spectral_index - 1.0)
}

/// Transfer function (simplified Eisenstein & Hu form)
/// 
/// T(k) suppresses power on small scales
pub fn transfer_function_eh(k: f64, omega_m: f64, omega_b: f64, h: f64) -> f64 {
 // Eisenstein & Hu (1998) fitting formula
 let omega_m_h2 = omega_m * h * h;
 let omega_b_h2 = omega_b * h * h;
 let theta_cmb = T_CMB / 2.7;
 
 // Sound horizon
 let s = 44.5 * (omega_m_h2 / theta_cmb.powi(4)).ln() 
 / (1.0 + 10.0 * omega_b_h2.powf(0.75)).sqrt();
 
 // Silk damping scale
 let k_silk = 1.6 * omega_b_h2.powf(0.52) * omega_m_h2.powf(0.73)
 * (1.0 + (10.4 * omega_m_h2).powf(-0.95));
 
 let q = k / (13.41 * omega_m_h2.powf(0.5));
 let l = (2.0 * 2.718 + 1.8 * q).ln();
 let c = 14.2 + 731.0 / (1.0 + 62.5 * q);
 
 let t = l / (l + c * q * q);
 
 t
}

/// Linear matter power spectrum P(k)
/// 
/// P(k) = A_s (k/k_pivot)^(n_s-1) T²(k) D²(z) (2π²/k³)
pub fn matter_power_spectrum(
 k: f64,
 z: f64,
 omega_m: f64,
 omega_b: f64,
 h: f64,
 amplitude: f64,
 spectral_index: f64,
) -> f64 {
 let p_prim = primordial_power_spectrum(k, amplitude, spectral_index, 0.05);
 let t = transfer_function_eh(k, omega_m, omega_b, h);
 
 // Growth factor (simplified for matter-dominated)
 let d = 1.0 / (1.0 + z);
 
 // Dimensional factor
 let factor = 2.0 * PI * PI / (k * k * k);
 
 p_prim * t * t * d * d * factor
}

/// Dimensionless power spectrum Δ²(k)
pub fn dimensionless_power(
 k: f64,
 z: f64,
 omega_m: f64,
 omega_b: f64,
 h: f64,
 amplitude: f64,
 spectral_index: f64,
) -> f64 {
 let p_k = matter_power_spectrum(k, z, omega_m, omega_b, h, amplitude, spectral_index);
 k * k * k * p_k / (2.0 * PI * PI)
}

use crate::constants::T_CMB;

#[cfg(test)]
mod tests {
 use super::*;
 
 #[test]
 fn test_primordial_power_at_pivot() {
 let p = primordial_power_spectrum(0.05, 2.1e-9, 0.96, 0.05);
 assert!((p - 2.1e-9).abs() < 1e-15);
 }
 
 #[test]
 fn test_transfer_function_limiting_cases() {
 // Should be ~1 at large scales (small k)
 let t_large = transfer_function_eh(1e-4, 0.3, 0.05, 0.7);
 assert!(t_large > 0.9 && t_large <= 1.0);
 
 // Should be <1 at small scales (large k)
 let t_small = transfer_function_eh(10.0, 0.3, 0.05, 0.7);
 assert!(t_small < 0.5);
 }
}
```

---

## Task 5: Interactive Plotting with Plotly

### Step 5.1: Create `src/visualization/plotly_plots.rs`

```rust
//! Interactive plotting using Plotly

use plotly::{
 common::{Mode, Title, Line, Marker},
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
```

Update `src/visualization/mod.rs`:
```rust
pub mod plots_2d;
pub mod plotly_plots;
pub mod colors;

pub use plots_2d::PlotConfig;
pub use plotly_plots::{create_interactive_plot, create_loglog_interactive};
```

---

## Task 6: Comprehensive Examples

### Step 6.1: Create `examples/distance_measures.rs`

```rust
//! Example: Compare different distance measures

use andam::dynamics::Universe;
use andam::observations::{
 luminosity_distance,
 angular_diameter_distance,
 comoving_distance,
};
use andam::visualization::plotly_plots::create_interactive_plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 let mut data_lum = Vec::new();
 let mut data_ang = Vec::new();
 let mut data_com = Vec::new();
 
 for i in 1..100 {
 let z = (i as f64) * 0.05; // z from 0.05 to 5.0
 let d_l = luminosity_distance(z, &universe);
 let d_a = angular_diameter_distance(z, &universe);
 let d_c = comoving_distance(z, &universe);
 
 data_lum.push((z, d_l));
 data_ang.push((z, d_a));
 data_com.push((z, d_c));
 }
 
 let datasets = vec![
 (data_lum.as_slice(), "Luminosity Distance"),
 (data_ang.as_slice(), "Angular Diameter Distance"),
 (data_com.as_slice(), "Comoving Distance"),
 ];
 
 create_interactive_plot(
 "distance_measures.html",
 datasets,
 "Cosmological Distance Measures",
 "Redshift z",
 "Distance [Mpc]",
 )?;
 
 println!("Created distance_measures.html");
 Ok(())
}
```

### Step 6.2: Create `examples/recombination.rs`

```rust
//! Example: Plot ionization fraction during recombination

use andam::dynamics::Universe;
use andam::cmb::recombination::{ionization_fraction, recombination_redshift};
use andam::visualization::plots_2d::{create_line_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 let mut data = Vec::new();
 
 // From z=2000 to z=500
 for i in 0..200 {
 let z = 2000.0 - (i as f64) * 7.5;
 let x_e = ionization_fraction(z, &universe);
 data.push((z, x_e));
 }
 
 let z_rec = recombination_redshift(&universe);
 println!("Recombination redshift: z = {:.1}", z_rec);
 
 let config = PlotConfig {
 title: "Ionization Fraction during Recombination".to_string(),
 x_label: "Redshift z".to_string(),
 y_label: "Ionization Fraction X_e".to_string(),
 ..Default::default()
 };
 
 create_line_plot("recombination.png", &data, &config)?;
 println!("Created recombination.png");
 
 Ok(())
}
```

### Step 6.3: Create `examples/matter_power_spectrum.rs`

```rust
//! Example: Plot matter power spectrum

use andam::structure::power_spectrum::dimensionless_power;
use andam::visualization::plotly_plots::create_loglog_interactive;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 // Planck 2018 parameters
 let omega_m = 0.3111;
 let omega_b = 0.0490;
 let h = 0.6766;
 let a_s = 2.1e-9;
 let n_s = 0.9665;
 
 let mut data_z0 = Vec::new();
 let mut data_z1 = Vec::new();
 let mut data_z10 = Vec::new();
 
 // k from 10^-4 to 1 h/Mpc
 for i in 0..200 {
 let log_k = -4.0 + (i as f64) * 0.025;
 let k = 10_f64.powf(log_k);
 
 let delta2_z0 = dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s);
 let delta2_z1 = dimensionless_power(k, 1.0, omega_m, omega_b, h, a_s, n_s);
 let delta2_z10 = dimensionless_power(k, 10.0, omega_m, omega_b, h, a_s, n_s);
 
 data_z0.push((k, delta2_z0));
 data_z1.push((k, delta2_z1));
 data_z10.push((k, delta2_z10));
 }
 
 let datasets = vec![
 (data_z0.as_slice(), "z = 0"),
 (data_z1.as_slice(), "z = 1"),
 (data_z10.as_slice(), "z = 10"),
 ];
 
 create_loglog_interactive(
 "matter_power_spectrum.html",
 datasets,
 "Matter Power Spectrum",
 "Wavenumber k [h/Mpc]",
 "Δ²(k)",
 )?;
 
 println!("Created matter_power_spectrum.html");
 Ok(())
}
```

---

## Task 7: Update Library Entry Point

Update `src/lib.rs`:

```rust
//! # Andam
//! 
//! A comprehensive Rust library for cosmological calculations and visualizations.

pub mod constants;
pub mod units;
pub mod dynamics;
pub mod observations;
pub mod cmb;
pub mod structure;
pub mod visualization;

// Re-export commonly used items
pub use constants::*;
pub use units::{Length, Mass, Time, Energy};
pub use dynamics::Universe;
pub use observations::{
 luminosity_distance,
 angular_diameter_distance,
 comoving_distance,
};
pub use cmb::recombination::{ionization_fraction, recombination_redshift};
pub use structure::power_spectrum::matter_power_spectrum;
```

---

## Task 8: Integration Tests

Create `tests/phase2_tests.rs`:

```rust
use andam::*;
use approx::assert_relative_eq;

#[test]
fn test_distance_consistency() {
 let universe = Universe::benchmark();
 let z = 1.0;
 
 let d_l = luminosity_distance(z, &universe);
 let d_a = angular_diameter_distance(z, &universe);
 
 // Etherington reciprocity relation
 assert_relative_eq!(d_l / d_a, (1.0 + z).powi(2), epsilon = 1e-6);
}

#[test]
fn test_recombination_physics() {
 let universe = Universe::benchmark();
 let z_rec = recombination_redshift(&universe);
 
 let x_e_rec = ionization_fraction(z_rec, &universe);
 
 // Should be around 0.5 at recombination
 assert!((x_e_rec - 0.5).abs() < 0.1);
}

#[test]
fn test_power_spectrum_shape() {
 use structure::power_spectrum::dimensionless_power;
 
 let omega_m = 0.3;
 let omega_b = 0.05;
 let h = 0.7;
 let a_s = 2.1e-9;
 let n_s = 0.96;
 
 // Power should increase with k on large scales
 let p1 = dimensionless_power(0.01, 0.0, omega_m, omega_b, h, a_s, n_s);
 let p2 = dimensionless_power(0.1, 0.0, omega_m, omega_b, h, a_s, n_s);
 
 assert!(p2 > p1);
}
```

---

## Phase 2 Completion Checklist

- [ ] All new modules compile
- [ ] All tests pass
- [ ] Distance measures work correctly
- [ ] Recombination physics implemented
- [ ] Power spectrum basics working
- [ ] Interactive plots generate correctly
- [ ] All examples run and produce output
- [ ] Documentation updated

---

## Expected Outputs

After Phase 2:

1. [DONE] `distance_measures.html` - Interactive distance comparison
2. [DONE] `recombination.png` - Ionization fraction plot
3. [DONE] `matter_power_spectrum.html` - Interactive power spectrum
4. [DONE] Working distance calculations
5. [DONE] CMB recombination physics
6. [DONE] Basic matter power spectrum

---

## Next Steps

Proceed to **Phase 3: Advanced Features** for:
- Full Boltzmann equation solver
- CMB angular power spectrum
- Structure formation simulations
- 3D visualizations
- Weak lensing calculations
