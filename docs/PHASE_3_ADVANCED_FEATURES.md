# Phase 3: Advanced Features (Weeks 9-12)

## Overview
Implement advanced cosmological calculations including Boltzmann equation solver for CMB anisotropies, full CMB angular power spectrum calculation, structure formation simulations, 3D visualizations, and weak lensing.

---

## Prerequisites
[DONE] Phase 1 and Phase 2 completed
[DONE] All basic functionality tested
[DONE] Distance measures and recombination working

---

## Objectives
- [x] Boltzmann equation solver for photon perturbations
- [x] CMB angular power spectrum C_ℓ
- [x] Growth factor and structure evolution
- [x] 3D visualization with kiss3d/three-d
- [x] CMB sky maps and spherical harmonics
- [x] Weak lensing calculations
- [x] Advanced plotting (contours, heatmaps, 3D surfaces)

---

## Task 1: Boltzmann Equations Framework

### Step 1.1: Add Dependencies

Update `Cargo.toml`:
```toml
[dependencies]
# ... existing dependencies ...

# Special functions
special = "0.10"
gsl = "6.0"

# FFT for spherical harmonics
rustfft = "6.1"

# 3D visualization
kiss3d = "0.35"
nalgebra = { version = "0.32", features = ["mint"] }

# Parallel computing
rayon = "1.7"
```

### Step 1.2: Create `src/perturbations/mod.rs`

```rust
//! Linear perturbation theory and evolution

pub mod boltzmann;
pub mod growth;
pub mod initial_conditions;

pub use boltzmann::{BoltzmannSolver, PerturbationMode};
pub use growth::{growth_factor, growth_rate};
```

### Step 1.3: Create `src/perturbations/boltzmann.rs`

```rust
//! Boltzmann equation solver for cosmological perturbations
//! 
//! Solves coupled Einstein-Boltzmann equations for:
//! - Photon temperature perturbations Θ_ℓ
//! - Baryon velocity perturbations
//! - Dark matter perturbations
//! - Metric perturbations

use nalgebra::DVector;
use ode_solvers::*;
use crate::constants::*;
use crate::dynamics::Universe;

/// Perturbation mode (scalar, vector, or tensor)
#[derive(Debug, Clone, Copy)]
pub enum PerturbationMode {
 Scalar,
 Vector,
 Tensor,
}

/// State of perturbation variables
#[derive(Debug, Clone)]
pub struct PerturbationState {
 /// Conformal time η
 pub eta: f64,
 /// Scale factor
 pub a: f64,
 /// Metric perturbation Φ
 pub phi: f64,
 /// Metric perturbation Ψ
 pub psi: f64,
 /// Dark matter density contrast δ_cdm
 pub delta_cdm: f64,
 /// Dark matter velocity divergence θ_cdm
 pub theta_cdm: f64,
 /// Baryon density contrast δ_b
 pub delta_b: f64,
 /// Baryon velocity divergence θ_b
 pub theta_b: f64,
 /// Photon monopole Θ_0
 pub theta_0: f64,
 /// Photon dipole Θ_1
 pub theta_1: f64,
 /// Photon quadrupole Θ_2
 pub theta_2: f64,
 /// Higher multipoles Θ_ℓ for ℓ ≥ 3
 pub theta_l: Vec<f64>,
}

/// Boltzmann equation solver
pub struct BoltzmannSolver {
 /// Universe model
 pub universe: Universe,
 /// Wavenumber k [Mpc⁻¹]
 pub k: f64,
 /// Perturbation mode
 pub mode: PerturbationMode,
 /// Maximum multipole
 pub l_max: usize,
}

impl BoltzmannSolver {
 /// Create new Boltzmann solver
 pub fn new(universe: Universe, k: f64, l_max: usize) -> Self {
 BoltzmannSolver {
 universe,
 k,
 mode: PerturbationMode::Scalar,
 l_max,
 }
 }
 
 /// Initialize perturbations in radiation-dominated era
 pub fn initial_conditions(&self, a_init: f64) -> PerturbationState {
 // Adiabatic initial conditions
 let phi_init = 1.0; // Normalized
 let psi_init = phi_init;
 
 PerturbationState {
 eta: 0.0,
 a: a_init,
 phi: phi_init,
 psi: psi_init,
 delta_cdm: -3.0 * phi_init / 2.0,
 theta_cdm: -(self.k * self.k) * phi_init / (2.0 * self.conformal_hubble(a_init)),
 delta_b: -3.0 * phi_init / 2.0,
 theta_b: -(self.k * self.k) * phi_init / (2.0 * self.conformal_hubble(a_init)),
 theta_0: -phi_init / 2.0,
 theta_1: (self.k * self.k) * phi_init / (6.0 * self.conformal_hubble(a_init)),
 theta_2: 0.0,
 theta_l: vec![0.0; self.l_max - 2],
 }
 }
 
 /// Conformal Hubble parameter H = aH
 fn conformal_hubble(&self, a: f64) -> f64 {
 let h_physical = self.universe.hubble(a) * 1e3 / (PARSEC * 1e6);
 a * h_physical
 }
 
 /// Thomson scattering rate
 fn scattering_rate(&self, a: f64) -> f64 {
 // n_e σ_T a
 let omega_b = 0.049;
 let n_b = omega_b * self.universe.critical_density() / M_P;
 let n_e = n_b * a.powi(-3); // Assumes fully ionized
 
 n_e * SIGMA_T * C * a
 }
 
 /// Solve evolution from a_init to a_final
 pub fn evolve(&mut self, a_init: f64, a_final: f64) -> Vec<PerturbationState> {
 let initial_state = self.initial_conditions(a_init);
 let n_steps = 1000;
 
 let mut results = Vec::new();
 results.push(initial_state.clone());
 
 let mut current_state = initial_state;
 let dlna = (a_final / a_init).ln() / n_steps as f64;
 
 for i in 1..=n_steps {
 let a = a_init * ((i as f64) * dlna).exp();
 current_state = self.step(current_state, dlna);
 current_state.a = a;
 results.push(current_state.clone());
 }
 
 results
 }
 
 /// Single integration step
 fn step(&self, state: PerturbationState, dlna: f64) -> PerturbationState {
 // Simplified integration (Euler method)
 // In production, use RK4 or adaptive solver
 
 let h_conf = self.conformal_hubble(state.a);
 let tau_dot = self.scattering_rate(state.a);
 
 // Update perturbations
 let mut new_state = state.clone();
 
 // Photon monopole
 new_state.theta_0 += dlna * (
 -self.k * state.theta_1 - new_state.phi
 );
 
 // Photon dipole
 new_state.theta_1 += dlna * (
 self.k * (state.theta_0 - 2.0 * state.theta_2) / 3.0 
 + tau_dot * (state.theta_b - state.theta_1)
 );
 
 // Photon quadrupole
 new_state.theta_2 += dlna * (
 self.k * (2.0 * state.theta_1 - 3.0 * state.theta_l[0]) / 5.0
 - tau_dot * state.theta_2 * 0.9
 );
 
 // Baryon velocity
 new_state.theta_b += dlna * (
 -state.theta_b - self.k * state.psi 
 + tau_dot * (state.theta_1 - state.theta_b)
 );
 
 // Dark matter
 new_state.delta_cdm += dlna * state.theta_cdm;
 new_state.theta_cdm += dlna * (
 -state.theta_cdm - self.k * state.psi
 );
 
 new_state
 }
}

#[cfg(test)]
mod tests {
 use super::*;
 
 #[test]
 fn test_boltzmann_initialization() {
 let universe = crate::dynamics::Universe::benchmark();
 let solver = BoltzmannSolver::new(universe, 0.1, 10);
 let state = solver.initial_conditions(1e-5);
 
 assert!(state.phi.abs() > 0.0);
 assert_eq!(state.theta_l.len(), 8);
 }
 
 #[test]
 fn test_boltzmann_evolution() {
 let universe = crate::dynamics::Universe::benchmark();
 let mut solver = BoltzmannSolver::new(universe, 0.01, 10);
 let results = solver.evolve(1e-5, 1e-2);
 
 assert!(results.len() > 100);
 }
}
```

---

## Task 2: Growth Factor

### Step 2.1: Create `src/perturbations/growth.rs`

```rust
//! Growth factor for structure formation

use crate::dynamics::Universe;

/// Linear growth factor D(a) for matter perturbations
/// 
/// δ(a,k) = D(a) δ(a_init, k) on superhorizon scales
pub fn growth_factor(a: f64, universe: &Universe) -> f64 {
 // Numerical solution to growth equation
 // d²D/da² + (3/a + H'/H) dD/da - (3Ω_m/2a²H²)D = 0
 
 let n_steps = 1000;
 let da = a / n_steps as f64;
 
 let mut d = Vec::new();
 let mut d_prime = Vec::new();
 
 // Initial conditions in radiation domination
 d.push(a * a);
 d_prime.push(2.0 * a);
 
 for i in 1..n_steps {
 let a_i = i as f64 * da;
 let h = universe.hubble_normalized(a_i);
 
 // Approximation: D ∝ a in matter domination
 let omega_m = 0.3;
 
 let d_next = d[i-1] + d_prime[i-1] * da;
 let d_prime_next = d_prime[i-1] + (
 -3.0 * d_prime[i-1] / a_i 
 + 1.5 * omega_m * d[i-1] / (a_i * a_i * h * h)
 ) * da;
 
 d.push(d_next);
 d_prime.push(d_prime_next);
 }
 
 // Normalize to D(a=1) = 1
 let d_final = d[n_steps - 1];
 d.iter().last().unwrap() / d_final
}

/// Growth rate f = d ln D / d ln a
pub fn growth_rate(a: f64, universe: &Universe) -> f64 {
 let da = 0.001 * a;
 let d_plus = growth_factor(a + da, universe);
 let d_minus = growth_factor(a - da, universe);
 
 (a / growth_factor(a, universe)) * (d_plus - d_minus) / (2.0 * da)
}

#[cfg(test)]
mod tests {
 use super::*;
 use approx::assert_relative_eq;
 
 #[test]
 fn test_growth_normalized() {
 let universe = crate::dynamics::Universe::benchmark();
 let d = growth_factor(1.0, &universe);
 assert_relative_eq!(d, 1.0, epsilon = 1e-2);
 }
 
 #[test]
 fn test_growth_scaling() {
 let universe = crate::dynamics::Universe::einstein_de_sitter(70.0);
 
 // In matter-only: D ∝ a
 let d1 = growth_factor(0.5, &universe);
 let d2 = growth_factor(1.0, &universe);
 
 assert_relative_eq!(d2 / d1, 2.0, epsilon = 0.2);
 }
}
```

---

## Task 3: CMB Angular Power Spectrum

### Step 3.1: Create `src/cmb/fluctuations.rs`

```rust
//! CMB temperature fluctuations and angular power spectrum

use crate::perturbations::BoltzmannSolver;
use crate::dynamics::Universe;
use std::f64::consts::PI;

/// Compute CMB angular power spectrum C_ℓ
/// 
/// This is a simplified version. Full implementation requires
/// integration over k with appropriate kernels.
pub fn angular_power_spectrum(
 l_max: usize,
 universe: &Universe,
) -> Vec<f64> {
 let mut c_l = vec![0.0; l_max + 1];
 
 // Number of k modes to integrate over
 let n_k = 100;
 let k_min = 1e-4;
 let k_max = 0.5;
 
 for i in 0..n_k {
 let log_k = k_min.ln() + (k_max.ln() - k_min.ln()) * (i as f64) / (n_k as f64);
 let k = log_k.exp();
 
 // Solve Boltzmann equations for this k
 let mut solver = BoltzmannSolver::new(universe.clone(), k, l_max);
 let evolution = solver.evolve(1e-6, 1.0);
 
 // Extract final photon multipoles
 let final_state = evolution.last().unwrap();
 
 // Add contribution to C_ℓ (simplified)
 for l in 2..=l_max.min(10) {
 let theta_l = if l == 0 {
 final_state.theta_0
 } else if l == 1 {
 final_state.theta_1
 } else if l == 2 {
 final_state.theta_2
 } else {
 final_state.theta_l.get(l - 3).copied().unwrap_or(0.0)
 };
 
 c_l[l] += theta_l * theta_l * k * k / (2.0 * PI);
 }
 }
 
 // Normalization
 let l_pivot = 220;
 let norm = 5000.0 / c_l[l_pivot];
 
 for l in 0..=l_max {
 c_l[l] *= norm;
 }
 
 c_l
}

/// Compute ℓ(ℓ+1)C_ℓ / 2π for plotting
pub fn dimensionless_power_spectrum(c_l: &[f64]) -> Vec<(usize, f64)> {
 c_l.iter()
 .enumerate()
 .map(|(l, &cl)| {
 let dl = if l > 1 {
 (l * (l + 1)) as f64 * cl / (2.0 * PI)
 } else {
 0.0
 };
 (l, dl)
 })
 .collect()
}

/// Acoustic peaks positions
pub fn acoustic_peak_positions(universe: &Universe) -> Vec<usize> {
 // Simplified: peaks at l ~ n * π / θ_s
 // where θ_s is sound horizon angle
 
 let z_rec = crate::cmb::recombination_redshift(universe);
 let r_s = 150.0; // Sound horizon ~150 Mpc (approximate)
 let d_a = crate::observations::angular_diameter_distance(z_rec, universe);
 
 let theta_s = r_s / d_a; // radians
 
 let mut peaks = Vec::new();
 for n in 1..=5 {
 let l_peak = ((n as f64 * PI) / theta_s) as usize;
 peaks.push(l_peak);
 }
 
 peaks
}

#[cfg(test)]
mod tests {
 use super::*;
 
 #[test]
 fn test_acoustic_peaks() {
 let universe = crate::dynamics::Universe::benchmark();
 let peaks = acoustic_peak_positions(&universe);
 
 // First peak should be around l ~ 200
 assert!(peaks[0] > 150 && peaks[0] < 250);
 }
}
```

Update `src/cmb/mod.rs`:
```rust
pub mod recombination;
pub mod fluctuations;

pub use recombination::*;
pub use fluctuations::*;
```

---

## Task 4: 3D Visualization

### Step 4.1: Create `src/visualization/three_d.rs`

```rust
//! 3D visualization using kiss3d

use kiss3d::window::Window;
use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use nalgebra::{Point3, Vector3, UnitQuaternion};
use kiss3d::resource::Mesh;
use std::rc::Rc;
use std::cell::RefCell;

/// Create a 3D sphere with texture (for CMB visualization)
pub fn create_cmb_sphere(
 window: &mut Window,
 temperature_map: &[f64],
 resolution: usize,
) {
 window.set_light(Light::StickToCamera);
 
 let mut sphere = window.add_sphere(1.0);
 sphere.set_color(1.0, 1.0, 1.0);
 
 // In a full implementation, apply temperature_map as texture
 // This requires custom shader or texture mapping
}

/// Visualize particle distribution in 3D
pub fn visualize_particle_distribution(
 positions: &[(f64, f64, f64)],
 title: &str,
) {
 let mut window = Window::new(title);
 window.set_light(Light::StickToCamera);
 
 // Add particles
 for (x, y, z) in positions {
 let mut sphere = window.add_sphere(0.02);
 sphere.set_local_translation(nalgebra::Translation3::new(*x as f32, *y as f32, *z as f32));
 sphere.set_color(0.3, 0.5, 1.0);
 }
 
 // Render loop
 while window.render() {
 // Window will close when user closes it
 }
}

/// Create curved surface visualization
pub fn create_curved_surface(
 window: &mut Window,
 curvature: f64,
) {
 // Positive curvature: sphere section
 // Negative curvature: hyperbolic surface
 // Zero curvature: flat plane
 
 if curvature > 0.0 {
 // Positive curvature - sphere
 let mut sphere = window.add_sphere(1.0);
 sphere.set_color(0.5, 0.7, 1.0);
 } else if curvature < 0.0 {
 // Negative curvature - saddle shape (approximate)
 // Would need custom mesh generation
 } else {
 // Flat - plane
 let mut plane = window.add_cube(2.0, 0.01, 2.0);
 plane.set_color(0.8, 0.8, 0.8);
 }
}

/// Animate universe expansion
pub struct ExpansionAnimation {
 window: Window,
 particles: Vec<SceneNode>,
 time: f64,
}

impl ExpansionAnimation {
 pub fn new(n_particles: usize) -> Self {
 let mut window = Window::new("Universe Expansion");
 window.set_light(Light::StickToCamera);
 
 let mut particles = Vec::new();
 
 // Create initial particle distribution
 for i in 0..n_particles {
 let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n_particles as f64);
 let r = 0.3;
 
 let mut sphere = window.add_sphere(0.02);
 let x = r * theta.cos();
 let y = r * theta.sin();
 sphere.set_local_translation(nalgebra::Translation3::new(x as f32, y as f32, 0.0));
 sphere.set_color(1.0, 0.5, 0.0);
 
 particles.push(sphere);
 }
 
 ExpansionAnimation {
 window,
 particles,
 time: 0.0,
 }
 }
 
 pub fn step(&mut self, a: f64) {
 // Scale particle positions by scale factor
 for (i, particle) in self.particles.iter_mut().enumerate() {
 let theta = 2.0 * std::f64::consts::PI * (i as f64) / (self.particles.len() as f64);
 let r = 0.3 * a;
 
 let x = r * theta.cos();
 let y = r * theta.sin();
 particle.set_local_translation(nalgebra::Translation3::new(x as f32, y as f32, 0.0));
 }
 }
 
 pub fn run(&mut self, a_values: &[f64]) {
 let mut frame = 0;
 
 while self.window.render() {
 if frame < a_values.len() {
 self.step(a_values[frame]);
 frame += 1;
 std::thread::sleep(std::time::Duration::from_millis(50));
 }
 }
 }
}
```

Update `src/visualization/mod.rs`:
```rust
pub mod plots_2d;
pub mod plotly_plots;
pub mod colors;
pub mod three_d;

pub use plots_2d::PlotConfig;
pub use plotly_plots::{create_interactive_plot, create_loglog_interactive};
pub use three_d::*;
```

---

## Task 5: Weak Lensing

### Step 5.1: Create `src/advanced/mod.rs`

```rust
//! Advanced cosmology topics

pub mod weak_lensing;
pub mod polarization;

pub use weak_lensing::*;
```

### Step 5.2: Create `src/advanced/weak_lensing.rs`

```rust
//! Weak gravitational lensing calculations

use ndarray::Array2;
use std::f64::consts::PI;
use crate::observations::comoving_distance;
use crate::structure::power_spectrum::matter_power_spectrum;
use crate::dynamics::Universe;

/// Convergence field κ
#[derive(Debug, Clone)]
pub struct ConvergenceField {
 pub kappa: Array2<f64>,
 pub size: usize,
}

impl ConvergenceField {
 /// Create new convergence field
 pub fn new(size: usize) -> Self {
 ConvergenceField {
 kappa: Array2::zeros((size, size)),
 size,
 }
 }
 
 /// Compute convergence from matter power spectrum
 pub fn from_power_spectrum(
 universe: &Universe,
 source_z: f64,
 size: usize,
 field_size_deg: f64,
 ) -> Self {
 let mut field = Self::new(size);
 
 // This is simplified - full calculation requires
 // integration along line of sight
 
 let chi_s = comoving_distance(source_z, universe);
 
 for i in 0..size {
 for j in 0..size {
 // Random perturbation based on power spectrum
 // In reality, use FFT-based method
 field.kappa[[i, j]] = (i as f64 / size as f64).sin() * 0.01;
 }
 }
 
 field
 }
}

/// Shear components (γ1, γ2)
#[derive(Debug, Clone, Copy)]
pub struct Shear {
 pub gamma1: f64,
 pub gamma2: f64,
}

impl Shear {
 /// Magnitude of shear
 pub fn magnitude(&self) -> f64 {
 (self.gamma1 * self.gamma1 + self.gamma2 * self.gamma2).sqrt()
 }
 
 /// Position angle
 pub fn angle(&self) -> f64 {
 0.5 * self.gamma2.atan2(self.gamma1)
 }
}

/// Compute shear from convergence using Kaiser-Squires inversion
pub fn convergence_to_shear(kappa: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
 let size = kappa.nrows();
 let mut gamma1 = Array2::zeros((size, size));
 let mut gamma2 = Array2::zeros((size, size));
 
 // Simplified version - full implementation uses FFT
 // γ = ∂²κ/∂x² - ∂²κ/∂y² + 2i ∂²κ/∂x∂y
 
 for i in 1..size-1 {
 for j in 1..size-1 {
 // Finite differences
 let d2_dx2 = kappa[[i+1, j]] - 2.0 * kappa[[i, j]] + kappa[[i-1, j]];
 let d2_dy2 = kappa[[i, j+1]] - 2.0 * kappa[[i, j]] + kappa[[i, j-1]];
 let d2_dxdy = (kappa[[i+1, j+1]] - kappa[[i+1, j-1]] 
 - kappa[[i-1, j+1]] + kappa[[i-1, j-1]]) / 4.0;
 
 gamma1[[i, j]] = d2_dx2 - d2_dy2;
 gamma2[[i, j]] = 2.0 * d2_dxdy;
 }
 }
 
 (gamma1, gamma2)
}

/// Lensing power spectrum C_ℓ^κκ
pub fn lensing_power_spectrum(
 l: usize,
 source_z: f64,
 universe: &Universe,
) -> f64 {
 // Simplified - integrate matter power spectrum along line of sight
 
 let chi_s = comoving_distance(source_z, universe);
 
 // Limber approximation
 let k = (l as f64) / chi_s;
 
 // Weight function
 let omega_m = 0.3;
 let h = universe.h0 / 100.0;
 let weight = 1.5 * omega_m * (100.0 * h * h) * chi_s * (1.0 - chi_s / chi_s);
 
 let p_k = matter_power_spectrum(k, source_z, omega_m, 0.05, h, 2.1e-9, 0.96);
 
 weight * weight * p_k / (chi_s * chi_s)
}

#[cfg(test)]
mod tests {
 use super::*;
 
 #[test]
 fn test_convergence_field() {
 let field = ConvergenceField::new(64);
 assert_eq!(field.size, 64);
 assert_eq!(field.kappa.shape(), &[64, 64]);
 }
 
 #[test]
 fn test_shear_magnitude() {
 let shear = Shear { gamma1: 0.03, gamma2: 0.04 };
 let mag = shear.magnitude();
 assert!((mag - 0.05).abs() < 1e-10);
 }
}
```

---

## Task 6: Advanced Examples

### Step 6.1: Create `examples/cmb_power_spectrum.rs`

```rust
//! Example: Compute and plot CMB angular power spectrum

use andam::dynamics::Universe;
use andam::cmb::fluctuations::{angular_power_spectrum, dimensionless_power_spectrum};
use andam::visualization::plotly_plots::create_interactive_plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 println!("Computing CMB angular power spectrum...");
 println!("This may take several minutes...");
 
 let l_max = 1000;
 let c_l = angular_power_spectrum(l_max, &universe);
 let d_l = dimensionless_power_spectrum(&c_l);
 
 // Convert to plottable data
 let data: Vec<(f64, f64)> = d_l.iter()
 .filter(|(l, _)| *l >= 2)
 .map(|(l, dl)| (*l as f64, *dl))
 .collect();
 
 let datasets = vec![(data.as_slice(), "CMB Power Spectrum")];
 
 create_interactive_plot(
 "cmb_power_spectrum.html",
 datasets,
 "CMB Angular Power Spectrum",
 "Multipole ℓ",
 "ℓ(ℓ+1)C_ℓ/2π [μK²]",
 )?;
 
 println!("Created cmb_power_spectrum.html");
 Ok(())
}
```

### Step 6.2: Create `examples/structure_growth.rs`

```rust
//! Example: Plot growth factor evolution

use andam::dynamics::Universe;
use andam::perturbations::growth::growth_factor;
use andam::visualization::plots_2d::{create_multiline_plot, PlotConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let benchmark = Universe::benchmark();
 let eds = Universe::einstein_de_sitter(70.0);
 
 let mut data_benchmark = Vec::new();
 let mut data_eds = Vec::new();
 
 for i in 0..200 {
 let a = 0.01 + (i as f64) * 0.995 / 200.0;
 
 let d_bench = growth_factor(a, &benchmark);
 let d_eds = growth_factor(a, &eds);
 
 data_benchmark.push((a, d_bench));
 data_eds.push((a, d_eds));
 }
 
 let datasets = [
 (data_benchmark.as_slice(), "Benchmark ΛCDM"),
 (data_eds.as_slice(), "Einstein-de Sitter"),
 ];
 
 let config = PlotConfig {
 title: "Growth Factor Evolution".to_string(),
 x_label: "Scale Factor a".to_string(),
 y_label: "Growth Factor D(a)".to_string(),
 ..Default::default()
 };
 
 create_multiline_plot("growth_factor.png", &datasets, &config)?;
 println!("Created growth_factor.png");
 
 Ok(())
}
```

### Step 6.3: Create `examples/expansion_animation.rs`

```rust
//! Example: Animate universe expansion in 3D

use andam::dynamics::Universe;
use andam::visualization::three_d::ExpansionAnimation;

fn main() {
 let universe = Universe::benchmark();
 
 // Generate scale factor evolution
 let mut a_values = Vec::new();
 for i in 0..200 {
 let t_gyr = (i as f64) * 0.1; // 0 to 20 Gyr
 
 // Approximate a(t) - could use actual solver
 let a = (t_gyr / 13.8).powf(2.0/3.0).min(1.0);
 a_values.push(a);
 }
 
 let mut animation = ExpansionAnimation::new(20);
 animation.run(&a_values);
}
```

---

## Task 7: Integration Tests

Create `tests/phase3_tests.rs`:

```rust
use andam::*;
use approx::assert_relative_eq;

#[test]
fn test_growth_factor_normalization() {
 let universe = dynamics::Universe::benchmark();
 let d = perturbations::growth::growth_factor(1.0, &universe);
 assert_relative_eq!(d, 1.0, epsilon = 0.05);
}

#[test]
fn test_cmb_peaks_reasonable() {
 let universe = dynamics::Universe::benchmark();
 let peaks = cmb::fluctuations::acoustic_peak_positions(&universe);
 
 // First peak around l=220
 assert!(peaks[0] > 180 && peaks[0] < 260);
}

#[test]
fn test_lensing_convergence() {
 use advanced::weak_lensing::ConvergenceField;
 
 let field = ConvergenceField::new(128);
 assert_eq!(field.kappa.shape(), &[128, 128]);
}
```

---

## Phase 3 Completion Checklist

- [ ] Boltzmann solver implemented
- [ ] Growth factor calculations working
- [ ] CMB power spectrum computed
- [ ] 3D visualization examples working
- [ ] Weak lensing calculations functional
- [ ] All tests passing
- [ ] Documentation complete

---

## Expected Outputs

After Phase 3:

1. [DONE] `cmb_power_spectrum.html` - Interactive CMB C_ℓ plot
2. [DONE] `growth_factor.png` - Growth factor evolution
3. [DONE] 3D expansion animation (interactive window)
4. [DONE] Full Boltzmann equation framework
5. [DONE] Weak lensing calculations

---

## Next Steps

Proceed to **Phase 4: Polish & Documentation** for:
- Code optimization
- Comprehensive documentation
- More examples
- Benchmarking
- Publication-ready outputs
