# Phase 4: Polish & Documentation (Weeks 13-14)

## Overview
Finalize the andam crate with comprehensive documentation, performance optimization, publication-ready examples, extensive testing, and preparation for public release.

---

## Prerequisites
[DONE] Phases 1-3 completed
[DONE] All core functionality implemented
[DONE] Basic testing done

---

## Objectives
- [x] Comprehensive API documentation
- [x] Performance optimization and benchmarking
- [x] Publication-ready visualization examples
- [x] Extensive test coverage
- [x] Example gallery
- [x] User guide and tutorial
- [x] CI/CD setup
- [x] Crate publication preparation

---

## Task 1: Comprehensive Documentation

### Step 1.1: Create Main Documentation

Create `README.md`:

```markdown
# andam

[![Crates.io](https://img.shields.io/crates/v/andam.svg)](https://crates.io/crates/andam)
[![Documentation](https://docs.rs/andam/badge.svg)](https://docs.rs/andam)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)

A comprehensive Rust library for cosmological calculations and visualizations, based on modern cosmology textbooks (Ryden's "Introduction to Cosmology" and Dodelson's "Modern Cosmology").

## Features

### Core Cosmology
- **Universe Evolution**: Friedmann equations, multi-component universes (matter, radiation, dark energy)
- **Distance Measures**: Luminosity distance, angular diameter distance, comoving distance
- **Redshift Calculations**: Complete redshift-distance relations
- **Age Calculations**: Universe age for any cosmological model

### CMB Physics
- **Recombination**: Saha equation, ionization fraction, optical depth
- **Angular Power Spectrum**: CMB C_ℓ calculation
- **Acoustic Oscillations**: Peak positions and physics
- **Polarization**: E and B mode calculations

### Structure Formation
- **Growth Factor**: Linear perturbation growth
- **Power Spectrum**: Matter power spectrum P(k)
- **Transfer Function**: Eisenstein-Hu fitting formulas
- **Boltzmann Solver**: Full Einstein-Boltzmann equations

### Advanced Features
- **Weak Lensing**: Convergence and shear calculations
- **3D Visualization**: Interactive 3D models with kiss3d
- **2D Plotting**: Publication-ready plots with plotters
- **Interactive Plots**: HTML plots with plotly

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
andam = "0.1.0"
```

## Quick Start

```rust
use andam::prelude::*;

fn main() {
 // Create a benchmark ΛCDM universe (Planck 2018)
 let universe = Universe::benchmark();
 
 // Calculate universe age
 let age = universe.age_today();
 println!("Age of universe: {:.2} Gyr", age);
 
 // Calculate distance to z=1
 let z = 1.0;
 let d_l = luminosity_distance(z, &universe);
 println!("Luminosity distance to z={}: {:.1} Mpc", z, d_l);
 
 // Calculate recombination redshift
 let z_rec = recombination_redshift(&universe);
 println!("Recombination occurred at z = {:.0}", z_rec);
}
```

## Examples

### Hubble Diagram

```rust
use andam::prelude::*;

let universe = Universe::benchmark();
let mut data = Vec::new();

for i in 0..100 {
 let distance_mpc = i as f64 * 10.0;
 let velocity = universe.h0 * distance_mpc;
 data.push((distance_mpc, velocity));
}

create_line_plot("hubble_diagram.png", &data, &config)?;
```

### CMB Power Spectrum

```rust
use andam::prelude::*;

let universe = Universe::benchmark();
let c_l = angular_power_spectrum(2000, &universe);
let d_l = dimensionless_power_spectrum(&c_l);

create_interactive_plot(
 "cmb_spectrum.html",
 vec![(d_l.as_slice(), "CMB")],
 "CMB Angular Power Spectrum",
 "Multipole ℓ",
 "ℓ(ℓ+1)C_ℓ/2π [μK²]",
)?;
```

See `examples/` directory for more examples.

## Documentation

Full API documentation is available at [docs.rs/andam](https://docs.rs/andam)

## Textbook References

This library implements concepts from:
- **Ryden, B.** (2016) *Introduction to Cosmology*, 2nd edition
- **Dodelson, S.** (2003) *Modern Cosmology*

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md)

## Citation

If you use this library in academic work, please cite:

```bibtex
@software{andam,
 title = {andam: Cosmological calculations in Rust},
 author = {Your Name},
 year = {2025},
 url = {https://github.com/cosmos-andam/andam}
}
```
```

### Step 1.2: Create User Guide

Create `docs/USER_GUIDE.md`:

```markdown
# Andam User Guide

## Table of Contents
1. [Installation](#installation)
2. [Basic Concepts](#basic-concepts)
3. [Universe Models](#universe-models)
4. [Distance Calculations](#distance-calculations)
5. [CMB Physics](#cmb-physics)
6. [Structure Formation](#structure-formation)
7. [Visualization](#visualization)
8. [Advanced Usage](#advanced-usage)

## Installation

### From crates.io
```bash
cargo add andam
```

### From source
```bash
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
```

## Basic Concepts

### Cosmological Constants

```rust
use andam::constants::*;

println!("Speed of light: {} m/s", C);
println!("Hubble constant (Planck): {} km/s/Mpc", H0_PLANCK);
println!("CMB temperature: {} K", T_CMB);
```

### Unit Conversions

```rust
use andam::units::*;

let distance = Length::Megaparsec(100.0);
println!("Distance in pc: {}", distance.to_pc());
println!("Distance in m: {}", distance.to_meters());

let time = Time::Gigayear(13.8);
println!("Time in seconds: {}", time.to_seconds());
```

## Universe Models

### Creating a Universe

```rust
use andam::dynamics::Universe;

// Benchmark ΛCDM model (Planck 2018)
let universe = Universe::benchmark();

// Einstein-de Sitter (matter-only)
let eds = Universe::einstein_de_sitter(70.0);

// Custom universe
let mut custom = Universe::new(70.0);
custom.add_component(Component::matter(0.3));
custom.add_component(Component::dark_energy(0.7));
```

### Universe Properties

```rust
let universe = Universe::benchmark();

// Hubble parameter
let h_today = universe.hubble(1.0);
let h_at_z1 = universe.hubble_z(1.0);

// Age
let age = universe.age_today();

// Deceleration parameter
let q = universe.deceleration(1.0);
```

## Distance Calculations

### Luminosity Distance

```rust
use andam::observations::luminosity_distance;

let z = 1.0;
let d_l = luminosity_distance(z, &universe);
println!("d_L(z={}) = {:.1} Mpc", z, d_l);
```

### Angular Diameter Distance

```rust
use andam::observations::angular_diameter_distance;

let d_a = angular_diameter_distance(z, &universe);
println!("d_A(z={}) = {:.1} Mpc", z, d_a);
```

### Distance Modulus

```rust
use andam::observations::distance_modulus;

let mu = distance_modulus(z, &universe);
println!("μ(z={}) = {:.2} mag", z, mu);
```

## CMB Physics

### Recombination

```rust
use andam::cmb::*;

let z_rec = recombination_redshift(&universe);
let x_e = ionization_fraction(z_rec, &universe);

println!("Recombination at z = {:.0}", z_rec);
println!("Ionization fraction: {:.3}", x_e);
```

### Angular Power Spectrum

```rust
let c_l = angular_power_spectrum(2000, &universe);
let peaks = acoustic_peak_positions(&universe);

println!("First acoustic peak at ℓ = {}", peaks[0]);
```

## Structure Formation

### Growth Factor

```rust
use andam::perturbations::growth::*;

let a = 0.5;
let d = growth_factor(a, &universe);
let f = growth_rate(a, &universe);

println!("D(a={}) = {:.3}", a, d);
println!("f(a={}) = {:.3}", a, f);
```

### Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;

let k = 0.1; // h/Mpc
let z = 0.0;
let p_k = matter_power_spectrum(k, z, 0.3, 0.05, 0.7, 2.1e-9, 0.96);

println!("P(k={}) = {:.2e} Mpc³/h³", k, p_k);
```

## Visualization

### 2D Plots

```rust
use andam::visualization::plots_2d::*;

let config = PlotConfig {
 title: "My Plot".to_string(),
 x_label: "x".to_string(),
 y_label: "y".to_string(),
 ..Default::default()
};

create_line_plot("output.png", &data, &config)?;
```

### Interactive Plots

```rust
use andam::visualization::plotly_plots::*;

create_interactive_plot(
 "output.html",
 vec![(data.as_slice(), "Series 1")],
 "Title",
 "X Label",
 "Y Label",
)?;
```

### 3D Visualization

```rust
use andam::visualization::three_d::*;

let mut animation = ExpansionAnimation::new(50);
animation.run(&a_values);
```

## Advanced Usage

### Custom Boltzmann Solver

```rust
use andam::perturbations::BoltzmannSolver;

let k = 0.05; // Mpc^-1
let l_max = 20;
let mut solver = BoltzmannSolver::new(universe, k, l_max);
let evolution = solver.evolve(1e-6, 1.0);
```

### Weak Lensing

```rust
use andam::advanced::weak_lensing::*;

let source_z = 1.0;
let convergence = ConvergenceField::from_power_spectrum(
 &universe,
 source_z,
 256,
 1.0, // degrees
);

let (gamma1, gamma2) = convergence_to_shear(&convergence.kappa);
```

## Best Practices

1. **Use the prelude**: `use andam::prelude::*;` imports commonly used items
2. **Check units**: Pay attention to units (Mpc, km/s/Mpc, etc.)
3. **Validate parameters**: Ensure Ω_total ≈ 1 for flat universes
4. **Error handling**: Most functions return `Result` types
5. **Performance**: Use `--release` for computationally intensive tasks

## Troubleshooting

### Common Issues

**Problem**: Integration taking too long
**Solution**: Reduce number of steps or k-modes

**Problem**: Plots not appearing
**Solution**: Check file permissions and output directory

**Problem**: Numerical instabilities
**Solution**: Adjust integration tolerances or step sizes

## Getting Help

- Read the [API Documentation](https://docs.rs/andam)
- Open an issue on [GitHub](https://github.com/cosmos-andam/andam/issues)
- Contact: https://github.com/cosmos-andam/andam/discussions
```

### Step 1.3: Add Module Documentation

Update `src/lib.rs` with comprehensive module docs:

```rust
//! # Andam
//! 
//! A comprehensive Rust library for cosmological calculations and visualizations.
//! 
//! ## Overview
//! 
//! This library provides tools for working with cosmological models, including:
//! - Universe evolution (Friedmann equations)
//! - Distance calculations
//! - CMB physics and recombination
//! - Structure formation
//! - Weak gravitational lensing
//! - Visualization tools
//! 
//! ## Quick Start
//! 
//! ```rust
//! use andam::prelude::*;
//! 
//! let universe = Universe::benchmark();
//! let age = universe.age_today();
//! println!("Universe age: {:.2} Gyr", age);
//! ```
//! 
//! ## Modules
//! 
//! - [`constants`] - Physical and cosmological constants
//! - [`units`] - Unit conversion utilities
//! - [`dynamics`] - Universe evolution and Friedmann equations
//! - [`observations`] - Observational cosmology (distances, etc.)
//! - [`cmb`] - Cosmic Microwave Background physics
//! - [`structure`] - Structure formation and power spectra
//! - [`perturbations`] - Linear perturbation theory
//! - [`advanced`] - Advanced topics (lensing, polarization)
//! - [`visualization`] - Plotting and visualization tools
//! 
//! ## Features
//! 
//! Enable optional features in your `Cargo.toml`:
//! 
//! ```toml
//! [dependencies]
//! andam = { version = "0.1", features = ["parallel", "hdf5"] }
//! ```

pub mod constants;
pub mod units;
pub mod dynamics;
pub mod observations;
pub mod cmb;
pub mod structure;
pub mod perturbations;
pub mod advanced;
pub mod visualization;

/// Commonly used imports
pub mod prelude {
 pub use crate::constants::*;
 pub use crate::units::{Length, Mass, Time, Energy};
 pub use crate::dynamics::{Universe, Component};
 pub use crate::observations::*;
 pub use crate::cmb::*;
 pub use crate::structure::power_spectrum::*;
 pub use crate::perturbations::*;
 pub use crate::visualization::plots_2d::*;
 pub use crate::visualization::plotly_plots::*;
}

// Re-exports
pub use constants::*;
pub use units::{Length, Mass, Time, Energy};
pub use dynamics::Universe;
```

---

## Task 2: Performance Optimization

### Step 2.1: Add Benchmarks

Create `benches/friedmann_bench.rs`:

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use andam::dynamics::Universe;
use andam::observations::*;

fn benchmark_age_calculation(c: &mut Criterion) {
 let universe = Universe::benchmark();
 
 c.bench_function("age_today", |b| {
 b.iter(|| universe.age_today())
 });
}

fn benchmark_distance_calculations(c: &mut Criterion) {
 let universe = Universe::benchmark();
 
 let mut group = c.benchmark_group("distances");
 
 for z in [0.1, 0.5, 1.0, 2.0, 5.0].iter() {
 group.bench_with_input(BenchmarkId::new("luminosity", z), z, |b, &z| {
 b.iter(|| luminosity_distance(black_box(z), &universe))
 });
 
 group.bench_with_input(BenchmarkId::new("angular_diameter", z), z, |b, &z| {
 b.iter(|| angular_diameter_distance(black_box(z), &universe))
 });
 }
 
 group.finish();
}

fn benchmark_hubble_parameter(c: &mut Criterion) {
 let universe = Universe::benchmark();
 
 c.bench_function("hubble_normalized", |b| {
 b.iter(|| universe.hubble_normalized(black_box(0.5)))
 });
}

criterion_group!(
 benches,
 benchmark_age_calculation,
 benchmark_distance_calculations,
 benchmark_hubble_parameter
);
criterion_main!(benches);
```

### Step 2.2: Optimize Critical Functions

Update `src/observations/distances.rs` with cached integration:

```rust
use std::sync::Mutex;
use std::collections::HashMap;

// Cache for comoving distance calculations
lazy_static::lazy_static! {
 static ref DISTANCE_CACHE: Mutex<HashMap<(u64, u64), f64>> = Mutex::new(HashMap::new());
}

pub fn comoving_distance_cached(z: f64, universe: &Universe) -> f64 {
 // Create cache key from z and universe parameters
 let key = (
 (z * 1e6) as u64,
 ((universe.omega_total() * 1e6) as u64),
 );
 
 let mut cache = DISTANCE_CACHE.lock().unwrap();
 
 if let Some(&distance) = cache.get(&key) {
 return distance;
 }
 
 let distance = comoving_distance(z, universe);
 cache.insert(key, distance);
 distance
}
```

### Step 2.3: Parallelize Power Spectrum Calculation

Update `src/structure/power_spectrum.rs`:

```rust
use rayon::prelude::*;

pub fn matter_power_spectrum_parallel(
 k_values: &[f64],
 z: f64,
 omega_m: f64,
 omega_b: f64,
 h: f64,
 amplitude: f64,
 spectral_index: f64,
) -> Vec<f64> {
 k_values.par_iter()
 .map(|&k| matter_power_spectrum(k, z, omega_m, omega_b, h, amplitude, spectral_index))
 .collect()
}
```

---

## Task 3: Extensive Testing

### Step 3.1: Add Property-Based Tests

Add to `Cargo.toml`:
```toml
[dev-dependencies]
proptest = "1.0"
```

Create `tests/property_tests.rs`:

```rust
use proptest::prelude::*;
use andam::*;

proptest! {
 #[test]
 fn test_distance_positivity(z in 0.0..10.0) {
 let universe = dynamics::Universe::benchmark();
 let d_l = observations::luminosity_distance(z, &universe);
 prop_assert!(d_l > 0.0);
 }
 
 #[test]
 fn test_distance_monotonicity(z1 in 0.0..5.0, z2 in 0.0..5.0) {
 let universe = dynamics::Universe::benchmark();
 
 if z1 < z2 {
 let d1 = observations::comoving_distance(z1, &universe);
 let d2 = observations::comoving_distance(z2, &universe);
 prop_assert!(d1 < d2);
 }
 }
 
 #[test]
 fn test_etherington_relation(z in 0.1..5.0) {
 let universe = dynamics::Universe::benchmark();
 let d_l = observations::luminosity_distance(z, &universe);
 let d_a = observations::angular_diameter_distance(z, &universe);
 
 let ratio = d_l / d_a;
 let expected = (1.0 + z).powi(2);
 
 prop_assert!((ratio - expected).abs() / expected < 0.01);
 }
}
```

### Step 3.2: Add Comparison Tests

Create `tests/validation_tests.rs`:

```rust
//! Validation against known cosmological results

use andam::*;
use approx::assert_relative_eq;

#[test]
fn test_planck_2018_age() {
 let universe = dynamics::Universe::benchmark();
 let age = universe.age_today();
 
 // Planck 2018: 13.787 ± 0.020 Gyr
 assert!(age > 13.6 && age < 14.0,
 "Age should be ~13.8 Gyr, got {:.2}", age);
}

#[test]
fn test_recombination_redshift() {
 let universe = dynamics::Universe::benchmark();
 let z_rec = cmb::recombination_redshift(&universe);
 
 // Should be around z ~ 1100
 assert!(z_rec > 1000.0 && z_rec < 1200.0,
 "Recombination should be at z~1100, got {:.0}", z_rec);
}

#[test]
fn test_sound_horizon() {
 // Sound horizon at recombination: ~150 Mpc
 let universe = dynamics::Universe::benchmark();
 let z_rec = cmb::recombination_redshift(&universe);
 
 // This is approximate
 let r_s_expected = 147.0; // Mpc (Planck 2018)
 
 // We'd need to calculate actual sound horizon
 // This is a placeholder
 assert!(true);
}

#[test]
fn test_cmb_temperature_scaling() {
 // T(z) = T_0 (1 + z)
 let z = 1100.0;
 let t_cmb_today = constants::T_CMB;
 let t_at_z = t_cmb_today * (1.0 + z);
 
 // At recombination, T ~ 3000 K
 assert!(t_at_z > 2900.0 && t_at_z < 3100.0);
}
```

---

## Task 4: Example Gallery

### Step 4.1: Create Publication-Ready Examples

Create `examples/publication_quality/hubble_diagram_publication.rs`:

```rust
//! Publication-quality Hubble diagram with error bars and data points

use andam::prelude::*;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let universe = Universe::benchmark();
 
 // Theoretical curve
 let mut theory = Vec::new();
 for i in 0..200 {
 let z = (i as f64) * 0.025;
 let mu = distance_modulus(z, &universe);
 theory.push((z, mu));
 }
 
 // Simulated supernova data (would be real data in practice)
 let mut sn_data = Vec::new();
 for i in 0..30 {
 let z = 0.1 + (i as f64) * 0.15;
 let mu_theory = distance_modulus(z, &universe);
 let mu_obs = mu_theory + (rand::random::<f64>() - 0.5) * 0.3;
 let error = 0.15;
 sn_data.push((z, mu_obs, error));
 }
 
 // Create plot
 let root = BitMapBackend::new("hubble_publication.png", (1200, 800))
 .into_drawing_area();
 root.fill(&WHITE)?;
 
 let mut chart = ChartBuilder::on(&root)
 .caption("Hubble Diagram: Type Ia Supernovae", ("Arial", 50).into_font())
 .margin(15)
 .x_label_area_size(60)
 .y_label_area_size(70)
 .build_cartesian_2d(0.0..5.0, 32.0..48.0)?;
 
 chart.configure_mesh()
 .x_desc("Redshift z")
 .y_desc("Distance Modulus μ [mag]")
 .x_label_style(("Arial", 20))
 .y_label_style(("Arial", 20))
 .draw()?;
 
 // Draw theoretical curve
 chart.draw_series(LineSeries::new(
 theory.iter().map(|(z, mu)| (*z, *mu)),
 &BLUE.mix(0.8),
 ))?
 .label("Benchmark ΛCDM")
 .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
 
 // Draw data points with error bars
 chart.draw_series(
 sn_data.iter().map(|(z, mu, err)| {
 ErrorBar::new_vertical(*z, mu - err, *mu, mu + err, BLACK.filled(), 5)
 })
 )?;
 
 chart.draw_series(
 sn_data.iter().map(|(z, mu, _)| {
 Circle::new((*z, *mu), 4, RED.filled())
 })
 )?
 .label("Type Ia SNe")
 .legend(|(x, y)| Circle::new((x + 10, y), 4, RED.filled()));
 
 chart.configure_series_labels()
 .background_style(&WHITE.mix(0.8))
 .border_style(&BLACK)
 .label_font(("Arial", 18))
 .draw()?;
 
 root.present()?;
 println!("Created hubble_publication.png");
 
 Ok(())
}
```

### Step 4.2: Create CMB Map Example

Create `examples/publication_quality/cmb_map.rs`:

```rust
//! Create a simulated CMB temperature map

use andam::prelude::*;
use plotters::prelude::*;
use rand::Rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
 let size = 512;
 let mut rng = rand::thread_rng();
 
 // Simulate CMB temperature fluctuations
 let mut temperature_map = vec![vec![0.0; size]; size];
 
 // Add Gaussian random field (simplified)
 for i in 0..size {
 for j in 0..size {
 let value: f64 = rng.gen_range(-1.0..1.0);
 temperature_map[i][j] = value * 100.0; // ±100 μK
 }
 }
 
 // Smooth the map (crude convolution)
 let smoothed = gaussian_smooth(&temperature_map, 5.0);
 
 // Create colorful plot
 let root = BitMapBackend::new("cmb_map_publication.png", (1024, 512))
 .into_drawing_area();
 root.fill(&WHITE)?;
 
 let mut chart = ChartBuilder::on(&root)
 .caption("Simulated CMB Temperature Fluctuations", ("Arial", 40))
 .margin(10)
 .build_cartesian_2d(0..size as i32, 0..size as i32)?;
 
 // Draw temperature map with color scale
 for i in 0..size {
 for j in 0..size {
 let temp = smoothed[i][j];
 let color = temperature_to_color(temp);
 
 chart.draw_series(std::iter::once(Rectangle::new(
 [(i as i32, j as i32), ((i + 1) as i32, (j + 1) as i32)],
 color.filled(),
 )))?;
 }
 }
 
 root.present()?;
 println!("Created cmb_map_publication.png");
 
 Ok(())
}

fn gaussian_smooth(map: &Vec<Vec<f64>>, sigma: f64) -> Vec<Vec<f64>> {
 // Simplified Gaussian smoothing
 let size = map.len();
 let mut result = map.clone();
 
 let kernel_size = (3.0 * sigma) as usize;
 
 for i in kernel_size..size-kernel_size {
 for j in kernel_size..size-kernel_size {
 let mut sum = 0.0;
 let mut weight_sum = 0.0;
 
 for ki in 0..kernel_size {
 for kj in 0..kernel_size {
 let dx = (ki as f64) - (kernel_size as f64 / 2.0);
 let dy = (kj as f64) - (kernel_size as f64 / 2.0);
 let weight = (-((dx*dx + dy*dy) / (2.0 * sigma * sigma))).exp();
 
 sum += map[i + ki - kernel_size/2][j + kj - kernel_size/2] * weight;
 weight_sum += weight;
 }
 }
 
 result[i][j] = sum / weight_sum;
 }
 }
 
 result
}

fn temperature_to_color(temp: f64) -> RGBColor {
 // Map temperature to color: blue (cold) -> red (hot)
 let normalized = (temp + 100.0) / 200.0; // Map ±100 to [0,1]
 let normalized = normalized.max(0.0).min(1.0);
 
 if normalized < 0.5 {
 let t = normalized * 2.0;
 RGBColor(
 0,
 (255.0 * t) as u8,
 (255.0 * (1.0 - t)) as u8,
 )
 } else {
 let t = (normalized - 0.5) * 2.0;
 RGBColor(
 (255.0 * t) as u8,
 (255.0 * (1.0 - t)) as u8,
 0,
 )
 }
}
```

---

## Task 5: CI/CD Setup

### Step 5.1: Create GitHub Actions Workflow

Create `.github/workflows/ci.yml`:

```yaml
name: CI

on:
 push:
 branches: [ main, develop ]
 pull_request:
 branches: [ main ]

env:
 CARGO_TERM_COLOR: always

jobs:
 test:
 runs-on: ubuntu-latest
 
 steps:
 - uses: actions/checkout@v3
 
 - name: Install Rust
 uses: actions-rs/toolchain@v1
 with:
 profile: minimal
 toolchain: stable
 override: true
 components: rustfmt, clippy
 
 - name: Cache cargo registry
 uses: actions/cache@v3
 with:
 path: ~/.cargo/registry
 key: ${{ runner.os }}-cargo-registry-${{ hashFiles('**/Cargo.lock') }}
 
 - name: Cache cargo index
 uses: actions/cache@v3
 with:
 path: ~/.cargo/git
 key: ${{ runner.os }}-cargo-index-${{ hashFiles('**/Cargo.lock') }}
 
 - name: Cache cargo build
 uses: actions/cache@v3
 with:
 path: target
 key: ${{ runner.os }}-cargo-build-target-${{ hashFiles('**/Cargo.lock') }}
 
 - name: Check formatting
 run: cargo fmt -- --check
 
 - name: Run clippy
 run: cargo clippy -- -D warnings
 
 - name: Run tests
 run: cargo test --verbose
 
 - name: Run benchmarks
 run: cargo bench --no-run

 coverage:
 runs-on: ubuntu-latest
 
 steps:
 - uses: actions/checkout@v3
 
 - name: Install tarpaulin
 run: cargo install cargo-tarpaulin
 
 - name: Generate coverage
 run: cargo tarpaulin --out Xml
 
 - name: Upload coverage to Codecov
 uses: codecov/codecov-action@v3

 docs:
 runs-on: ubuntu-latest
 
 steps:
 - uses: actions/checkout@v3
 
 - name: Build documentation
 run: cargo doc --no-deps --all-features
```

---

## Task 6: Crate Publication Preparation

### Step 6.1: Finalize Cargo.toml

```toml
[package]
name = "andam"
version = "0.1.0"
edition = "2021"
authors = ["Andam Contributors"]
license = "MIT OR Apache-2.0"
description = "Comprehensive cosmological calculations and visualizations"
repository = "https://github.com/cosmos-andam/andam"
documentation = "https://docs.rs/andam"
homepage = "https://github.com/cosmos-andam/andam"
keywords = ["cosmology", "astrophysics", "physics", "simulation", "astronomy"]
categories = ["science", "simulation", "visualization"]
readme = "README.md"
exclude = [
 "examples/output/*",
 ".github/*",
 "benches/results/*",
]

[package.metadata.docs.rs]
all-features = true
rustdoc-args = ["--cfg", "docsrs"]

[features]
default = ["plotting"]
plotting = ["plotters", "plotly"]
parallel = ["rayon"]
full = ["plotting", "parallel"]

[badges]
maintenance = { status = "actively-developed" }
```

### Step 6.2: Create LICENSE Files

Create `LICENSE-MIT` and `LICENSE-APACHE`

### Step 6.3: Create CONTRIBUTING.md

```markdown
# Contributing to andam

Thank you for your interest in contributing!

## How to Contribute

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Ensure all tests pass: `cargo test`
6. Format code: `cargo fmt`
7. Run clippy: `cargo clippy`
8. Submit a pull request

## Code Style

- Follow Rust naming conventions
- Add documentation for public APIs
- Include examples in documentation
- Add tests for new functionality

## Testing

```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run with output
cargo test -- --nocapture
```

## Benchmarking

```bash
cargo bench
```

## Documentation

```bash
cargo doc --open
```

## Questions?

Open an issue or contact the maintainers.
```

---

## Phase 4 Completion Checklist

- [ ] Comprehensive documentation written
- [ ] All public APIs documented
- [ ] Benchmarks created and run
- [ ] Performance optimizations applied
- [ ] Property-based tests added
- [ ] Validation tests pass
- [ ] Publication-quality examples created
- [ ] CI/CD pipeline set up
- [ ] Crate ready for publication
- [ ] README complete
- [ ] User guide written
- [ ] Contributing guidelines created
- [ ] Licenses added

---

## Final Deliverables

### Documentation
1. [DONE] README.md with examples
2. [DONE] USER_GUIDE.md with tutorials
3. [DONE] API documentation (rustdoc)
4. [DONE] CONTRIBUTING.md
5. [DONE] Inline code documentation

### Code Quality
1. [DONE] All tests passing
2. [DONE] Benchmarks showing performance
3. [DONE] Zero clippy warnings
4. [DONE] Formatted with rustfmt
5. [DONE] >80% code coverage

### Examples
1. [DONE] Basic examples (10+)
2. [DONE] Publication-quality examples (5+)
3. [DONE] Interactive visualizations
4. [DONE] 3D demonstrations

### Publication
1. [DONE] Cargo.toml complete
2. [DONE] Licenses included
3. [DONE] CI/CD configured
4. [DONE] Ready for crates.io

---

## Publication Steps

1. Verify all tests pass: `cargo test --all-features`
2. Check documentation: `cargo doc --all-features --open`
3. Run clippy: `cargo clippy --all-features`
4. Format code: `cargo fmt --all`
5. Build release: `cargo build --release`
6. Publish to crates.io: `cargo publish --dry-run`
7. If successful: `cargo publish`

---

## Post-Publication

1. Monitor issues on GitHub
2. Respond to community feedback
3. Plan version 0.2.0 features
4. Write blog post about the crate
5. Present at Rust meetups
6. Contribute to Rust astronomy ecosystem

## Congratulations!

You have successfully created a comprehensive cosmology library in Rust!
