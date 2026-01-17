# Andam User Guide

Complete guide to using the Andam cosmology library for cosmological calculations and visualizations.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Core Concepts](#core-concepts)
4. [Module Reference](#module-reference)
5. [Examples](#examples)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)

---

## Installation

### As a Library Dependency

Add to your `Cargo.toml`:

```toml
[dependencies]
andam = "0.1.0"
```

### From Source

```bash
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
```

### Running Tests

```bash
cargo test
```

### Running Examples

```bash
# List all examples
cargo run --example

# Run specific example
cargo run --example hubble_diagram
cargo run --example universe_evolution
cargo run --example cmb_power_spectrum
```

---

## Quick Start

### Basic Universe Calculations

```rust
use andam::prelude::*;

fn main() {
    // Create a standard ΛCDM universe with Planck 2018 parameters
    let universe = Universe::benchmark();

    // Calculate the age of the universe
    let age_gyr = universe.age_today();
    println!("Universe age: {:.2} Gyr", age_gyr);

    // Calculate Hubble parameter at different scale factors
    let h_today = universe.hubble(1.0);
    let h_half = universe.hubble(0.5);
    println!("H(a=1.0) = {:.2} km/s/Mpc", h_today);
    println!("H(a=0.5) = {:.2} km/s/Mpc", h_half);
}
```

### Distance Calculations

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    // Calculate distances to redshift z=1.0
    let z = 1.0;
    let d_l = luminosity_distance(z, &universe);
    let d_a = angular_diameter_distance(z, &universe);
    let d_c = comoving_distance(z, &universe);

    println!("Distances to z = {}:", z);
    println!("  Luminosity distance: {:.1} Mpc", d_l);
    println!("  Angular diameter distance: {:.1} Mpc", d_a);
    println!("  Comoving distance: {:.1} Mpc", d_c);
}
```

### Simple Visualization

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Generate scale factors from 0.1 to 1.0
    let scale_factors: Vec<f64> = (10..=100)
        .map(|i| i as f64 / 100.0)
        .collect();

    // Calculate Hubble parameter at each scale factor
    let hubble_values: Vec<f64> = scale_factors
        .iter()
        .map(|&a| universe.hubble(a))
        .collect();

    // Create plot
    let data: Vec<(f64, f64)> = scale_factors
        .iter()
        .zip(hubble_values.iter())
        .map(|(&a, &h)| (a, h))
        .collect();

    create_line_plot(
        &data,
        "hubble_evolution.png",
        "Hubble Parameter Evolution",
        "Scale Factor a",
        "H(a) [km/s/Mpc]",
    ).expect("Failed to create plot");

    println!("Created hubble_evolution.png");
}
```

---

## Core Concepts

### Universe Models

The `Universe` struct represents a cosmological model with various components:

```rust
use andam::dynamics::*;

// Planck 2018 ΛCDM parameters
let universe = Universe::benchmark();

// Einstein-de Sitter (matter only)
let eds = Universe::einstein_de_sitter(70.0); // H0 = 70 km/s/Mpc

// Custom universe
let mut custom = Universe::new(67.4); // H0 = 67.4 km/s/Mpc
custom.add_component(Component::matter(0.3)); // Ω_m = 0.3
custom.add_component(Component::lambda(0.7)); // Ω_Λ = 0.7
```

### Components

Universe components are characterized by their equation of state parameter w:

- **Matter**: w = 0, ρ ∝ a^(-3)
- **Radiation**: w = 1/3, ρ ∝ a^(-4)
- **Dark Energy (Λ)**: w = -1, ρ = constant

```rust
use andam::dynamics::Component;

let matter = Component::matter(0.3111);     // Ω_m = 0.3111
let radiation = Component::radiation(9.3e-5); // Ω_r = 9.3e-5
let lambda = Component::lambda(0.6889);     // Ω_Λ = 0.6889

// Custom equation of state
let custom = Component::new(0.05, -0.8); // Ω = 0.05, w = -0.8
```

### Scale Factor and Redshift

The scale factor a and redshift z are related by:

```
a = 1 / (1 + z)
z = (1 / a) - 1
```

- Today: a = 1.0, z = 0
- Recombination: a ≈ 0.0009, z ≈ 1100
- High redshift: a → 0, z → ∞

---

## Module Reference

### 1. Constants and Units

#### Physical Constants

```rust
use andam::constants::*;

println!("Speed of light: {} m/s", C);
println!("Hubble constant: {} km/s/Mpc", H0_PLANCK);
println!("CMB temperature: {} K", T_CMB);

// Helper functions
let t_h = hubble_time(); // Hubble time in seconds
let rho_c = critical_density(); // Critical density in kg/m^3
```

#### Unit Conversions

```rust
use andam::units::*;

// Length conversions
let meters = Length::Mpc.to_meters(1.0);
let km = Length::Mpc.to_km(1.0);
let ly = Length::Mpc.to_lightyears(1.0);

// Time conversions
let seconds = Time::Gyr.to_seconds(13.8);
let years = Time::Gyr.to_years(13.8);

// Temperature to energy
let energy_ev = temperature_to_ev(2.725); // CMB temperature
```

### 2. Fundamental Observations

#### Hubble's Law

```rust
use andam::fundamental::*;

let universe = Universe::benchmark();

// Hubble parameter at any scale factor
let h_today = hubble_parameter(1.0, &universe);
let h_early = hubble_parameter(0.1, &universe);

// Deceleration parameter
let q = deceleration_parameter(1.0, &universe);
println!("Deceleration parameter today: {:.3}", q);
```

### 3. Dynamics and Evolution

#### Age Calculation

```rust
use andam::dynamics::Universe;

let universe = Universe::benchmark();

// Age at different scale factors
let age_today = universe.age(1.0);  // Age today
let age_rec = universe.age(1.0 / 1091.0); // Age at recombination

println!("Age today: {:.2} Gyr", age_today);
println!("Age at recombination: {:.3} Myr", age_rec * 1000.0);
```

#### Universe Evolution

```rust
use andam::dynamics::Universe;

let universe = Universe::benchmark();

// Evolution from early times to today
for i in 1..=10 {
    let a = i as f64 / 10.0;
    let t = universe.age(a);
    let h = universe.hubble(a);

    println!("a = {:.1}: t = {:.2} Gyr, H = {:.1} km/s/Mpc", a, t, h);
}
```

### 4. Distance Measures

#### Comoving Distance

```rust
use andam::observations::*;

let universe = Universe::benchmark();

// Comoving distance to various redshifts
for z in [0.5, 1.0, 2.0, 5.0] {
    let d_c = comoving_distance(z, &universe);
    println!("z = {}: d_c = {:.1} Mpc", z, d_c);
}
```

#### Luminosity Distance

```rust
use andam::observations::*;

let universe = Universe::benchmark();

// Used for standard candles (supernovae)
let z = 1.0;
let d_l = luminosity_distance(z, &universe);
println!("Luminosity distance to z={}: {:.1} Mpc", z, d_l);

// Distance modulus
let mu = 5.0 * (d_l * 1e6 / 10.0).log10();
println!("Distance modulus: {:.2} mag", mu);
```

#### Angular Diameter Distance

```rust
use andam::observations::*;

let universe = Universe::benchmark();

// Used for standard rulers (BAO, CMB)
let z = 1100.0; // CMB
let d_a = angular_diameter_distance(z, &universe);
println!("Angular diameter distance to CMB: {:.1} Mpc", d_a);

// Angular size of 1 Mpc at this redshift
let theta_rad = 1.0 / d_a;
let theta_deg = theta_rad * 180.0 / std::f64::consts::PI;
println!("1 Mpc subtends {:.4} degrees", theta_deg);
```

### 5. CMB Physics

#### Recombination

```rust
use andam::cmb::*;

let universe = Universe::benchmark();

// Recombination redshift
let z_rec = recombination_redshift(&universe);
println!("Recombination redshift: {:.0}", z_rec);

// Ionization fraction at various redshifts
for z in [500.0, 800.0, 1000.0, 1200.0, 1500.0] {
    let x_e = ionization_fraction(z, &universe);
    println!("z = {}: x_e = {:.4}", z, x_e);
}
```

#### CMB Temperature Fluctuations

```rust
use andam::cmb::fluctuations::*;

let universe = Universe::benchmark();

// Compute angular power spectrum
let l_max = 2000;
let c_l = angular_power_spectrum(l_max, &universe);

// Get dimensionless power for plotting
let d_l = dimensionless_power_spectrum(&c_l);

// Find peak positions
let peaks = acoustic_peak_positions(&universe);
println!("Acoustic peak positions: {:?}", peaks);
```

### 6. Structure Formation

#### Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;

// Planck 2018 cosmological parameters
let omega_m = 0.3111;
let omega_b = 0.049;
let h = 0.6766;
let a_s = 2.1e-9;
let n_s = 0.9665;

// Power spectrum at various scales
let k_values = vec![0.01, 0.1, 1.0, 10.0]; // h/Mpc

for k in k_values {
    let p_k = matter_power_spectrum(k, 0.0, omega_m, omega_b, h, a_s, n_s);
    let delta_sq = dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s);

    println!("k = {:.2} h/Mpc: P(k) = {:.2e}, Δ²(k) = {:.4}",
             k, p_k, delta_sq);
}
```

#### Growth Factor

```rust
use andam::perturbations::growth::*;

let universe = Universe::benchmark();

// Growth factor evolution
for i in 1..=10 {
    let a = i as f64 / 10.0;
    let d = growth_factor(a, &universe);
    let f = growth_rate(a, &universe);

    println!("a = {:.1}: D = {:.4}, f = {:.3}", a, d, f);
}
```

### 7. Perturbation Theory

#### Boltzmann Solver

```rust
use andam::perturbations::*;

let universe = Universe::benchmark();
let k = 0.1; // Wavenumber in h/Mpc
let l_max = 10; // Maximum multipole

// Create solver
let mut solver = BoltzmannSolver::new(
    universe,
    k,
    PerturbationMode::Scalar,
    l_max
);

// Evolve from radiation to matter domination
let a_initial = 1e-5;
let a_final = 1e-2;
let results = solver.evolve(a_initial, a_final);

println!("Computed {} timesteps", results.len());
```

### 8. Weak Lensing

#### Convergence Field

```rust
use andam::advanced::weak_lensing::*;

let universe = Universe::benchmark();
let source_z = 1.0;
let size = 128;
let field_size_deg = 1.0;

// Generate convergence field
let field = ConvergenceField::from_power_spectrum(
    &universe,
    source_z,
    size,
    field_size_deg
);

println!("Convergence field: {} x {}", field.size, field.size);
```

#### Shear Calculation

```rust
use andam::advanced::weak_lensing::*;

let field = ConvergenceField::new(128);

// Convert convergence to shear
let (gamma1, gamma2) = convergence_to_shear(&field.kappa);

// Compute shear at a point
let shear = Shear {
    gamma1: gamma1[[64, 64]],
    gamma2: gamma2[[64, 64]],
};

println!("Shear magnitude: {:.4}", shear.magnitude());
println!("Shear angle: {:.2} rad", shear.angle());
```

#### Lensing Power Spectrum

```rust
use andam::advanced::weak_lensing::*;

let universe = Universe::benchmark();
let source_z = 1.0;

// Lensing power at various multipoles
for l in [10, 100, 500, 1000, 2000] {
    let c_l = lensing_power_spectrum(l, source_z, &universe);
    println!("l = {}: C_l^κκ = {:.4e}", l, c_l);
}
```

### 9. Visualization

#### Static 2D Plots

```rust
use andam::visualization::plots_2d::*;

// Simple line plot
let data: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
create_line_plot(
    &data,
    "output.png",
    "Plot Title",
    "X Axis",
    "Y Axis"
).unwrap();

// Log-log plot
let log_data: Vec<(f64, f64)> = vec![(1.0, 1.0), (10.0, 100.0), (100.0, 10000.0)];
create_loglog_plot(
    &log_data,
    "loglog.png",
    "Log-Log Plot",
    "X (log)",
    "Y (log)"
).unwrap();

// Multiple series
let series1: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 1.0)];
let series2: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 2.0)];
create_multiline_plot(
    &[series1, series2],
    &["Linear", "Quadratic"],
    "multiline.png",
    "Multi-Series Plot",
    "X",
    "Y"
).unwrap();
```

#### Interactive HTML Plots

```rust
use andam::visualization::plotly_plots::*;

// Interactive plot with hover
let data: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
create_interactive_plot(
    &data,
    "interactive.html",
    "Interactive Plot",
    "X Axis",
    "Y Axis"
).unwrap();

// Interactive log-log plot
let log_data: Vec<(f64, f64)> = vec![(1.0, 1.0), (10.0, 100.0)];
create_loglog_interactive(
    &log_data,
    "interactive_log.html",
    "Interactive Log Plot",
    "X (log)",
    "Y (log)"
).unwrap();
```

---

## Examples

### Example 1: Hubble Diagram

Calculate luminosity distance for Type Ia supernovae:

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Redshift range for SNe Ia
    let redshifts: Vec<f64> = (1..=200)
        .map(|i| i as f64 / 100.0)
        .collect();

    // Calculate distance modulus
    let distance_moduli: Vec<f64> = redshifts
        .iter()
        .map(|&z| {
            let d_l = luminosity_distance(z, &universe);
            5.0 * (d_l * 1e6 / 10.0).log10()
        })
        .collect();

    // Create plot data
    let data: Vec<(f64, f64)> = redshifts
        .iter()
        .zip(distance_moduli.iter())
        .map(|(&z, &mu)| (z, mu))
        .collect();

    create_line_plot(
        &data,
        "hubble_diagram.png",
        "Hubble Diagram",
        "Redshift z",
        "Distance Modulus μ [mag]"
    ).expect("Failed to create plot");

    println!("Created hubble_diagram.png");
}
```

### Example 2: Universe Evolution

Plot the evolution of the universe from Big Bang to today:

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Scale factors from 0.001 to 1.0
    let scale_factors: Vec<f64> = (1..=1000)
        .map(|i| i as f64 / 1000.0)
        .collect();

    // Calculate age at each scale factor
    let ages: Vec<f64> = scale_factors
        .iter()
        .map(|&a| universe.age(a))
        .collect();

    // Convert to Gyr and create data
    let data: Vec<(f64, f64)> = scale_factors
        .iter()
        .zip(ages.iter())
        .map(|(&a, &t)| (a, t))
        .collect();

    create_line_plot(
        &data,
        "universe_evolution.png",
        "Universe Evolution",
        "Scale Factor a",
        "Age [Gyr]"
    ).expect("Failed to create plot");

    println!("Universe age today: {:.2} Gyr", universe.age_today());
}
```

### Example 3: CMB Recombination

Calculate ionization fraction during recombination:

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Redshift range around recombination
    let redshifts: Vec<f64> = (500..=2000)
        .step_by(10)
        .map(|z| z as f64)
        .collect();

    // Calculate ionization fraction
    let ionization: Vec<f64> = redshifts
        .iter()
        .map(|&z| ionization_fraction(z, &universe))
        .collect();

    let data: Vec<(f64, f64)> = redshifts
        .iter()
        .zip(ionization.iter())
        .map(|(&z, &x_e)| (z, x_e))
        .collect();

    create_line_plot(
        &data,
        "recombination.png",
        "CMB Recombination",
        "Redshift z",
        "Ionization Fraction x_e"
    ).expect("Failed to create plot");

    let z_rec = recombination_redshift(&universe);
    println!("Recombination occurs at z = {:.0}", z_rec);
}
```

### Example 4: Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;
use andam::visualization::plotly_plots::*;

fn main() {
    let omega_m = 0.3111;
    let omega_b = 0.049;
    let h = 0.6766;
    let a_s = 2.1e-9;
    let n_s = 0.9665;

    // Wavenumber range
    let k_values: Vec<f64> = (0..=300)
        .map(|i| 10_f64.powf(-3.0 + i as f64 * 4.0 / 300.0))
        .collect();

    // Calculate dimensionless power
    let power: Vec<f64> = k_values
        .iter()
        .map(|&k| dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s))
        .collect();

    let data: Vec<(f64, f64)> = k_values
        .iter()
        .zip(power.iter())
        .map(|(&k, &p)| (k, p))
        .collect();

    create_loglog_interactive(
        &data,
        "matter_power_spectrum.html",
        "Matter Power Spectrum",
        "k [h/Mpc]",
        "Δ²(k)"
    ).expect("Failed to create plot");

    println!("Created matter_power_spectrum.html");
}
```

### Example 5: CMB Angular Power Spectrum

```rust
use andam::prelude::*;
use andam::cmb::fluctuations::*;
use andam::visualization::plotly_plots::*;

fn main() {
    let universe = Universe::benchmark();

    // Compute C_l up to l=2000
    let l_max = 2000;
    let c_l = angular_power_spectrum(l_max, &universe);

    // Get dimensionless power
    let d_l_data = dimensionless_power_spectrum(&c_l);

    // Filter out l < 2
    let plot_data: Vec<(f64, f64)> = d_l_data
        .into_iter()
        .filter(|(l, _)| *l >= 2)
        .map(|(l, d)| (l as f64, d))
        .collect();

    create_interactive_plot(
        &plot_data,
        "cmb_power_spectrum.html",
        "CMB Angular Power Spectrum",
        "Multipole l",
        "l(l+1)C_l / 2π [μK²]"
    ).expect("Failed to create plot");

    let peaks = acoustic_peak_positions(&universe);
    println!("Acoustic peaks at l = {:?}", peaks);
}
```

### Example 6: Growth Factor Evolution

```rust
use andam::prelude::*;
use andam::perturbations::growth::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Scale factors from early matter domination to today
    let scale_factors: Vec<f64> = (1..=100)
        .map(|i| i as f64 / 100.0)
        .collect();

    // Calculate growth factor
    let growth: Vec<f64> = scale_factors
        .iter()
        .map(|&a| growth_factor(a, &universe))
        .collect();

    let data: Vec<(f64, f64)> = scale_factors
        .iter()
        .zip(growth.iter())
        .map(|(&a, &d)| (a, d))
        .collect();

    create_line_plot(
        &data,
        "growth_factor.png",
        "Linear Growth Factor",
        "Scale Factor a",
        "Growth Factor D(a)"
    ).expect("Failed to create plot");

    println!("Created growth_factor.png");
}
```

---

## Advanced Usage

### Custom Cosmology

Create a universe with custom parameters:

```rust
use andam::dynamics::*;

// Create empty universe with H0 = 70 km/s/Mpc
let mut universe = Universe::new(70.0);

// Add components
universe.add_component(Component::matter(0.25));
universe.add_component(Component::radiation(0.0001));
universe.add_component(Component::lambda(0.75));

// Custom dark energy with w = -0.9
let dark_energy = Component::new(0.75, -0.9);
```

### Parallel Computation

Use Rayon for parallel calculations:

```rust
use rayon::prelude::*;
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    // Parallel distance calculation
    let redshifts: Vec<f64> = (1..=1000).map(|i| i as f64 / 100.0).collect();

    let distances: Vec<f64> = redshifts
        .par_iter()
        .map(|&z| luminosity_distance(z, &universe))
        .collect();

    println!("Computed {} distances in parallel", distances.len());
}
```

### Custom Integration

For higher precision, adjust integration parameters:

```rust
use andam::observations::*;

// The comoving_distance function uses Simpson's rule
// with default 1000 steps. For higher precision,
// you may need to modify the source or use multiple
// calculations to verify convergence.
```

### Error Handling

```rust
use andam::visualization::plots_2d::*;

fn create_plot_safely() -> Result<(), Box<dyn std::error::Error>> {
    let data: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 1.0)];

    create_line_plot(
        &data,
        "output.png",
        "Title",
        "X",
        "Y"
    )?;

    Ok(())
}

fn main() {
    match create_plot_safely() {
        Ok(_) => println!("Plot created successfully"),
        Err(e) => eprintln!("Error creating plot: {}", e),
    }
}
```

---

## Troubleshooting

### Common Issues

#### 1. Plot Files Not Created

**Problem**: Plot functions don't show errors but files aren't created.

**Solution**: Check file permissions and ensure output directory exists:

```rust
use std::fs;

// Create output directory if it doesn't exist
fs::create_dir_all("output").expect("Failed to create directory");

create_line_plot(
    &data,
    "output/plot.png",
    "Title",
    "X",
    "Y"
).expect("Failed to create plot");
```

#### 2. Numerical Instabilities

**Problem**: Functions return NaN or infinite values.

**Solution**: Check input ranges. Most functions expect:
- Scale factor: 0 < a ≤ 1
- Redshift: z ≥ 0
- Density parameters: Ω > 0

```rust
use andam::prelude::*;

let universe = Universe::benchmark();

// Bad: scale factor too small
let bad = universe.age(1e-10); // May be unstable

// Good: reasonable scale factor
let good = universe.age(0.001); // Stable
```

#### 3. Test Failures

**Problem**: Tests fail with tolerance errors.

**Solution**: This is normal for numerical methods. Check if failures are within acceptable range:

```bash
cargo test -- --nocapture
```

#### 4. Slow Performance

**Problem**: Calculations take too long.

**Solution**: Use release mode for production:

```bash
cargo build --release
cargo run --release --example matter_power_spectrum
```

### Getting Help

1. Check the API documentation:
   ```bash
   cargo doc --open
   ```

2. Review example code in `examples/` directory

3. Run tests to verify installation:
   ```bash
   cargo test --all
   ```

4. For bugs or feature requests, file an issue on GitHub

---

## Best Practices

### 1. Use Benchmark Universe

For standard calculations, always start with the benchmark universe:

```rust
let universe = Universe::benchmark(); // Planck 2018 parameters
```

### 2. Check Physical Constraints

Verify that density parameters sum to approximately 1:

```rust
let omega_total = universe.omega_total();
assert!((omega_total - 1.0).abs() < 0.01);
```

### 3. Validate Results

Compare against known values:

```rust
let universe = Universe::benchmark();
let age = universe.age_today();
assert!((age - 13.8).abs() < 0.1); // Should be ~13.8 Gyr
```

### 4. Use Appropriate Precision

For visualization, high precision isn't necessary:

```rust
// 100 points is usually sufficient for smooth plots
let scale_factors: Vec<f64> = (1..=100).map(|i| i as f64 / 100.0).collect();
```

For science, verify convergence:

```rust
// Compare results with different resolutions
let coarse = calculate_with_n_points(100);
let fine = calculate_with_n_points(1000);
let difference = (fine - coarse).abs() / fine;
assert!(difference < 0.01); // 1% convergence
```

### 5. Document Your Assumptions

Always document cosmological parameters used:

```rust
// Using Planck 2018 TT,TE,EE+lowE+lensing parameters:
// H0 = 67.66 km/s/Mpc
// Ω_m = 0.3111
// Ω_Λ = 0.6889
// Ω_b h² = 0.02242
let universe = Universe::benchmark();
```

---

## Performance Notes

### Computation Times (Approximate)

On a modern CPU (single core):

- Universe age calculation: ~1 μs
- Comoving distance (one redshift): ~10 μs
- Matter power spectrum (one k): ~5 μs
- Growth factor (one scale factor): ~100 μs
- CMB power spectrum (l_max=2000): ~5 ms
- Boltzmann evolution (k=0.1, 1000 steps): ~10 ms

### Memory Usage

Typical memory requirements:

- Universe struct: ~1 KB
- CMB power spectrum (l_max=2000): ~16 KB
- Convergence field (128x128): ~128 KB
- Boltzmann solver state: ~1 KB per timestep

### Optimization Tips

1. Reuse Universe instances:
   ```rust
   let universe = Universe::benchmark();
   for z in redshifts {
       let d = luminosity_distance(z, &universe); // Reuse universe
   }
   ```

2. Use parallel iterators for independent calculations:
   ```rust
   use rayon::prelude::*;
   let results: Vec<_> = inputs.par_iter().map(|x| compute(x)).collect();
   ```

3. Compile with optimizations:
   ```bash
   cargo build --release
   ```

---

## API Stability

### Current Version: 0.1.0

This is an initial release. API may change in future versions.

### Stability Guarantees

- Core modules (dynamics, observations) are relatively stable
- Visualization APIs may change to add features
- Advanced modules (perturbations, weak lensing) are experimental

### Deprecation Policy

Breaking changes will be announced with at least one minor version warning period.

---

## Contributing

Contributions are welcome! See CONTRIBUTING.md for guidelines.

Areas where contributions are especially valuable:

1. Additional cosmological models
2. Performance optimizations
3. More visualization options
4. Additional examples
5. Documentation improvements
6. Bug fixes

---

## References

### Textbooks

1. Ryden, B. (2016). Introduction to Cosmology (2nd ed.)
2. Dodelson, S. (2003). Modern Cosmology

### Data Sources

1. Planck Collaboration (2018). Planck 2018 results
2. Eisenstein & Hu (1998). ApJ 496:605 (Transfer function)

### Numerical Methods

1. Simpson's rule for integration
2. Euler method for ODEs
3. Limber approximation for lensing

---

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT license

at your option.

---

## Changelog

### Version 0.1.0 (2025)

Initial release with:
- Core cosmology calculations
- Distance measures
- CMB recombination and fluctuations
- Structure formation
- Perturbation theory
- Weak lensing (basic)
- 2D and interactive visualizations
- Comprehensive examples and tests
