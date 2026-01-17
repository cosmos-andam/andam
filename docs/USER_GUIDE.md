# Andam User Guide

Comprehensive guide to using the Andam cosmology library.

## Table of Contents

1. [Introduction](#introduction)
2. [Core Concepts](#core-concepts)
3. [Module Guide](#module-guide)
4. [Advanced Topics](#advanced-topics)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

## Introduction

Andam is a Rust library implementing cosmological calculations from Ryden's "Introduction to Cosmology" and Dodelson's "Modern Cosmology". It provides numerical methods for:

- Universe evolution (Friedmann equations)
- Distance calculations in curved spacetime
- CMB physics and recombination
- Structure formation and power spectra
- Perturbation theory
- Weak gravitational lensing

See [Quick Start](QUICK_START.md) for installation and first examples.

## Core Concepts

### Universe Models

The `Universe` struct represents a cosmological model:

```rust
use andam::dynamics::*;

// Standard ΛCDM with Planck 2018 parameters
let universe = Universe::benchmark();

// Einstein-de Sitter (Ω_m = 1)
let eds = Universe::einstein_de_sitter(70.0);

// Custom components
let mut custom = Universe::new(67.4);
custom.add_component(Component::matter(0.3));
custom.add_component(Component::radiation(9.3e-5));
custom.add_component(Component::lambda(0.7));
```

### Components and Equation of State

Components are characterized by w where p = wρc²:

```rust
Component::matter(0.3)      // w = 0,  ρ ∝ a^-3
Component::radiation(9e-5)  // w = 1/3, ρ ∝ a^-4
Component::lambda(0.7)      // w = -1,  ρ = constant
Component::new(0.1, -0.8)   // Custom w
```

### Scale Factor and Redshift

```rust
// Today: a = 1, z = 0
// CMB: a ≈ 0.001, z ≈ 1090
// General: a = 1/(1+z)

let a = 0.5;
let z = (1.0 / a) - 1.0;  // z = 1.0
```

## Module Guide

### 1. Constants

Physical and cosmological constants in SI units:

```rust
use andam::constants::*;

// Physical constants
C            // Speed of light: 2.998e8 m/s
G            // Gravitational constant
K_B          // Boltzmann constant
M_SUN        // Solar mass

// Cosmological constants (Planck 2018)
H0_PLANCK        // 67.66 km/s/Mpc
OMEGA_M_PLANCK   // 0.3111
OMEGA_B_PLANCK   // 0.0490
T_CMB            // 2.7255 K

// Derived quantities
let t_h = hubble_time();         // Hubble time
let rho_c = critical_density();  // Critical density
```

### 2. Units

Conversion between cosmological units:

```rust
use andam::units::*;

// Length
let meters = Length::Mpc.to_meters(1.0);      // Mpc to m
let ly = Length::Mpc.to_lightyears(1.0);      // Mpc to ly

// Time
let seconds = Time::Gyr.to_seconds(13.8);     // Gyr to s
let years = Time::Gyr.to_years(13.8);         // Gyr to yr

// Temperature-energy
let e_ev = temperature_to_ev(2.7255);         // K to eV
let t_k = ev_to_temperature(13.6);            // eV to K
```

### 3. Dynamics

Universe evolution via Friedmann equations:

```rust
use andam::dynamics::Universe;

let universe = Universe::benchmark();

// Hubble parameter
let h_today = universe.hubble(1.0);           // km/s/Mpc
let h_half = universe.hubble(0.5);

// Age calculation
let age = universe.age(1.0);                  // Gyr
let age_rec = universe.age(1.0/1091.0);

// Deceleration parameter
let q = universe.deceleration(1.0);

// Total density
let omega = universe.omega_total();           // Should ≈ 1.0
```

### 4. Observations - Distance Measures

Three fundamental distances in cosmology:

```rust
use andam::observations::*;

let universe = Universe::benchmark();
let z = 1.0;

// Comoving distance (line-of-sight)
let d_c = comoving_distance(z, &universe);

// Luminosity distance (standard candles)
let d_l = luminosity_distance(z, &universe);
let mu = 5.0 * (d_l * 1e6 / 10.0).log10();    // Distance modulus

// Angular diameter distance (standard rulers)
let d_a = angular_diameter_distance(z, &universe);
let theta = 1.0 / d_a;                        // Angular size of 1 Mpc

// Etherington relation: d_L = (1+z)² d_A
assert!((d_l / ((1.0 + z).powi(2) * d_a) - 1.0).abs() < 1e-10);
```

### 5. CMB Physics

#### Recombination

```rust
use andam::cmb::*;

let universe = Universe::benchmark();

// Find recombination redshift (x_e = 0.5)
let z_rec = recombination_redshift(&universe);  // ~1091

// Ionization fraction vs redshift
for z in [500.0, 800.0, 1000.0, 1200.0] {
    let x_e = ionization_fraction(z, &universe);
    println!("z={}: x_e={:.4}", z, x_e);
}
```

#### Angular Power Spectrum

```rust
use andam::cmb::fluctuations::*;

let universe = Universe::benchmark();

// Compute C_ℓ
let c_l = angular_power_spectrum(2000, &universe);

// Dimensionless power for plotting
let d_l = dimensionless_power_spectrum(&c_l);

// Peak positions
let peaks = acoustic_peak_positions(&universe);
println!("First acoustic peak: l = {}", peaks[0]);  // ~219
```

### 6. Structure Formation

#### Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;

// Planck 2018 parameters
let omega_m = 0.3111;
let omega_b = 0.049;
let h = 0.6766;
let a_s = 2.1e-9;
let n_s = 0.9665;

// Power at wavenumber k
let k = 0.1;  // h/Mpc
let p_k = matter_power_spectrum(k, 0.0, omega_m, omega_b, h, a_s, n_s);

// Dimensionless power
let delta_sq = dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s);
```

#### Transfer Function

```rust
use andam::structure::power_spectrum::*;

// Eisenstein-Hu transfer function
let t_k = transfer_function_eh(0.1, 0.3111, 0.049, 0.6766);

// Primordial power
let p_r = primordial_power_spectrum(0.05, 2.1e-9, 0.9665);
```

### 7. Perturbation Theory

#### Linear Growth Factor

```rust
use andam::perturbations::growth::*;

let universe = Universe::benchmark();

// Growth factor normalized to D(1) = 1
let d_half = growth_factor(0.5, &universe);
let d_today = growth_factor(1.0, &universe);  // = 1.0

// Growth rate f = d ln D / d ln a
let f = growth_rate(1.0, &universe);
println!("Growth rate today: {:.3}", f);
```

#### Boltzmann Solver

```rust
use andam::perturbations::*;

let universe = Universe::benchmark();

// Create solver for mode k
let mut solver = BoltzmannSolver::new(
    universe,
    0.1,                          // k in h/Mpc
    PerturbationMode::Scalar,
    10                            // l_max
);

// Evolve from radiation to matter era
let results = solver.evolve(1e-5, 1e-2);
```

### 8. Weak Lensing

#### Convergence Field

```rust
use andam::advanced::weak_lensing::*;

let universe = Universe::benchmark();

// Generate convergence field
let field = ConvergenceField::from_power_spectrum(
    &universe,
    1.0,      // source redshift
    128,      // grid size
    1.0       // field size in degrees
);
```

#### Shear Calculation

```rust
use andam::advanced::weak_lensing::*;

let field = ConvergenceField::new(128);

// Convert to shear
let (gamma1, gamma2) = convergence_to_shear(&field.kappa);

// Shear magnitude and angle
let shear = Shear { gamma1: 0.03, gamma2: 0.04 };
println!("Shear: |γ| = {:.4}", shear.magnitude());
println!("Angle: φ = {:.2} rad", shear.angle());
```

#### Lensing Power Spectrum

```rust
use andam::advanced::weak_lensing::*;

let universe = Universe::benchmark();

// Convergence power spectrum
for l in [100, 500, 1000] {
    let c_l = lensing_power_spectrum(l, 1.0, &universe);
    println!("l={}: C_l^κκ = {:.4e}", l, c_l);
}
```

### 9. Visualization

#### Static Plots (PNG)

```rust
use andam::visualization::plots_2d::*;

// Line plot
create_line_plot(&data, "output.png", "Title", "X", "Y")?;

// Log-log plot
create_loglog_plot(&data, "loglog.png", "Title", "X", "Y")?;

// Multiple series
create_multiline_plot(
    &[series1, series2],
    &["Label 1", "Label 2"],
    "multi.png",
    "Title",
    "X",
    "Y"
)?;
```

#### Interactive Plots (HTML)

```rust
use andam::visualization::plotly_plots::*;

// Interactive plot with zoom and hover
create_interactive_plot(&data, "plot.html", "Title", "X", "Y")?;

// Interactive log-log
create_loglog_interactive(&data, "log.html", "Title", "X", "Y")?;
```

## Advanced Topics

### Custom Cosmology

```rust
let mut universe = Universe::new(70.0);

// Early dark energy
let ede = Component::new(0.05, -0.9);
universe.add_component(ede);

// Massive neutrinos
let nu = Component::new(0.002, 1.0/3.0);
universe.add_component(nu);
```

### Parallel Computation

```rust
use rayon::prelude::*;

let redshifts: Vec<f64> = (1..=1000).map(|i| i as f64 / 100.0).collect();

// Parallel distance calculation
let distances: Vec<f64> = redshifts
    .par_iter()
    .map(|&z| luminosity_distance(z, &universe))
    .collect();
```

### Error Handling

```rust
use std::error::Error;

fn create_plot() -> Result<(), Box<dyn Error>> {
    let data = vec![(0.0, 0.0), (1.0, 1.0)];
    create_line_plot(&data, "plot.png", "Title", "X", "Y")?;
    Ok(())
}
```

## Best Practices

### 1. Reuse Universe Instances

```rust
// Good: create once
let universe = Universe::benchmark();
for z in redshifts {
    let d = luminosity_distance(z, &universe);
}

// Bad: creating multiple times
for z in redshifts {
    let universe = Universe::benchmark();  // Inefficient!
    let d = luminosity_distance(z, &universe);
}
```

### 2. Validate Physical Constraints

```rust
let universe = Universe::benchmark();

// Check flatness
let omega_total = universe.omega_total();
assert!((omega_total - 1.0).abs() < 0.01);

// Verify age
let age = universe.age_today();
assert!((age - 13.8).abs() < 0.2);  // Within 200 Myr
```

### 3. Use Appropriate Precision

```rust
// For plots: 100 points is sufficient
let data: Vec<(f64, f64)> = (1..=100)
    .map(|i| calculate(i as f64 / 100.0))
    .collect();

// For science: verify convergence
let n_points = vec![100, 1000, 10000];
for n in n_points {
    let result = calculate_with_n_points(n);
    println!("n={}: result={:.6}", n, result);
}
```

### 4. Document Assumptions

```rust
// Using Planck 2018 TT,TE,EE+lowE+lensing:
// H0 = 67.66 ± 0.42 km/s/Mpc
// Ω_m = 0.3111 ± 0.0056
// Ω_b h² = 0.02242 ± 0.00014
let universe = Universe::benchmark();
```

## Troubleshooting

### Numerical Instabilities

**Problem**: Functions return NaN or Inf.

**Solution**: Check input ranges:
```rust
// Valid ranges
assert!(a > 0.0 && a <= 1.0);  // Scale factor
assert!(z >= 0.0);              // Redshift
assert!(omega > 0.0);           // Density parameters
```

### Plot Files Not Created

**Problem**: No error but file doesn't exist.

**Solution**: Create output directory:
```rust
use std::fs;
fs::create_dir_all("output")?;
create_line_plot(&data, "output/plot.png", ...)?;
```

### Slow Performance

**Problem**: Calculations take too long.

**Solution**: Use release mode:
```bash
cargo build --release
cargo run --release --example matter_power_spectrum
```

### Test Failures

**Problem**: Tests fail with small tolerance errors.

**Solution**: This is normal for numerical methods. Check if error is reasonable:
```bash
cargo test -- --nocapture
# Look at actual vs expected values
```

## Performance Notes

### Typical Computation Times

On modern CPU (single core):
- Universe age: ~1 μs
- Distance (one redshift): ~10 μs
- Matter power (one k): ~5 μs
- Growth factor (one a): ~100 μs
- CMB spectrum (l_max=2000): ~5 ms

### Memory Usage

- Universe struct: ~1 KB
- CMB spectrum (l=2000): ~16 KB
- Convergence field (128×128): ~128 KB

### Optimization

1. Use `--release` for production
2. Reuse Universe instances
3. Use Rayon for parallel work
4. Avoid unnecessary allocations

## Additional Resources

- [API Reference](API_REFERENCE.md) - Detailed API documentation
- [Quick Start](QUICK_START.md) - Getting started guide
- [Project Overview](PROJECT_OVERVIEW.md) - Architecture and phases
- Examples directory - 12+ complete programs
- `cargo doc --open` - Generated documentation

## Scientific References

1. Ryden, B. (2016). Introduction to Cosmology (2nd ed.)
2. Dodelson, S. (2003). Modern Cosmology
3. Planck Collaboration (2018). Cosmological parameters
4. Eisenstein & Hu (1998). ApJ 496:605
