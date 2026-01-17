# Andam Quick Start Guide

Get up and running with Andam cosmology calculations in 5 minutes.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
andam = "0.1.0"
```

Or clone from source:

```bash
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
cargo test
```

## Your First Program

Create `src/main.rs`:

```rust
use andam::prelude::*;

fn main() {
    // Create standard ΛCDM universe (Planck 2018)
    let universe = Universe::benchmark();

    // Calculate age of the universe
    let age = universe.age_today();
    println!("Universe age: {:.2} Gyr", age);

    // Calculate distance to redshift z=1
    let z = 1.0;
    let distance = luminosity_distance(z, &universe);
    println!("Distance to z={}: {:.1} Mpc", z, distance);
}
```

Run it:

```bash
cargo run
```

Output:
```
Universe age: 13.80 Gyr
Distance to z=1: 6601.8 Mpc
```

## Common Tasks

### 1. Calculate Universe Age at Different Times

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    println!("Universe age at different epochs:");
    println!("  Today (z=0):        {:.2} Gyr", universe.age(1.0));
    println!("  z=1:                {:.2} Gyr", universe.age(0.5));
    println!("  Recombination (z≈1100): {:.2} Myr",
             universe.age(1.0/1091.0) * 1000.0);
}
```

### 2. Plot Hubble Diagram

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Generate redshift-distance data
    let data: Vec<(f64, f64)> = (1..=200)
        .map(|i| {
            let z = i as f64 / 100.0;
            let d_l = luminosity_distance(z, &universe);
            let mu = 5.0 * (d_l * 1e6 / 10.0).log10();
            (z, mu)
        })
        .collect();

    // Create plot
    create_line_plot(
        &data,
        "hubble_diagram.png",
        "Hubble Diagram",
        "Redshift z",
        "Distance Modulus μ [mag]"
    ).unwrap();

    println!("Created hubble_diagram.png");
}
```

### 3. Calculate CMB Properties

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    // Find recombination redshift
    let z_rec = recombination_redshift(&universe);
    println!("Recombination redshift: {:.0}", z_rec);

    // Calculate angular diameter distance to CMB
    let d_a = angular_diameter_distance(z_rec, &universe);
    println!("Angular diameter distance: {:.1} Mpc", d_a);

    // Sound horizon angle
    let r_s = 145.0; // Mpc
    let theta_s = r_s / d_a;
    println!("Sound horizon angle: {:.4} radians", theta_s);
}
```

### 4. Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;

fn main() {
    // Planck 2018 parameters
    let omega_m = 0.3111;
    let omega_b = 0.049;
    let h = 0.6766;
    let a_s = 2.1e-9;
    let n_s = 0.9665;

    println!("Matter power spectrum:");
    for k in [0.01, 0.1, 1.0] {
        let delta_sq = dimensionless_power(
            k, 0.0, omega_m, omega_b, h, a_s, n_s
        );
        println!("  k = {:.2} h/Mpc: Δ²(k) = {:.4}", k, delta_sq);
    }
}
```

### 5. Interactive Visualization

```rust
use andam::prelude::*;
use andam::visualization::plotly_plots::*;

fn main() {
    let universe = Universe::benchmark();

    // Calculate scale factor vs time
    let data: Vec<(f64, f64)> = (1..=1000)
        .map(|i| {
            let a = i as f64 / 1000.0;
            let t = universe.age(a);
            (a, t)
        })
        .collect();

    // Create interactive HTML plot
    create_interactive_plot(
        &data,
        "universe_evolution.html",
        "Universe Evolution",
        "Scale Factor a",
        "Age [Gyr]"
    ).unwrap();

    println!("Open universe_evolution.html in your browser");
}
```

## Running Examples

Andam includes many example programs:

```bash
# Basic examples
cargo run --example hubble_diagram
cargo run --example universe_evolution
cargo run --example distance_measures

# Advanced examples
cargo run --example cmb_power_spectrum
cargo run --example matter_power_spectrum
cargo run --example structure_growth
```

## Key Concepts

### Universe Models

```rust
// Standard ΛCDM (Planck 2018)
let universe = Universe::benchmark();

// Einstein-de Sitter (matter only)
let eds = Universe::einstein_de_sitter(70.0);

// Custom universe
let mut custom = Universe::new(70.0);
custom.add_component(Component::matter(0.3));
custom.add_component(Component::lambda(0.7));
```

### Distance Measures

Three important distances in cosmology:

```rust
let z = 1.0;
let universe = Universe::benchmark();

// Comoving distance
let d_c = comoving_distance(z, &universe);

// Luminosity distance (for standard candles)
let d_l = luminosity_distance(z, &universe);

// Angular diameter distance (for standard rulers)
let d_a = angular_diameter_distance(z, &universe);

// They're related: d_L = (1+z)² d_A
```

### Scale Factor vs Redshift

```rust
// Convert between scale factor and redshift
let z = 1.0;
let a = 1.0 / (1.0 + z);  // a = 0.5

// Or the reverse
let a = 0.5;
let z = (1.0 / a) - 1.0;  // z = 1.0
```

## Prelude Module

The prelude imports commonly used items:

```rust
use andam::prelude::*;

// This gives you:
// - Universe, Component
// - comoving_distance, luminosity_distance, angular_diameter_distance
// - ionization_fraction, recombination_redshift
```

For specialized functions, import specific modules:

```rust
use andam::structure::power_spectrum::*;
use andam::perturbations::growth::*;
use andam::cmb::fluctuations::*;
```

## Common Patterns

### Calculate Values at Multiple Redshifts

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    let redshifts = vec![0.5, 1.0, 2.0, 5.0];

    for z in redshifts {
        let d_c = comoving_distance(z, &universe);
        let age = universe.age(1.0 / (1.0 + z));
        println!("z = {}: d_c = {:.1} Mpc, age = {:.2} Gyr",
                 z, d_c, age);
    }
}
```

### Generate Plot Data

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    // Create vectors of x and y values
    let scale_factors: Vec<f64> = (10..=100)
        .map(|i| i as f64 / 100.0)
        .collect();

    let hubble_values: Vec<f64> = scale_factors
        .iter()
        .map(|&a| universe.hubble(a))
        .collect();

    // Combine into (x, y) pairs
    let data: Vec<(f64, f64)> = scale_factors
        .iter()
        .zip(hubble_values.iter())
        .map(|(&a, &h)| (a, h))
        .collect();

    // Now data can be plotted
}
```

### Error Handling

```rust
use andam::visualization::plots_2d::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let data = vec![(0.0, 0.0), (1.0, 1.0)];

    create_line_plot(
        &data,
        "plot.png",
        "Title",
        "X",
        "Y"
    )?;

    Ok(())
}
```

## Performance Tips

1. Use release mode for production:
   ```bash
   cargo run --release
   ```

2. Reuse `Universe` instances:
   ```rust
   let universe = Universe::benchmark(); // Create once
   for z in redshifts {
       let d = luminosity_distance(z, &universe); // Reuse
   }
   ```

3. For parallel calculations, use Rayon:
   ```rust
   use rayon::prelude::*;

   let distances: Vec<f64> = redshifts
       .par_iter()
       .map(|&z| luminosity_distance(z, &universe))
       .collect();
   ```

## Typical Values

For quick sanity checks:

```rust
let universe = Universe::benchmark();

// Age of universe: ~13.8 Gyr
assert!((universe.age_today() - 13.8).abs() < 0.1);

// H0: 67.66 km/s/Mpc
assert_eq!(universe.h0, 67.66);

// Comoving distance to z=1: ~3300 Mpc
let d = comoving_distance(1.0, &universe);
assert!((d - 3300.0).abs() < 100.0);

// Recombination: z ~ 1090
let z_rec = recombination_redshift(&universe);
assert!((z_rec - 1090.0).abs() < 10.0);
```

## Next Steps

1. Read the [User Guide](USER_GUIDE.md) for comprehensive documentation
2. Check [API Reference](API_REFERENCE.md) for detailed function signatures
3. Explore examples in the `examples/` directory
4. Run `cargo doc --open` for generated documentation

## Common Issues

### Plot Not Created

Check that output directory exists:
```rust
use std::fs;
fs::create_dir_all("output")?;
create_line_plot(&data, "output/plot.png", ...)?;
```

### Numerical Instabilities

Ensure inputs are in valid ranges:
- Scale factor: 0 < a ≤ 1
- Redshift: z ≥ 0
- Density parameters: Ω > 0

### Slow Performance

Always compile with --release for production code.

## Help and Support

- Documentation: `cargo doc --open`
- Examples: `cargo run --example`
- Tests: `cargo test`
- GitHub: https://github.com/cosmos-andam/andam

## Complete Example

Here's a complete program that demonstrates multiple features:

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create universe
    let universe = Universe::benchmark();

    // Print basic info
    println!("ΛCDM Universe (Planck 2018)");
    println!("  H0 = {:.2} km/s/Mpc", universe.h0);
    println!("  Age = {:.2} Gyr", universe.age_today());
    println!();

    // CMB properties
    let z_rec = recombination_redshift(&universe);
    let d_a_cmb = angular_diameter_distance(z_rec, &universe);
    println!("CMB Properties:");
    println!("  Recombination redshift: {:.0}", z_rec);
    println!("  Angular diameter distance: {:.1} Mpc", d_a_cmb);
    println!();

    // Calculate distances to various redshifts
    println!("Distances:");
    for z in [0.1, 0.5, 1.0, 2.0] {
        let d_c = comoving_distance(z, &universe);
        let d_l = luminosity_distance(z, &universe);
        let d_a = angular_diameter_distance(z, &universe);
        println!("  z = {:.1}: d_c = {:.0} Mpc, d_L = {:.0} Mpc, d_A = {:.0} Mpc",
                 z, d_c, d_l, d_a);
    }
    println!();

    // Create Hubble diagram
    let data: Vec<(f64, f64)> = (1..=200)
        .map(|i| {
            let z = i as f64 / 100.0;
            let d_l = luminosity_distance(z, &universe);
            let mu = 5.0 * (d_l * 1e6 / 10.0).log10();
            (z, mu)
        })
        .collect();

    create_line_plot(
        &data,
        "hubble_diagram.png",
        "Hubble Diagram",
        "Redshift z",
        "Distance Modulus μ [mag]"
    )?;

    println!("Created hubble_diagram.png");

    Ok(())
}
```

Run this with:
```bash
cargo run --release
```

Output:
```
ΛCDM Universe (Planck 2018)
  H0 = 67.66 km/s/Mpc
  Age = 13.80 Gyr

CMB Properties:
  Recombination redshift: 1091
  Angular diameter distance: 13.7 Mpc

Distances:
  z = 0.1: d_c = 421 Mpc, d_L = 463 Mpc, d_A = 383 Mpc
  z = 0.5: d_c = 1902 Mpc, d_L = 2854 Mpc, d_A = 1268 Mpc
  z = 1.0: d_c = 3298 Mpc, d_L = 6595 Mpc, d_A = 1649 Mpc
  z = 2.0: d_c = 5279 Mpc, d_L = 15838 Mpc, d_A = 1760 Mpc

Created hubble_diagram.png
```

Happy cosmology computing!
