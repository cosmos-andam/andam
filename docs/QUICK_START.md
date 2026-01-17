# Andam Quick Start Guide

Get up and running with Andam in 5 minutes.

## Installation

```bash
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
cargo test
```

Or add to your `Cargo.toml`:
```toml
[dependencies]
andam = "0.1.0"
```

## Your First Program

```rust
use andam::prelude::*;

fn main() {
    // Create standard ΛCDM universe
    let universe = Universe::benchmark();

    // Calculate universe age
    let age = universe.age_today();
    println!("Universe age: {:.2} Gyr", age);

    // Calculate distance to redshift z=1
    let distance = luminosity_distance(1.0, &universe);
    println!("Distance to z=1: {:.0} Mpc", distance);
}
```

Output:
```
Universe age: 13.80 Gyr
Distance to z=1: 6602 Mpc
```

## Three Essential Examples

### 1. Distance Calculations

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();
    let z = 1.0;

    let d_c = comoving_distance(z, &universe);
    let d_l = luminosity_distance(z, &universe);
    let d_a = angular_diameter_distance(z, &universe);

    println!("Comoving distance: {:.0} Mpc", d_c);
    println!("Luminosity distance: {:.0} Mpc", d_l);
    println!("Angular diameter distance: {:.0} Mpc", d_a);
}
```

### 2. Create a Plot

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

fn main() {
    let universe = Universe::benchmark();

    // Generate data
    let data: Vec<(f64, f64)> = (10..=100)
        .map(|i| {
            let a = i as f64 / 100.0;
            (a, universe.hubble(a))
        })
        .collect();

    // Create plot
    create_line_plot(
        &data,
        "hubble.png",
        "Hubble Parameter",
        "Scale Factor a",
        "H(a) [km/s/Mpc]"
    ).unwrap();
}
```

### 3. CMB Properties

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();

    let z_rec = recombination_redshift(&universe);
    println!("Recombination at z = {:.0}", z_rec);

    let x_e = ionization_fraction(1100.0, &universe);
    println!("Ionization fraction at z=1100: {:.4}", x_e);
}
```

## Running Examples

```bash
# Basic examples
cargo run --example hubble_diagram
cargo run --example universe_evolution

# Advanced examples
cargo run --example cmb_power_spectrum
cargo run --example matter_power_spectrum
```

## Key Concepts

### Universe Models

```rust
// Planck 2018 parameters
let universe = Universe::benchmark();

// Matter-only universe
let eds = Universe::einstein_de_sitter(70.0);

// Custom universe
let mut custom = Universe::new(70.0);
custom.add_component(Component::matter(0.3));
custom.add_component(Component::lambda(0.7));
```

### Scale Factor and Redshift

```rust
// Related by: a = 1/(1+z)
let z = 1.0;
let a = 1.0 / (1.0 + z);  // a = 0.5
```

### Prelude Module

The prelude imports commonly used items:

```rust
use andam::prelude::*;
// Gives you: Universe, Component, distance functions, CMB functions
```

For specialized features:
```rust
use andam::structure::power_spectrum::*;
use andam::perturbations::growth::*;
```

## Next Steps

1. Read the [User Guide](USER_GUIDE.md) for comprehensive documentation
2. Check the [API Reference](API_REFERENCE.md) for detailed function signatures
3. Explore examples in the `examples/` directory
4. Run `cargo doc --open` for generated documentation

## Common Issues

### Plot not created?
Check that the output directory exists and you have write permissions.

### Numerical instabilities?
Ensure inputs are in valid ranges:
- Scale factor: 0 < a ≤ 1
- Redshift: z ≥ 0

### Slow performance?
Always use `cargo run --release` for production code.

## Getting Help

- Documentation: `cargo doc --open`
- Examples: `examples/` directory
- Issues: https://github.com/cosmos-andam/andam/issues
- Discussions: https://github.com/cosmos-andam/andam/discussions
