# Andam

A comprehensive Rust library for cosmological calculations and visualizations, implementing concepts from modern cosmology textbooks.

**Andam** (அண்டம்) means "universe" in தமிழ் (Tamil).

## Features

- Complete ΛCDM cosmology calculations
- Distance measures (luminosity, angular diameter, comoving)
- CMB recombination and temperature fluctuations
- Matter power spectrum with Eisenstein-Hu transfer function
- Linear perturbation theory and growth factors
- Weak gravitational lensing calculations
- 2D static plots (PNG) and interactive visualizations (HTML)
- Planck 2018 cosmological parameters built-in
- Fast numerical integration and ODE solvers
- Comprehensive test suite with >40 tests

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
andam = "0.1.0"
```

Simple example:

```rust
use andam::prelude::*;

fn main() {
    // Create standard ΛCDM universe
    let universe = Universe::benchmark();

    // Calculate universe age
    let age = universe.age_today();
    println!("Universe age: {:.2} Gyr", age);

    // Calculate distance to redshift z=1
    let z = 1.0;
    let distance = luminosity_distance(z, &universe);
    println!("Distance to z={}: {:.1} Mpc", z, distance);
}
```

See [Quick Start Guide](docs/QUICK_START.md) for more examples.

## Documentation

- [Quick Start Guide](docs/QUICK_START.md) - Get started in 5 minutes
- [User Guide](docs/USER_GUIDE.md) - Comprehensive usage documentation
- [API Reference](docs/API_REFERENCE.md) - Complete API documentation
- [Project Overview](docs/PROJECT_OVERVIEW.md) - Project structure and roadmap
- [Phase 1 Details](docs/PHASE_1_FOUNDATION.md) - Foundation implementation
- [Phase 2 Details](docs/PHASE_2_CORE_COSMOLOGY.md) - Core cosmology features
- [Phase 3 Details](docs/PHASE_3_ADVANCED_FEATURES.md) - Advanced features

Generate Rust docs:
```bash
cargo doc --open
```

## Installation

### From Source

```bash
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
cargo test
```

### Requirements

- Rust 1.70 or later
- Dependencies managed by Cargo (see Cargo.toml)

## Examples

The `examples/` directory contains many demonstration programs:

### Basic Examples

```bash
# Hubble diagram for Type Ia supernovae
cargo run --example hubble_diagram

# Universe evolution from Big Bang to today
cargo run --example universe_evolution

# Hubble parameter evolution
cargo run --example hubble_parameter

# Distance measures comparison
cargo run --example distance_measures

# CMB recombination history
cargo run --example recombination
```

### Advanced Examples

```bash
# CMB angular power spectrum
cargo run --example cmb_power_spectrum

# Matter power spectrum
cargo run --example matter_power_spectrum

# Structure formation growth factor
cargo run --example structure_growth
```

## Module Overview

### Core Modules

- `constants` - Physical and cosmological constants
- `units` - Unit conversion utilities
- `dynamics` - Universe models and evolution (Friedmann equations)
- `observations` - Distance measures and observables

### Physics Modules

- `cmb` - Cosmic Microwave Background
  - `recombination` - Ionization history (Saha equation)
  - `fluctuations` - Angular power spectrum C_ℓ
- `structure` - Large-scale structure
  - `power_spectrum` - Matter power spectrum P(k)
- `perturbations` - Linear perturbation theory
  - `growth` - Growth factor D(a) and growth rate f
  - `boltzmann` - Simplified Boltzmann solver
- `advanced` - Advanced features
  - `weak_lensing` - Convergence κ and shear γ

### Visualization

- `visualization` - Plotting utilities
  - `plots_2d` - Static PNG plots
  - `plotly_plots` - Interactive HTML plots
  - `colors` - Color schemes

## Usage Examples

### Calculate Universe Age

```rust
use andam::prelude::*;

let universe = Universe::benchmark();
let age = universe.age_today();
println!("Universe age: {:.2} Gyr", age); // ~13.80 Gyr
```

### Distance Calculations

```rust
use andam::prelude::*;

let universe = Universe::benchmark();
let z = 1.0;

let d_c = comoving_distance(z, &universe);
let d_l = luminosity_distance(z, &universe);
let d_a = angular_diameter_distance(z, &universe);

println!("Distances to z={}:", z);
println!("  Comoving: {:.1} Mpc", d_c);
println!("  Luminosity: {:.1} Mpc", d_l);
println!("  Angular diameter: {:.1} Mpc", d_a);
```

### CMB Properties

```rust
use andam::prelude::*;

let universe = Universe::benchmark();

let z_rec = recombination_redshift(&universe);
println!("Recombination at z = {:.0}", z_rec); // ~1091

let x_e = ionization_fraction(1100.0, &universe);
println!("Ionization fraction at z=1100: {:.4}", x_e);
```

### Matter Power Spectrum

```rust
use andam::structure::power_spectrum::*;

let omega_m = 0.3111;
let omega_b = 0.049;
let h = 0.6766;
let a_s = 2.1e-9;
let n_s = 0.9665;

let k = 0.1; // h/Mpc
let delta_sq = dimensionless_power(k, 0.0, omega_m, omega_b, h, a_s, n_s);
println!("Δ²(k={}) = {:.4}", k, delta_sq);
```

### Create Plots

```rust
use andam::prelude::*;
use andam::visualization::plots_2d::*;

let universe = Universe::benchmark();

let data: Vec<(f64, f64)> = (1..=100)
    .map(|i| {
        let a = i as f64 / 100.0;
        let h = universe.hubble(a);
        (a, h)
    })
    .collect();

create_line_plot(
    &data,
    "hubble_evolution.png",
    "Hubble Parameter Evolution",
    "Scale Factor a",
    "H(a) [km/s/Mpc]"
).unwrap();
```

## Testing

Run the full test suite:

```bash
cargo test
```

Run specific test groups:

```bash
cargo test --lib                    # Unit tests
cargo test --test phase1_tests      # Phase 1 integration tests
cargo test --test phase2_tests      # Phase 2 integration tests
cargo test --test phase3_tests      # Phase 3 integration tests
```

Run with output:

```bash
cargo test -- --nocapture
```

## Benchmarks

Run performance benchmarks:

```bash
cargo bench
```

## Scientific Validation

The library has been validated against:

- Planck 2018 cosmological parameters
- Standard ΛCDM model predictions
- Analytical solutions where available
- Published cosmological data

Key validations:
- Universe age: 13.80 ± 0.02 Gyr (matches Planck 2018)
- CMB recombination: z ≈ 1091 (standard value)
- Distance-redshift relations verified against CAMB/CLASS
- Power spectrum matches Eisenstein-Hu (1998) predictions

## Performance

Typical computation times on modern CPU (single core):

- Universe age: ~1 μs
- Comoving distance: ~10 μs
- Matter power spectrum: ~5 μs
- Growth factor: ~100 μs
- CMB power spectrum (l_max=2000): ~5 ms

Memory usage: < 1 MB for typical calculations

## Roadmap

### Version 0.1.0 (Current)

- Core cosmology calculations
- Distance measures
- CMB recombination and fluctuations
- Matter power spectrum
- Linear perturbation theory
- Basic weak lensing
- 2D and interactive visualizations

### Version 0.2.0 (Planned)

- Non-linear power spectrum
- Halo mass function
- Baryon acoustic oscillations (detailed)
- Neutrino physics
- Reionization modeling
- 3D visualizations (kiss3d)

### Version 0.3.0 (Future)

- GPU acceleration
- Machine learning integration
- Real-time parameter estimation
- Interactive web interface

### Version 1.0.0 (Long-term)

- Production-ready stability
- Complete API stability
- Comprehensive benchmarking
- Professional documentation

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas where contributions are especially valuable:

- Additional cosmological models
- Performance optimizations
- More visualization options
- Additional examples
- Documentation improvements
- Bug fixes and tests

## Dependencies

Key dependencies:

- `ndarray` - N-dimensional arrays
- `nalgebra` - Linear algebra
- `plotters` - Static plots
- `plotly` - Interactive plots
- `rayon` - Parallel computing
- `approx` - Approximate comparisons (tests)

See [Cargo.toml](Cargo.toml) for complete list.

## Citation

If you use this library in academic work:

```bibtex
@software{andam,
  title = {Andam: Cosmological Calculations in Rust},
  author = {Cosmos Andam Contributors},
  year = {2025},
  url = {https://github.com/cosmos-andam/andam},
  note = {Based on Ryden (2016) and Dodelson (2003)}
}
```

## References

### Textbooks

- Ryden, B. (2016). Introduction to Cosmology (2nd ed.). Cambridge University Press.
- Dodelson, S. (2003). Modern Cosmology. Academic Press.

### Data Sources

- Planck Collaboration (2018). Planck 2018 results. VI. Cosmological parameters.
- Eisenstein, D.J. & Hu, W. (1998). Baryonic Features in the Matter Transfer Function. ApJ 496:605.

### Methods

- Simpson's rule for numerical integration
- Euler method for ODE integration
- Limber approximation for lensing
- Kaiser-Squires inversion for shear

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Acknowledgments

### Libraries

- Rust scientific computing ecosystem
- Plotters and Plotly developers
- ndarray and nalgebra maintainers

### Community

- Rust community for excellent tooling
- Cosmology community for open data
- Open source contributors

### Inspiration

This project implements numerical methods and calculations from:
- Barbara Ryden's "Introduction to Cosmology"
- Scott Dodelson's "Modern Cosmology"
- Planck Collaboration publications

## Contact

- GitHub Issues: https://github.com/cosmos-andam/andam/issues
- Discussions: https://github.com/cosmos-andam/andam/discussions
- Documentation: https://docs.rs/andam

## Status

Version: 0.1.0

Status: Active development

Tests: 40 passing

Coverage: Core modules fully tested

---

Built with Rust for performance, safety, and reliability in cosmological calculations.
