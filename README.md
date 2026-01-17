# Andam

A comprehensive Rust library for cosmological calculations and visualizations.

**Andam** (அண்டம்) means "universe" in தமிழ் (Tamil).

## Features

- ΛCDM cosmology with Planck 2018 parameters
- Distance measures (luminosity, angular diameter, comoving)
- CMB recombination and angular power spectrum
- Matter power spectrum and structure formation
- Linear perturbation theory and growth factors
- Weak gravitational lensing
- Static (PNG) and interactive (HTML) visualizations

## Quick Example

```rust
use andam::prelude::*;

fn main() {
    let universe = Universe::benchmark();
    println!("Universe age: {:.2} Gyr", universe.age_today());

    let distance = luminosity_distance(1.0, &universe);
    println!("Distance to z=1: {:.0} Mpc", distance);
}
```

## Installation

```bash
# From source
git clone https://github.com/cosmos-andam/andam
cd andam
cargo build --release
cargo test
```

Add to `Cargo.toml`:
```toml
[dependencies]
andam = "0.1.0"
```

## Documentation

- [Quick Start](docs/QUICK_START.md) - Get started in 5 minutes
- [User Guide](docs/USER_GUIDE.md) - Comprehensive documentation
- [API Reference](docs/API_REFERENCE.md) - Technical API details
- [Project Overview](docs/PROJECT_OVERVIEW.md) - Architecture and roadmap

Generate API docs: `cargo doc --open`

## Examples

Run built-in examples:

```bash
cargo run --example hubble_diagram
cargo run --example universe_evolution
cargo run --example cmb_power_spectrum
cargo run --example matter_power_spectrum
```

See `examples/` directory for 12+ demonstration programs.

## Testing

```bash
cargo test                  # Run all tests
cargo test --lib           # Unit tests only
cargo bench                # Performance benchmarks
```

## Project Structure

```
andam/
├── src/
│   ├── constants.rs       # Physical constants
│   ├── units.rs           # Unit conversions
│   ├── dynamics/          # Universe evolution
│   ├── observations/      # Distance measures
│   ├── cmb/              # CMB physics
│   ├── structure/        # Power spectrum
│   ├── perturbations/    # Growth theory
│   ├── advanced/         # Weak lensing
│   └── visualization/    # Plotting tools
├── examples/             # Demo programs
├── tests/               # Integration tests
└── docs/                # Documentation
```

## Scientific Validation

Validated against:
- Planck 2018 cosmological parameters
- Standard ΛCDM predictions
- Published transfer functions (Eisenstein-Hu 1998)

Key results:
- Universe age: 13.80 Gyr (matches Planck)
- CMB recombination: z ≈ 1091
- Distance measures verified against CAMB/CLASS

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md).

Areas of interest:
- Performance optimizations
- Additional cosmological models
- More examples and documentation
- Bug fixes and tests

## License

Licensed under either:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Citation

```bibtex
@software{andam,
  title = {Andam: Cosmological Calculations in Rust},
  author = {Cosmos Andam Contributors},
  year = {2025},
  url = {https://github.com/cosmos-andam/andam}
}
```

## References

- Ryden, B. (2016). Introduction to Cosmology (2nd ed.)
- Dodelson, S. (2003). Modern Cosmology
- Planck Collaboration (2018). Cosmological parameters

## Contact

- Issues: https://github.com/cosmos-andam/andam/issues
- Discussions: https://github.com/cosmos-andam/andam/discussions
- Documentation: https://docs.rs/andam

---

**Status**: Version 0.1.0 | 40+ tests passing | Active development
