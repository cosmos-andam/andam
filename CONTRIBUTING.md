# Contributing to Andam

Thank you for your interest in contributing to Andam! This document provides guidelines and instructions for contributing.

## Code of Conduct

Be respectful, inclusive, and professional in all interactions.

## Getting Started

### Prerequisites

- Rust 1.70 or later
- Git
- Familiarity with cosmology concepts is helpful but not required

### Development Setup

1. Fork and clone the repository:
```bash
git clone https://github.com/YOUR_USERNAME/andam.git
cd andam
```

2. Build the project:
```bash
cargo build
```

3. Run tests:
```bash
cargo test
```

4. Run tests with all features:
```bash
# Default features (plotting)
cargo test --features plotting

# With HDF5 (requires system library)
cargo test --features hdf5-storage
```

## Development Workflow

### 1. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

Use descriptive branch names:
- `feature/` for new features
- `fix/` for bug fixes
- `docs/` for documentation
- `refactor/` for code refactoring
- `test/` for test improvements

### 2. Make Changes

- Write clear, documented code
- Follow Rust conventions
- Add tests for new functionality
- Update documentation as needed

### 3. Run Quality Checks

Before committing, ensure your code passes all checks:

```bash
# Format code
cargo fmt

# Run clippy
cargo clippy -- -D warnings

# Run tests
cargo test

# Build documentation
cargo doc --no-deps --features plotting
```

### 4. Commit Changes

Write clear commit messages:
```bash
git commit -m "Add feature: brief description"
```

Good commit message format:
```
Add feature: Calculate neutrino density evolution

- Implement neutrino decoupling calculations
- Add tests for temperature evolution
- Update documentation with examples
```

### 5. Push and Create PR

```bash
git push origin your-branch-name
```

Then create a pull request on GitHub.

## Code Style Guidelines

### Rust Style

- Follow [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Use `rustfmt` for formatting (run `cargo fmt`)
- Use `clippy` for linting (run `cargo clippy`)
- Maximum line length: 100 characters (enforced by rustfmt)

### Documentation

- All public items must have doc comments
- Include examples in doc comments where applicable
- Use `///` for item documentation
- Use `//!` for module documentation
- **NO EMOJIS** in documentation (per project standards)

Example:
```rust
/// Calculate the Hubble parameter at redshift z
///
/// The Hubble parameter describes the expansion rate of the universe
/// at a given redshift.
///
/// # Arguments
///
/// * `z` - Redshift value (must be >= 0)
/// * `universe` - Cosmological parameters
///
/// # Returns
///
/// Hubble parameter in km/s/Mpc
///
/// # Examples
///
/// ```
/// use andam::dynamics::Universe;
/// use andam::observations::hubble_parameter;
///
/// let universe = Universe::benchmark();
/// let h_z = hubble_parameter(1.0, &universe);
/// assert!(h_z > 0.0);
/// ```
pub fn hubble_parameter(z: f64, universe: &Universe) -> f64 {
    // Implementation
}
```

### Testing

- Write unit tests for all functions
- Use `#[test]` for test functions
- Use `approx::assert_relative_eq!` for floating-point comparisons
- Group related tests in test modules

Example:
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hubble_parameter_at_z_zero() {
        let universe = Universe::benchmark();
        let h0 = hubble_parameter(0.0, &universe);
        assert_relative_eq!(h0, universe.h0, epsilon = 1e-10);
    }
}
```

## Project Structure

```
andam/
├── src/
│   ├── lib.rs              # Library root
│   ├── constants.rs        # Physical constants
│   ├── dynamics/           # Friedmann equations
│   ├── observations/       # Observable calculations
│   ├── cmb/               # CMB physics
│   ├── structure/         # Structure formation
│   ├── early_universe/    # BBN and early physics
│   ├── statistics/        # MCMC and parameter estimation
│   ├── beyond_lcdm/       # Dark energy, neutrinos
│   ├── visualization/     # Plotting utilities
│   └── storage/           # HDF5 storage (optional)
├── examples/              # Example programs
├── docs/                  # Documentation
└── tests/                 # Integration tests
```

## Adding New Features

### 1. Plan Your Feature

- Open an issue to discuss the feature
- Get feedback before implementing
- Consider how it fits with existing code

### 2. Implement the Feature

- Add code in appropriate module
- Follow existing patterns
- Keep functions focused and small
- Use descriptive names

### 3. Add Tests

- Unit tests in the same file
- Integration tests in `tests/`
- Aim for >80% code coverage

### 4. Document

- Add doc comments
- Update relevant documentation in `docs/`
- Add examples if applicable
- Update CHANGELOG.md

### 5. Feature Flags

If your feature requires optional dependencies:

1. Add dependency to `Cargo.toml`:
```toml
[dependencies]
your-crate = { version = "1.0", optional = true }

[features]
your-feature = ["your-crate"]
```

2. Use conditional compilation:
```rust
#[cfg(feature = "your-feature")]
pub mod your_module;
```

3. Update CI to test the feature:
```yaml
- name: Test your feature
  run: cargo test --features your-feature
```

## Fixing Bugs

### 1. Reproduce the Bug

- Create a minimal test case
- Document the expected vs actual behavior

### 2. Fix the Bug

- Make the smallest change necessary
- Ensure the fix doesn't break existing tests

### 3. Add Regression Test

- Add a test that would fail without your fix
- Prevents the bug from reoccurring

## Documentation

### Types of Documentation

1. **Code documentation** (`///` comments)
2. **Module documentation** (`//!` at the top of files)
3. **User guides** (in `docs/` directory)
4. **Examples** (in `examples/` directory)

### Writing Good Examples

- Keep examples self-contained
- Show realistic use cases
- Include error handling
- Add comments explaining non-obvious steps

## CI/CD Pipeline

All pull requests automatically run:
- Tests on Linux, macOS, Windows
- Clippy linting
- Code formatting check
- Documentation build
- Security audit

Ensure your changes pass all checks before requesting review.

## Review Process

1. Automated CI checks must pass
2. At least one maintainer approval required
3. Address review feedback
4. Squash commits if requested

## Release Process

Maintainers handle releases. The process:
1. Update version in `Cargo.toml`
2. Update `CHANGELOG.md`
3. Create git tag
4. Push tag (triggers automated release)

## Questions?

- Open an issue for questions
- Join discussions in existing issues
- Check documentation in `docs/`

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (MIT OR Apache-2.0).
