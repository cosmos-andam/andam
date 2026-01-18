# Justfile for running CI/CD tasks locally
# Install just: https://github.com/casey/just

# Default recipe (show available commands)
default:
    @just --list

# Run all CI checks locally (equivalent to full CI pipeline)
ci-local: format-check clippy test test-hdf5 doc

# Check code formatting
format-check:
    @echo "==> Checking code formatting..."
    cargo fmt --all -- --check

# Run clippy lints
clippy:
    @echo "==> Running clippy (plotting feature)..."
    cargo clippy --all-targets --no-default-features --features plotting -- -D warnings
    @echo "==> Running clippy (no features)..."
    cargo clippy --lib --no-default-features -- -D warnings

# Run all tests
test:
    @echo "==> Running tests (default features)..."
    cargo test --verbose
    @echo "==> Running tests (no default features)..."
    cargo test --no-default-features --verbose
    @echo "==> Running tests (plotting feature)..."
    cargo test --no-default-features --features plotting --verbose
    @echo "==> Running doc tests..."
    cargo test --doc --no-default-features --features plotting

# Run HDF5 tests (requires libhdf5-dev)
test-hdf5:
    @echo "==> Running HDF5 tests..."
    cargo test --no-default-features --features hdf5-storage --verbose

# Check documentation
doc:
    @echo "==> Checking documentation..."
    RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --no-default-features --features plotting --verbose

# Run security audit
security:
    @echo "==> Running cargo audit..."
    cargo audit --ignore RUSTSEC-2020-0071 --ignore RUSTSEC-2024-0384 --ignore RUSTSEC-2024-0436
    @echo "==> Running cargo deny..."
    -cargo deny check

# Check MSRV (Minimum Supported Rust Version)
msrv:
    @echo "==> Checking MSRV (Rust 1.70)..."
    cargo +1.70 check --no-default-features --features plotting --verbose

# Run code coverage
coverage:
    @echo "==> Generating code coverage..."
    cargo tarpaulin --no-default-features --features plotting --out xml --verbose

# Quick check (format, clippy, test - no HDF5)
quick: format-check clippy
    @echo "==> Running quick tests..."
    cargo test --no-default-features --features plotting

# Clean build artifacts
clean:
    cargo clean

# Format code
format:
    cargo fmt --all

# Build with all feature combinations
build-all:
    @echo "==> Building (no features)..."
    cargo build --no-default-features
    @echo "==> Building (plotting)..."
    cargo build --no-default-features --features plotting
    @echo "==> Building (hdf5-storage)..."
    cargo build --no-default-features --features hdf5-storage
    @echo "==> Building (full)..."
    cargo build --no-default-features --features full

# Run benchmarks
bench:
    cargo bench --features plotting

# Install development dependencies
install-deps:
    @echo "Installing development dependencies..."
    cargo install cargo-tarpaulin
    cargo install cargo-audit
    cargo install cargo-deny
