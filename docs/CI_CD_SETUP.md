# CI/CD Pipeline Documentation

This document describes the CI/CD pipeline for the Andam project.

## GitHub Actions Workflows

### 1. CI Workflow (ci.yml)

**Triggers**: Push to main/develop, Pull requests

**Jobs**:
- **Test Suite**: Runs on Linux, macOS, Windows with stable and nightly Rust
  - Tests with default features
  - Tests with no default features
  - Tests with plotting only
  - Runs doc tests

- **Test with HDF5**: Ubuntu only, installs HDF5 system library
  - Tests with hdf5-storage feature

- **Clippy**: Linting with clippy
  - Checks default features
  - Checks no default features
  - Fails on warnings

- **Rustfmt**: Code formatting check
  - Ensures consistent code style

- **Documentation**: Verifies docs build without warnings
  - Simulates docs.rs build

- **Coverage**: Code coverage with tarpaulin
  - Uploads to Codecov

- **MSRV**: Minimum Supported Rust Version check
  - Currently set to Rust 1.70

### 2. Release Workflow (release.yml)

**Triggers**: Git tags matching `v*.*.*` (e.g., v0.1.2)

**Jobs**:
- **Publish**: Publishes to crates.io
  - Verifies tag matches Cargo.toml version
  - Runs tests
  - Checks documentation
  - Publishes with CARGO_REGISTRY_TOKEN secret

- **Create Release**: Creates GitHub release
  - Extracts changelog for the version
  - Creates release notes from CHANGELOG.md
  - Attaches to tag

### 3. Security Workflow (security.yml)

**Triggers**:
- Weekly schedule (Monday 00:00 UTC)
- Push to main
- Pull requests to main

**Jobs**:
- **Audit**: Runs cargo-audit for known vulnerabilities
- **Deny**: Runs cargo-deny for dependency licensing and security

### 4. Documentation Workflow (docs.yml)

**Triggers**: Push to main, Pull requests

**Jobs**:
- **Build and Deploy**: Builds documentation and deploys to GitHub Pages
  - Creates index.html redirect
  - Deploys to gh-pages branch (main branch only)

### 5. Benchmark Workflow (benchmark.yml)

**Triggers**: Push to main, Pull requests

**Jobs**:
- **Benchmark**: Conditionally runs benchmarks
  - Checks if `[[bench]]` exists in Cargo.toml
  - Only runs if benchmarks are configured
  - Stores results for comparison (on main branch only)
  - Auto-pushes to benchmark branch
  - Shows friendly message if no benchmarks exist

## Configuration Files

### dependabot.yml

Automated dependency updates:
- **Cargo dependencies**: Weekly updates
- **GitHub Actions**: Weekly updates
- Groups all Rust dependencies together
- Limits concurrent PRs

### deny.toml

Cargo-deny configuration (v2 format):
- **Advisories**: Checks for known vulnerabilities in dependencies
- **Licenses**: Allows common open-source licenses (MIT, Apache-2.0, BSD, etc.)
- **Bans**: Warns on multiple versions of the same crate
- **Sources**: Allows crates.io only

Note: Configuration uses cargo-deny v2 format with `version = 2` in each section.

## GitHub Templates

### Pull Request Template

Standard template requiring:
- Description
- Type of change
- Checklist (code style, tests, documentation)
- Testing details
- Related issues

### Issue Templates

1. **Bug Report**: Structured bug reporting with environment details
2. **Feature Request**: Feature proposals with use cases

## Required Secrets

To enable full CI/CD functionality, configure these secrets in GitHub repository settings:

### CARGO_REGISTRY_TOKEN

Required for publishing to crates.io.

**How to get it**:
1. Log in to https://crates.io
2. Go to Account Settings
3. Create a new API token
4. Copy the token
5. Add to GitHub: Settings → Secrets and variables → Actions → New repository secret
   - Name: `CARGO_REGISTRY_TOKEN`
   - Value: [your token]

### GITHUB_TOKEN

Automatically provided by GitHub Actions (no setup needed).

### CODECOV_TOKEN (Optional)

For code coverage uploads to Codecov.

**How to get it**:
1. Sign up at https://codecov.io with your GitHub account
2. Add your repository
3. Copy the upload token
4. Add to GitHub: Settings → Secrets and variables → Actions → New repository secret
   - Name: `CODECOV_TOKEN`
   - Value: [your token]

## Publishing a Release

### Step 1: Update Version

Update version in `Cargo.toml`:
```toml
[package]
version = "0.1.3"  # New version
```

### Step 2: Update CHANGELOG.md

Add release notes:
```markdown
## [0.1.3] - 2026-01-XX

### Added
- Feature description

### Fixed
- Bug fix description
```

### Step 3: Commit and Tag

```bash
git add Cargo.toml CHANGELOG.md
git commit -m "Bump version to 0.1.3"
git tag v0.1.3
git push origin main
git push origin v0.1.3
```

### Step 4: Automated Process

The release workflow will automatically:
1. Verify tag matches version
2. Run all tests
3. Build documentation
4. Publish to crates.io
5. Create GitHub release with changelog

## Feature Flags in CI

The CI pipeline tests multiple feature combinations:

| Configuration | Use Case |
|---------------|----------|
| Default | Standard user installation |
| No default features | Minimal installation |
| Plotting only | Documentation (docs.rs) |
| HDF5 storage | Users with HDF5 library |

## Caching Strategy

The CI uses GitHub Actions cache for:
- Cargo registry (~/.cargo/registry)
- Cargo git index (~/.cargo/git)
- Build artifacts (target/)

**Cache key**: Based on Cargo.lock hash for reproducibility

## Platform Coverage

| Platform | Rust Version | Features Tested |
|----------|-------------|-----------------|
| Linux | stable, nightly | All combinations |
| macOS | stable | Default + plotting |
| Windows | stable | Default + plotting |

## Performance

Typical CI run times:
- Test suite: ~5-10 minutes
- Clippy: ~2-3 minutes
- Documentation: ~2-3 minutes
- Full pipeline: ~15-20 minutes

## System Dependencies

The CI requires system libraries to be installed:

### Linux (Ubuntu)
- **libfontconfig1-dev**: Required by plotters for font rendering
- **pkg-config**: Required for library detection (CRITICAL - must be installed)
- **libhdf5-dev**: Required for HDF5 storage tests (test-hdf5 job only)

These are automatically installed in GitHub Actions workflows.

IMPORTANT: Both libfontconfig1-dev AND pkg-config must be installed together. Missing pkg-config will cause build failures even if libfontconfig1-dev is installed.

### Local Development

On Ubuntu/Debian:
```bash
sudo apt-get install libfontconfig1-dev pkg-config
```

On macOS:
```bash
# fontconfig usually comes with the system
# or install via Homebrew if needed
brew install fontconfig
```

## Troubleshooting

### CI Failing on Feature Tests

Check that new features are properly configured:
1. Added to `Cargo.toml` [features]
2. Conditional compilation in code
3. CI updated to test the feature

### Release Workflow Failing

Common issues:
- Tag doesn't match Cargo.toml version
- CARGO_REGISTRY_TOKEN not set
- Tests failing
- Documentation has warnings

### HDF5 Tests Failing

HDF5 tests only run on Ubuntu with libhdf5-dev installed. If failing:
- Check HDF5 version compatibility
- Verify pkg-config setup
- Check feature flag configuration

### Fontconfig Build Failures

If you see "fontconfig required by yeslogic-fontconfig-sys was not found" or "The file `fontconfig.pc` needs to be installed":
- CRITICAL: Install BOTH libfontconfig1-dev AND pkg-config together:
  ```bash
  sudo apt-get install -y libfontconfig1-dev pkg-config
  ```
- Verify the job includes system dependency installation step
- Check PKG_CONFIG_PATH environment variable if still failing
- Common cause: Installing libfontconfig1-dev without pkg-config

### Cargo-Deny Configuration Errors

If you see "unexpected-value" or "deprecated" errors from cargo-deny:
- Ensure deny.toml uses v2 format with `version = 2` in each section
- Remove deprecated keys like `unlicensed`, `copyleft`, `default`, `allow-osi-fsf-free`
- Security workflow uses EmbarkStudios/cargo-deny-action@v1 for better compatibility
- See deny.toml for the correct modern configuration format

### Security Advisory Ignores

The project currently ignores three security advisories due to dependency constraints:

**RUSTSEC-2020-0071** (time 0.1.45 vulnerability):
- Transitive dependency from plotly 0.8
- Updating plotly to 0.14+ requires Rust 1.88 (not yet stable)
- Will be fixed when Rust 1.88 is released

**RUSTSEC-2024-0384** (instant unmaintained):
- Transitive from hdf5 0.8 via parking_lot
- Latest hdf5 version still uses this dependency
- Monitoring for hdf5 updates

**RUSTSEC-2024-0436** (paste unmaintained):
- Transitive from simba/hdf5 used by nalgebra/ode_solvers
- Monitoring for upstream migrations

These are documented in deny.toml with TODOs and will be addressed when updates become available.

## Adding Benchmarks

Currently no benchmarks are configured. To add benchmarks:

1. Uncomment the benchmark section in `Cargo.toml`:
```toml
[[bench]]
name = "friedmann_bench"
harness = false
```

2. Create benchmark file in `benches/`:
```bash
mkdir benches
# Create benchmark files
```

3. The benchmark workflow will automatically run them on the next push

## Future Enhancements

Potential additions to CI/CD:
- [ ] Performance regression testing (when benchmarks added)
- [ ] Automated benchmark comparisons
- [ ] Multi-architecture testing (ARM, etc.)
- [ ] Container image publishing
- [ ] Automated changelogs
- [ ] Release candidate builds

## Badge Status

Add these badges to README.md (already configured):

```markdown
[![CI](https://github.com/cosmos-andam/andam/workflows/CI/badge.svg)](https://github.com/cosmos-andam/andam/actions/workflows/ci.yml)
[![Security Audit](https://github.com/cosmos-andam/andam/workflows/Security%20Audit/badge.svg)](https://github.com/cosmos-andam/andam/actions/workflows/security.yml)
[![Crates.io](https://img.shields.io/crates/v/andam.svg)](https://crates.io/crates/andam)
[![Documentation](https://docs.rs/andam/badge.svg)](https://docs.rs/andam)
```
