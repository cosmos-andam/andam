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
- **Benchmark**: Runs cargo bench if benchmarks configured
  - Stores results for comparison
  - Auto-pushes to benchmark branch

## Configuration Files

### dependabot.yml

Automated dependency updates:
- **Cargo dependencies**: Weekly updates
- **GitHub Actions**: Weekly updates
- Groups all Rust dependencies together
- Limits concurrent PRs

### deny.toml

Cargo-deny configuration:
- **Advisories**: Denies known vulnerabilities
- **Licenses**: Allows common open-source licenses
- **Bans**: Warns on multiple versions
- **Sources**: Allows crates.io only

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

## Future Enhancements

Potential additions to CI/CD:
- [ ] Performance regression testing
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
