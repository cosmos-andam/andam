# Quick Start: CI/CD Setup

Get the Andam CI/CD pipeline running in your fork or clone.

## Prerequisites

- GitHub account
- Repository access
- crates.io account (for publishing)

## Step 1: Fork or Clone

```bash
# Fork on GitHub first, then:
git clone https://github.com/YOUR_USERNAME/andam.git
cd andam
```

## Step 2: Verify Local Setup

Ensure everything works locally:

```bash
# Build
cargo build

# Run tests
cargo test

# Check formatting
cargo fmt --check

# Run clippy
cargo clippy -- -D warnings

# Build docs
cargo doc --no-deps --features plotting
```

All should pass without errors.

## Step 3: Enable GitHub Actions

GitHub Actions should be enabled by default. Verify:

1. Go to repository Settings
2. Click Actions → General
3. Ensure "Allow all actions" is selected
4. Click "Save"

## Step 4: Configure Secrets (Optional)

Only needed if you plan to publish or use code coverage.

### For Publishing to crates.io

1. Get your crates.io API token:
   - Log in to https://crates.io
   - Go to Account Settings → API Tokens
   - Click "New Token"
   - Name it "GitHub Actions - andam"
   - Copy the token

2. Add to GitHub:
   - Go to repository Settings → Secrets and variables → Actions
   - Click "New repository secret"
   - Name: `CARGO_REGISTRY_TOKEN`
   - Value: [paste your token]
   - Click "Add secret"

### For Code Coverage (Optional)

1. Sign up at https://codecov.io with GitHub
2. Add your repository
3. Copy the upload token
4. Add to GitHub as `CODECOV_TOKEN` secret (same process as above)

## Step 5: First Push

Make a small change and push:

```bash
# Create branch
git checkout -b test-ci

# Make a small change (e.g., update README)
echo "\n## Testing CI" >> README.md

# Commit and push
git add README.md
git commit -m "Test: Verify CI pipeline"
git push origin test-ci
```

## Step 6: Verify CI is Running

1. Go to your repository on GitHub
2. Click "Actions" tab
3. You should see workflows running
4. Click on a workflow to see details

Expected workflows:
- CI (Test Suite, Clippy, Rustfmt, Docs)
- Security Audit

## Step 7: Create Pull Request

1. GitHub will suggest creating a PR
2. Click "Compare & pull request"
3. Fill in the PR template
4. Click "Create pull request"

The CI will run all checks on your PR.

## Step 8: Monitor Results

Watch the CI checks:
- Green checkmarks = passing
- Red X = failing
- Click "Details" to see logs

All checks must pass before merging.

## Common Issues

### Tests Failing

```bash
# Run locally to reproduce
cargo test --features plotting --verbose
```

### Clippy Warnings

```bash
# Fix clippy issues
cargo clippy --fix --allow-dirty
cargo clippy -- -D warnings
```

### Format Issues

```bash
# Auto-format
cargo fmt
```

### Documentation Warnings

```bash
# Check docs
cargo doc --no-deps --features plotting
```

## Publishing a Release

### Prerequisites
- `CARGO_REGISTRY_TOKEN` secret configured
- All CI checks passing
- Version updated in Cargo.toml
- CHANGELOG.md updated

### Process

1. Update version:
```toml
# Cargo.toml
version = "0.1.3"
```

2. Update changelog:
```markdown
# CHANGELOG.md
## [0.1.3] - 2026-01-XX

### Added
- New feature description
```

3. Commit and tag:
```bash
git add Cargo.toml CHANGELOG.md
git commit -m "Bump version to 0.1.3"
git tag v0.1.3
git push origin main
git push origin v0.1.3
```

4. Watch automated release:
   - Go to Actions tab
   - Find "Release" workflow
   - Monitor progress

The workflow will:
- Run all tests
- Publish to crates.io
- Create GitHub release

## Workflow Overview

### On Every Push/PR
✓ Build and test on 3 platforms
✓ Run clippy linting
✓ Check code formatting
✓ Verify documentation builds
✓ Security audit
✓ Code coverage (optional)

### On Release Tag
✓ Publish to crates.io
✓ Create GitHub release
✓ Extract changelog

### Weekly
✓ Dependency security audit
✓ License compliance check

## Next Steps

1. Review workflow configurations in `.github/workflows/`
2. Read full CI/CD documentation in `docs/CI_CD_SETUP.md`
3. Check CONTRIBUTING.md for development guidelines
4. Enable branch protection rules (optional but recommended)

## Branch Protection (Recommended)

Protect your main branch:

1. Go to Settings → Branches
2. Add rule for `main` branch
3. Enable:
   - Require pull request reviews
   - Require status checks to pass
   - Require branches to be up to date
4. Select required checks:
   - Test Suite
   - Clippy
   - Rustfmt
   - Documentation

This ensures all code is reviewed and tested before merging.

## Getting Help

- Check workflow logs in Actions tab
- Review `.github/workflows/` files
- Read `docs/CI_CD_SETUP.md`
- Open an issue if stuck

## Summary

You now have:
- ✓ Automated testing on 3 platforms
- ✓ Code quality checks (clippy, rustfmt)
- ✓ Documentation verification
- ✓ Security audits
- ✓ Automated publishing (when configured)
- ✓ Dependency updates via Dependabot

Your CI/CD pipeline is ready!
