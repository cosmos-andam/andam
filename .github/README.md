# GitHub Configuration

This directory contains GitHub-specific configuration for the Andam project.

## Structure

```
.github/
├── workflows/              # GitHub Actions workflows
│   ├── ci.yml             # Main CI pipeline
│   ├── release.yml        # Automated releases
│   ├── security.yml       # Security audits
│   ├── docs.yml           # Documentation deployment
│   └── benchmark.yml      # Performance benchmarks
├── ISSUE_TEMPLATE/        # Issue templates
│   ├── bug_report.md      # Bug report template
│   └── feature_request.md # Feature request template
├── PULL_REQUEST_TEMPLATE.md # PR template
├── dependabot.yml         # Dependency updates
└── README.md             # This file
```

## Workflows

### CI (ci.yml)
- Runs on every push and PR
- Tests on Linux, macOS, Windows
- Multiple Rust versions (stable, nightly)
- Feature flag combinations
- Clippy, rustfmt, docs checks
- Code coverage

### Release (release.yml)
- Triggered by version tags (v*)
- Publishes to crates.io
- Creates GitHub releases
- Extracts changelog

### Security (security.yml)
- Weekly security audits
- cargo-audit and cargo-deny
- Runs on main branch changes

### Documentation (docs.yml)
- Builds and deploys docs
- GitHub Pages deployment
- Runs on main branch

### Benchmark (benchmark.yml)
- Performance regression tracking
- Stores benchmark history

## Issue Templates

### Bug Report
Structured template for reporting bugs with:
- Description
- Reproduction steps
- Environment details
- Expected vs actual behavior

### Feature Request
Template for proposing new features with:
- Feature description
- Motivation and use case
- Proposed solution
- Implementation notes

## Pull Request Template

Standardized PR template requiring:
- Change description
- Type of change checklist
- Testing details
- Related issues

## Dependabot

Automated dependency updates:
- Cargo dependencies: weekly
- GitHub Actions: weekly
- Grouped updates for related packages

## Required Secrets

Configure these in repository settings (Settings → Secrets and variables → Actions):

### CARGO_REGISTRY_TOKEN
Required for publishing to crates.io
- Get from: https://crates.io/settings/tokens
- Required for: release.yml

### CODECOV_TOKEN (optional)
For code coverage uploads
- Get from: https://codecov.io
- Required for: ci.yml coverage job

## Badges

Available badges for README:

```markdown
[![CI](https://github.com/cosmos-andam/andam/workflows/CI/badge.svg)](https://github.com/cosmos-andam/andam/actions/workflows/ci.yml)
[![Security Audit](https://github.com/cosmos-andam/andam/workflows/Security%20Audit/badge.svg)](https://github.com/cosmos-andam/andam/actions/workflows/security.yml)
[![Crates.io](https://img.shields.io/crates/v/andam.svg)](https://crates.io/crates/andam)
[![Documentation](https://docs.rs/andam/badge.svg)](https://docs.rs/andam)
[![License](https://img.shields.io/crates/l/andam.svg)](https://github.com/cosmos-andam/andam#license)
[![Downloads](https://img.shields.io/crates/d/andam.svg)](https://crates.io/crates/andam)
```

## Local Testing

Test CI configuration locally:

```bash
# Format check
cargo fmt --all -- --check

# Clippy
cargo clippy --all-targets --features plotting -- -D warnings

# Tests
cargo test --features plotting

# Documentation
cargo doc --no-deps --features plotting
```

## Troubleshooting

### Workflow not running
- Check workflow triggers in YAML
- Verify branch names match
- Check repository settings for Actions

### Secret not working
- Verify secret name matches YAML
- Check secret scope (repository vs environment)
- Secrets are not shown in logs

### Tests failing in CI but passing locally
- Check Rust version differences
- Verify feature flags
- Check platform-specific code
- Review CI logs for details

## References

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Rust CI Best Practices](https://doc.rust-lang.org/cargo/guide/continuous-integration.html)
- [Dependabot Configuration](https://docs.github.com/en/code-security/dependabot)
