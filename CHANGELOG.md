# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-01-17

### Fixed
- docs.rs build failure caused by plotly_kaleido requiring network access
- Documentation warnings from unescaped square brackets in unit descriptions

### Changed
- Made kaleido an optional feature flag (`plotting-kaleido`)
- Separated plotting into base (`plotting`) and kaleido-enhanced (`plotting-kaleido`)
- Updated `full` feature to include `plotting-kaleido` instead of just `plotting`
- Improved README with clearer feature flag documentation

### Added
- New `plotting-kaleido` feature flag for static image export
- Documentation explaining kaleido network requirements
- DOCSRS_FIX.md explaining the build issue and solution

## [0.1.0] - 2026-01-17 (yanked due to docs.rs build failure)

### Added
- Initial release with complete Phases 1-9 implementation
- Core cosmology calculations (ΛCDM, Friedmann equations)
- CMB physics (recombination, power spectra, polarization)
- Structure formation (P(k), halos, correlation functions, cosmic web)
- Big Bang Nucleosynthesis (BBN) with primordial abundances
- Statistical analysis (MCMC, Fisher matrix, parameter estimation)
- Beyond-ΛCDM models (dark energy, massive neutrinos)
- Visualization capabilities (plotly, plotters)
- HDF5 data storage (optional feature)
- 60+ unit tests
- 15 comprehensive examples
- Full documentation

### Notes
- v0.1.0 was published but docs failed to build on docs.rs
- Issue was fixed in v0.1.1 within hours of discovery
