# docs.rs Build Fix: Kaleido Network Requirement

## Problem

The initial v0.1.0 release failed to build on docs.rs with this error:

```
error: failed to run custom build command for `plotly_kaleido v0.8.4`

Caused by:
  warning: spurious network error (3 tries remaining): [6] Could not resolve hostname
  error: download of config.json failed
  Could not resolve host: index.crates.io
```

## Root Cause

The `plotly_kaleido` crate requires network access during its build script to download platform-specific binaries. However, docs.rs runs all builds in an isolated Docker container with `--network none` for security and reproducibility.

## Solution

Changed the plotly dependency structure to make kaleido optional:

### Before (v0.1.0)
```toml
plotly = { version = "0.8", features = ["kaleido"], optional = true }

[features]
plotting = ["plotly", "plotters", "plotters-backend"]
```

### After (v0.1.1+)
```toml
plotly = { version = "0.8", default-features = false, optional = true }

[features]
plotting = ["plotly", "plotters", "plotters-backend"]
plotting-kaleido = ["plotting", "plotly/kaleido"]
full = ["plotting-kaleido", "hdf5-storage"]
```

## Feature Flags Explained

### `plotting` (default)
- Interactive HTML plots with plotly
- Static plots with plotters (PNG, SVG)
- No network required during build
- Works on docs.rs

### `plotting-kaleido`
- Adds static image export for plotly (PNG, SVG, PDF)
- Requires network access during build
- Not compatible with docs.rs
- Only needed if you want to export plotly figures as static images

### docs.rs Configuration
```toml
[package.metadata.docs.rs]
features = ["plotting"]  # Without kaleido
rustdoc-args = ["--cfg", "docsrs"]
```

## Impact on Users

### No Impact
Most users are unaffected:
- Interactive HTML plots work by default
- Plotters provides static PNG/SVG export
- All visualization examples work

### Users Needing Kaleido
If you specifically need plotly static image export:
```toml
[dependencies]
andam = { version = "0.1", features = ["plotting-kaleido"] }
```

## Verification

Test that docs.rs build succeeds:
```bash
# Simulate docs.rs build environment
cargo clean
cargo doc --no-deps --features plotting --offline
```

## Related Issues

- plotly_kaleido issue: https://github.com/plotly/plotly.rs/issues/156
- docs.rs network isolation: https://docs.rs/about/builds

## Summary

This is a build-time dependency issue that only affects offline environments like docs.rs. The fix maintains full functionality while ensuring documentation can be built in isolated environments.
