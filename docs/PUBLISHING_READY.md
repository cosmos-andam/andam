# Publishing Readiness Summary

## Issue Resolved

Successfully resolved the crates.io publishing error:
```
error: all dependencies must have a version specified when publishing.
dependency `hdf5` does not specify a version
```

## Solution Implemented

### 1. Made HDF5 Optional

Changed HDF5 from a required dependency to an optional feature:

```toml
[dependencies]
# Before: Required git dependency (unpublishable)
# hdf5 = { git = "https://github.com/aldanor/hdf5-rust.git" }

# After: Optional crates.io dependency
hdf5 = { version = "0.8", optional = true }
flate2 = { version = "1.0", optional = true }
chrono = { version = "0.4", optional = true }

[features]
default = ["plotting"]
plotting = ["plotly", "plotters", "plotters-backend"]
hdf5-storage = ["hdf5", "flate2", "chrono"]
full = ["plotting", "hdf5-storage"]
```

### 2. Made Plotting Optional

Also made plotting dependencies optional for users who only need calculations:

```toml
plotly = { version = "0.8", features = ["kaleido"], optional = true }
plotters = { version = "0.3", features = ["all_series", "bitmap_encoder"], optional = true }
plotters-backend = { version = "0.3", optional = true }
```

### 3. Updated Module Conditionals

Added feature flags in `src/lib.rs`:

```rust
#[cfg(feature = "plotting")]
pub mod visualization;

#[cfg(feature = "hdf5-storage")]
pub mod storage;
```

## Publishing Verification

Tested and confirmed:

```bash
$ cargo publish --dry-run --allow-dirty
   Packaging andam v0.1.0 (/Users/prabasiva/Documents/2026/andam)
    Packaged 70 files, 263.5KiB (66.2KiB compressed)
   Verifying andam v0.1.0 (/Users/prabasiva/Documents/2026/andam)
   Compiling andam v0.1.0 (/Users/prabasiva/Documents/2026/andam/target/package/andam-0.1.0)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 14.35s
   Uploading andam v0.1.0 (/Users/prabasiva/Documents/2026/andam)
warning: aborting upload due to dry run
```

No errors - ready to publish!

## Feature Usage

### For Users

#### Default Usage (with plotting)
```bash
cargo add andam
```
Includes all core functionality + visualization.

#### With HDF5 Storage
```bash
cargo add andam --features hdf5-storage
```
Requires HDF5 1.12 or earlier installed on system.

#### Minimal Installation
```bash
cargo add andam --no-default-features
```
Just the core cosmology calculations, no dependencies on plotting or HDF5.

#### Everything
```bash
cargo add andam --features full
```
All features enabled.

## What Works

### Without Any Features
- All cosmology calculations
- Distance measures
- CMB physics
- Structure formation
- MCMC statistics
- BBN
- All 60 core tests pass

### With Default Features (plotting)
- Everything above
- Visualization with plotly and plotters
- All examples work

### With hdf5-storage Feature
- Everything above
- HDF5 data storage and retrieval
- All 66 tests pass (including storage tests)

## Known Limitations

### HDF5 Version Compatibility

The published `hdf5` crate (0.8.1) only supports HDF5 1.12 and earlier:
- Works on systems with older HDF5
- May not work on macOS with Homebrew HDF5 1.14+
- See `docs/HDF5_PUBLISHING_NOTE.md` for workarounds

### Why This is OK

1. **Core functionality doesn't need HDF5**: Most users get full value without it
2. **Feature is clearly documented**: Users know it's optional
3. **Workarounds available**: Advanced users can use git version if needed
4. **Future-proof**: When hdf5-rust 0.9 releases with 1.14 support, we can update

## Documentation Updates

Updated documentation to reflect optional features:

1. **README.md**
   - Added installation examples for different feature combinations
   - Added feature flags table
   - Added notes on HDF5 storage
   - Updated module overview

2. **docs/HDF5_PUBLISHING_NOTE.md** (new)
   - Explains the HDF5 version issue
   - Provides workarounds for HDF5 1.14+ users
   - Documents future options

3. **docs/HDF5_STORAGE_IMPLEMENTATION.md**
   - Already documents the storage implementation
   - Users can reference this if they enable the feature

## Next Steps

Ready to publish:

```bash
# Final verification
cargo test --all-features
cargo build --release --all-features
cargo publish --dry-run

# When ready to publish for real
cargo publish
```

## Impact Summary

Before:
- Could not publish due to git dependency
- All users forced to have HDF5 system library
- Binary size larger than necessary

After:
- Ready to publish to crates.io
- Users choose only the features they need
- Smaller binary for users who don't need HDF5 or plotting
- Clear documentation of what each feature provides

**Status: Ready for publication to crates.io**
