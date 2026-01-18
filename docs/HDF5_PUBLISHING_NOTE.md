# HDF5 Storage and Publishing Notes

## Publishing Issue

When publishing to crates.io, the HDF5 dependency presents a challenge:

### Problem
- The published `hdf5` crate version 0.8.1 only supports HDF5 1.12.x and earlier
- Many systems (macOS via Homebrew, recent Linux distributions) have HDF5 1.14.x installed
- The git version of hdf5-rust supports HDF5 1.14.x but cannot be used for crates.io publishing

### Solution for Publishing

The HDF5 storage feature is now **optional**:

```toml
[features]
default = ["plotting"]
hdf5-storage = ["hdf5", "flate2", "chrono"]  # Optional feature
```

## For Users

### Option 1: Use Without HDF5 (Recommended for most users)

The core Andam library works perfectly without HDF5:

```bash
cargo add andam
```

All cosmology calculations, plotting, and analysis work without HDF5.

### Option 2: Enable HDF5 Storage (Advanced)

If you need HDF5 storage and have HDF5 1.12 or earlier:

```bash
cargo add andam --features hdf5-storage
```

### Option 3: Use Git Version for HDF5 1.14+ (Development)

If you have HDF5 1.14.x and need storage:

```toml
[dependencies]
andam = { version = "0.1", default-features = false, features = ["plotting"] }

# Add HDF5 support manually with git version
hdf5 = { git = "https://github.com/aldanor/hdf5-rust.git" }
flate2 = "1.0"
chrono = "0.4"
```

Then enable the storage module in your code:
```rust
#[cfg(feature = "hdf5")]
use andam::storage::*;
```

## Recommendation

For the initial crates.io release:
- Publish **without** HDF5 in default features
- Core functionality (cosmology, plotting, statistics) works perfectly
- Users who need HDF5 can add it separately following the documentation

## Future Options

1. **Wait for hdf5-rust 0.9**: The next version may support HDF5 1.14
2. **Alternative storage**: Consider adding support for other formats:
   - Apache Arrow/Parquet for columnar data
   - Simple JSON/TOML for small datasets
   - Custom binary format
3. **External crate**: Create a separate `andam-hdf5` crate for storage

## Impact on Publishing

With the current setup:
- `cargo publish` will work ✓
- Default users get full functionality except HDF5 ✓
- Advanced users can enable HDF5 if their system supports it ✓
- Documentation explains the situation ✓

## Checking Before Publish

```bash
# Verify no git dependencies
cargo publish --dry-run

# Should see no errors about git dependencies
```
