# Project Rename: cosmology-rs → Andam

## What Changed

All implementation documents have been updated to reflect the new project name **"Andam"**.

### Project Name
- **Old**: cosmology-rs / cosmology_rs
- **New**: andam (అండం - "universe" or "cosmos" in Telugu)

## Updated Elements

### 1. Package Names
```toml
# Before
name = "cosmology-rs"

# After
name = "andam"
```

### 2. Rust Imports
```rust
// Before
use cosmology_rs::dynamics::Universe;
use cosmology_rs::prelude::*;

// After
use andam::dynamics::Universe;
use andam::prelude::*;
```

### 3. Directory Structure
```
# Before
cosmology-rs/

# After
andam/
```

### 4. URLs and Links
- **GitHub**: `github.com/yourusername/andam`
- **Docs**: `docs.rs/andam`
- **Crates.io**: `crates.io/crates/andam`

### 5. Documentation Headers
- Updated all markdown file titles
- Updated project descriptions
- Added Telugu etymology note

## Files Updated

All 5 implementation documents have been updated:
1. [DONE] `PROJECT_OVERVIEW.md`
2. [DONE] `PHASE_1_FOUNDATION.md`
3. [DONE] `PHASE_2_CORE_COSMOLOGY.md`
4. [DONE] `PHASE_3_ADVANCED_FEATURES.md`
5. [DONE] `PHASE_4_POLISH_DOCUMENTATION.md`

## Quick Start Commands (Updated)

```bash
# Clone repository
git clone https://github.com/yourusername/andam
cd andam

# Build
cargo build --release

# Run tests
cargo test

# Generate documentation
cargo doc --open

# Run example
cargo run --example hubble_diagram
```

## Installation (Updated)

```toml
[dependencies]
andam = "0.1.0"
```

```rust
use andam::prelude::*;

fn main() {
 let universe = Universe::benchmark();
 println!("Universe age: {:.2} Gyr", universe.age_today());
}
```

## Next Steps

1. Download all updated markdown files
2. Start implementation with Phase 1
3. Follow the same structured approach
4. Create your project with `cargo new andam --lib`

---

**Note**: All code examples, imports, and references throughout the documentation have been consistently updated to use "andam" instead of "cosmology-rs".
