# Andam API Reference

Complete API documentation for the Andam cosmology library.

## Table of Contents

1. [Module Overview](#module-overview)
2. [Constants](#constants)
3. [Units](#units)
4. [Dynamics](#dynamics)
5. [Observations](#observations)
6. [CMB](#cmb)
7. [Structure](#structure)
8. [Perturbations](#perturbations)
9. [Advanced](#advanced)
10. [Visualization](#visualization)

---

## Module Overview

### Module Hierarchy

```
andam
 constants Physical and cosmological constants
 units Unit conversion utilities
 dynamics Universe models and evolution
 components Cosmological components (matter, radiation, Λ)
 friedmann Friedmann equations and solvers
 observations Observable quantities
 distances Distance measures (d_L, d_A, d_c)
 dark_matter Dark matter physics (future)
 cmb Cosmic Microwave Background
 recombination Ionization history
 fluctuations Angular power spectrum
 structure Large-scale structure
 power_spectrum Matter power spectrum P(k)
 transfer_function Transfer function T(k)
 perturbations Perturbation theory
 growth Linear growth factor
 boltzmann Boltzmann equations
 advanced Advanced features
 weak_lensing Convergence and shear
 polarization CMB polarization (future)
 visualization Plotting and visualization
 plots_2d Static PNG plots
 plotly_plots Interactive HTML plots
 three_d 3D visualization (future)
 colors Color schemes
```

---

## Constants

### Module: `andam::constants`

Physical and cosmological constants in SI units.

#### Physical Constants

```rust
pub const C: f64 = 299792458.0; // Speed of light [m/s]
pub const G: f64 = 6.67430e-11; // Gravitational constant [m³/kg/s²]
pub const H_PLANCK: f64 = 6.62607015e-34; // Planck constant [J⋅s]
pub const K_B: f64 = 1.380649e-23; // Boltzmann constant [J/K]
pub const SIGMA_SB: f64 = 5.670374419e-8; // Stefan-Boltzmann [W/m²/K⁴]
pub const M_P: f64 = 1.22089e19; // Planck mass [GeV/c²]
pub const M_SUN: f64 = 1.98847e30; // Solar mass [kg]
```

#### Cosmological Constants

```rust
pub const H0_PLANCK: f64 = 67.66; // Hubble constant [km/s/Mpc]
pub const OMEGA_M_PLANCK: f64 = 0.3111; // Matter density parameter
pub const OMEGA_LAMBDA_PLANCK: f64 = 0.6889; // Dark energy density
pub const OMEGA_B_PLANCK: f64 = 0.0490; // Baryon density
pub const T_CMB: f64 = 2.7255; // CMB temperature [K]
pub const RHO_CRIT_0: f64 = 1.8788e-26; // Critical density [kg/m³]
```

#### Unit Conversion Constants

```rust
pub const MPC_TO_M: f64 = 3.085677581e22; // Megaparsec to meters
pub const GYR_TO_S: f64 = 3.15576e16; // Gigayear to seconds
pub const EV_TO_J: f64 = 1.602176634e-19; // Electronvolt to joules
```

#### Functions

##### `hubble_time() -> f64`

Returns the Hubble time t_H = 1/H0 in seconds.

**Returns**: Hubble time in seconds (~4.55 × 10^17 s)

**Example**:
```rust
use andam::constants::hubble_time;

let t_h = hubble_time();
println!("Hubble time: {:.2e} s", t_h);
```

##### `critical_density() -> f64`

Returns the critical density ρ_c = 3H₀²/(8πG) in kg/m³.

**Returns**: Critical density in kg/m³ (~1.88 × 10^-26 kg/m³)

**Example**:
```rust
use andam::constants::critical_density;

let rho_c = critical_density();
println!("Critical density: {:.2e} kg/m³", rho_c);
```

---

## Units

### Module: `andam::units`

Unit conversion utilities for length, mass, time, and energy.

#### Length Enum

```rust
pub enum Length {
 Meter,
 Kilometer,
 Megaparsec,
 Lightyear,
}
```

**Methods**:

- `to_meters(value: f64) -> f64`: Convert to meters
- `to_km(value: f64) -> f64`: Convert to kilometers
- `to_mpc(value: f64) -> f64`: Convert to megaparsecs
- `to_lightyears(value: f64) -> f64`: Convert to light years

**Example**:
```rust
use andam::units::Length;

let mpc = 1.0;
let meters = Length::Mpc.to_meters(mpc);
let km = Length::Mpc.to_km(mpc);
let ly = Length::Mpc.to_lightyears(mpc);
```

#### Time Enum

```rust
pub enum Time {
 Second,
 Year,
 Megayear,
 Gigayear,
}
```

**Methods**:

- `to_seconds(value: f64) -> f64`: Convert to seconds
- `to_years(value: f64) -> f64`: Convert to years
- `to_myr(value: f64) -> f64`: Convert to megayears
- `to_gyr(value: f64) -> f64`: Convert to gigayears

**Example**:
```rust
use andam::units::Time;

let gyr = 13.8;
let seconds = Time::Gyr.to_seconds(gyr);
let years = Time::Gyr.to_years(gyr);
```

#### Temperature Functions

##### `temperature_to_ev(temp_kelvin: f64) -> f64`

Convert temperature in Kelvin to energy in electronvolts.

**Parameters**:
- `temp_kelvin`: Temperature in Kelvin

**Returns**: Energy in eV (E = k_B T)

**Example**:
```rust
use andam::units::temperature_to_ev;

let t_cmb = 2.7255; // K
let e_cmb = temperature_to_ev(t_cmb);
println!("CMB photon energy: {:.4e} eV", e_cmb);
```

##### `ev_to_temperature(energy_ev: f64) -> f64`

Convert energy in electronvolts to temperature in Kelvin.

**Parameters**:
- `energy_ev`: Energy in eV

**Returns**: Temperature in Kelvin (T = E / k_B)

**Example**:
```rust
use andam::units::ev_to_temperature;

let e = 13.6; // eV (ionization energy of hydrogen)
let t = ev_to_temperature(e);
println!("Temperature: {:.2e} K", t);
```

---

## Dynamics

### Module: `andam::dynamics`

Universe models, components, and evolution equations.

#### Struct: `Component`

Represents a cosmological component with equation of state.

```rust
pub struct Component {
 pub omega_0: f64, // Density parameter today
 pub w: f64, // Equation of state parameter
}
```

**Constructors**:

##### `Component::new(omega_0: f64, w: f64) -> Component`

Create custom component.

**Parameters**:
- `omega_0`: Density parameter at a=1
- `w`: Equation of state parameter (p = wρc²)

**Returns**: New Component instance

##### `Component::matter(omega_0: f64) -> Component`

Create matter component (w = 0).

**Parameters**:
- `omega_0`: Matter density parameter

**Returns**: Matter component with w = 0

##### `Component::radiation(omega_0: f64) -> Component`

Create radiation component (w = 1/3).

**Parameters**:
- `omega_0`: Radiation density parameter

**Returns**: Radiation component with w = 1/3

##### `Component::lambda(omega_0: f64) -> Component`

Create dark energy component (w = -1).

**Parameters**:
- `omega_0`: Dark energy density parameter

**Returns**: Dark energy component with w = -1

**Methods**:

##### `density(&self, a: f64) -> f64`

Calculate density parameter at scale factor a.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)

**Returns**: Density parameter Ω(a) = Ω₀ a^(-3(1+w))

**Example**:
```rust
use andam::dynamics::Component;

let matter = Component::matter(0.3);
let omega_a05 = matter.density(0.5);
println!("Ω_m(a=0.5) = {:.4}", omega_a05);
```

#### Struct: `Universe`

Complete cosmological model with multiple components.

```rust
pub struct Universe {
 pub h0: f64, // Hubble constant [km/s/Mpc]
 pub components: Vec<Component>, // List of components
}
```

**Constructors**:

##### `Universe::new(h0: f64) -> Universe`

Create empty universe with given Hubble constant.

**Parameters**:
- `h0`: Hubble constant in km/s/Mpc

**Returns**: Empty Universe (no components)

**Example**:
```rust
use andam::dynamics::Universe;

let mut universe = Universe::new(70.0);
```

##### `Universe::benchmark() -> Universe`

Create standard ΛCDM universe with Planck 2018 parameters.

**Returns**: Universe with:
- H₀ = 67.66 km/s/Mpc
- Ω_m = 0.3111
- Ω_r = 9.3×10^-5
- Ω_Λ = 0.6889

**Example**:
```rust
use andam::dynamics::Universe;

let universe = Universe::benchmark();
```

##### `Universe::einstein_de_sitter(h0: f64) -> Universe`

Create Einstein-de Sitter universe (matter only, Ω_m = 1).

**Parameters**:
- `h0`: Hubble constant in km/s/Mpc

**Returns**: Matter-only universe with Ω_m = 1.0

**Example**:
```rust
use andam::dynamics::Universe;

let eds = Universe::einstein_de_sitter(70.0);
```

**Methods**:

##### `add_component(&mut self, component: Component)`

Add a component to the universe.

**Parameters**:
- `component`: Component to add

**Example**:
```rust
use andam::dynamics::{Universe, Component};

let mut universe = Universe::new(70.0);
universe.add_component(Component::matter(0.3));
universe.add_component(Component::lambda(0.7));
```

##### `hubble(&self, a: f64) -> f64`

Calculate Hubble parameter H(a) in km/s/Mpc.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)

**Returns**: Hubble parameter in km/s/Mpc

**Formula**: H²(a) = H₀² Σ Ωᵢ a^(-3(1+wᵢ))

**Example**:
```rust
let universe = Universe::benchmark();
let h = universe.hubble(0.5);
println!("H(a=0.5) = {:.2} km/s/Mpc", h);
```

##### `hubble_normalized(&self, a: f64) -> f64`

Calculate dimensionless Hubble parameter E(a) = H(a)/H₀.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)

**Returns**: Dimensionless H(a)/H₀

**Example**:
```rust
let universe = Universe::benchmark();
let e = universe.hubble_normalized(0.5);
println!("E(a=0.5) = {:.4}", e);
```

##### `age(&self, a: f64) -> f64`

Calculate age of universe at scale factor a in Gyr.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)

**Returns**: Age in gigayears

**Formula**: t(a) = ∫₀ᵃ da'/(a' H(a'))

**Example**:
```rust
let universe = Universe::benchmark();
let age = universe.age(0.5);
println!("Age at a=0.5: {:.2} Gyr", age);
```

##### `age_today(&self) -> f64`

Calculate current age of universe (at a=1) in Gyr.

**Returns**: Current age in gigayears

**Example**:
```rust
let universe = Universe::benchmark();
let age = universe.age_today();
println!("Universe age: {:.2} Gyr", age);
```

##### `deceleration(&self, a: f64) -> f64`

Calculate deceleration parameter q(a).

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)

**Returns**: Deceleration parameter

**Formula**: q = -1 - d ln H / d ln a

**Example**:
```rust
let universe = Universe::benchmark();
let q = universe.deceleration(1.0);
println!("Deceleration parameter today: {:.3}", q);
```

##### `omega_total(&self) -> f64`

Calculate total density parameter today.

**Returns**: Ω_total = Σ Ωᵢ

**Example**:
```rust
let universe = Universe::benchmark();
let omega = universe.omega_total();
println!("Ω_total = {:.6}", omega); // Should be ~1.0
```

---

## Observations

### Module: `andam::observations`

Observable quantities including distance measures.

#### Distance Measures

##### `comoving_distance(z: f64, universe: &Universe) -> f64`

Calculate comoving distance to redshift z.

**Parameters**:
- `z`: Redshift (z ≥ 0)
- `universe`: Universe model

**Returns**: Comoving distance in Mpc

**Formula**: d_c = c ∫₀ᶻ dz'/H(z')

**Integration**: Simpson's rule with 1000 steps

**Example**:
```rust
use andam::observations::comoving_distance;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let d_c = comoving_distance(1.0, &universe);
println!("Comoving distance to z=1: {:.1} Mpc", d_c);
```

##### `luminosity_distance(z: f64, universe: &Universe) -> f64`

Calculate luminosity distance to redshift z.

**Parameters**:
- `z`: Redshift (z ≥ 0)
- `universe`: Universe model

**Returns**: Luminosity distance in Mpc

**Formula**: d_L = (1 + z) d_c

**Used for**: Standard candles (supernovae)

**Example**:
```rust
use andam::observations::luminosity_distance;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let d_l = luminosity_distance(1.0, &universe);

// Distance modulus
let mu = 5.0 * (d_l * 1e6 / 10.0).log10();
println!("Distance modulus: {:.2} mag", mu);
```

##### `angular_diameter_distance(z: f64, universe: &Universe) -> f64`

Calculate angular diameter distance to redshift z.

**Parameters**:
- `z`: Redshift (z ≥ 0)
- `universe`: Universe model

**Returns**: Angular diameter distance in Mpc

**Formula**: d_A = d_c / (1 + z)

**Used for**: Standard rulers (BAO, CMB)

**Relation**: d_L = (1 + z)² d_A (Etherington relation)

**Example**:
```rust
use andam::observations::angular_diameter_distance;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let z_cmb = 1090.0;
let d_a = angular_diameter_distance(z_cmb, &universe);
println!("Angular diameter distance to CMB: {:.1} Mpc", d_a);
```

---

## CMB

### Module: `andam::cmb`

Cosmic Microwave Background physics.

#### Recombination

##### `ionization_fraction(z: f64, universe: &Universe) -> f64`

Calculate ionization fraction using Saha equation.

**Parameters**:
- `z`: Redshift (z ≥ 0)
- `universe`: Universe model

**Returns**: Ionization fraction x_e (0 ≤ x_e ≤ 1)

**Formula**: Saha equation for hydrogen ionization

**Example**:
```rust
use andam::cmb::ionization_fraction;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let x_e = ionization_fraction(1100.0, &universe);
println!("Ionization at z=1100: {:.4}", x_e);
```

##### `recombination_redshift(universe: &Universe) -> f64`

Find redshift where ionization fraction = 0.5.

**Parameters**:
- `universe`: Universe model

**Returns**: Recombination redshift (typically ~1090)

**Method**: Bisection on ionization_fraction

**Example**:
```rust
use andam::cmb::recombination_redshift;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let z_rec = recombination_redshift(&universe);
println!("Recombination at z = {:.0}", z_rec);
```

#### Fluctuations

##### `angular_power_spectrum(l_max: usize, universe: &Universe) -> Vec<f64>`

Calculate CMB angular power spectrum C_ℓ.

**Parameters**:
- `l_max`: Maximum multipole
- `universe`: Universe model

**Returns**: Vector of C_ℓ values (length l_max + 1)

**Model**: Phenomenological with acoustic peaks

**Example**:
```rust
use andam::cmb::fluctuations::angular_power_spectrum;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let c_l = angular_power_spectrum(2000, &universe);
println!("C_ℓ at l=220: {:.4e}", c_l[220]);
```

##### `dimensionless_power_spectrum(c_l: &[f64]) -> Vec<(usize, f64)>`

Convert C_ℓ to dimensionless power D_ℓ = ℓ(ℓ+1)C_ℓ/(2π).

**Parameters**:
- `c_l`: Angular power spectrum

**Returns**: Vector of (ℓ, D_ℓ) pairs

**Used for**: Standard plotting convention

**Example**:
```rust
use andam::cmb::fluctuations::*;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let c_l = angular_power_spectrum(2000, &universe);
let d_l = dimensionless_power_spectrum(&c_l);

for (l, d) in d_l.iter().skip(2).take(5) {
 println!("l = {}: D_l = {:.2e}", l, d);
}
```

##### `acoustic_peak_positions(universe: &Universe) -> Vec<usize>`

Calculate positions of first 5 acoustic peaks.

**Parameters**:
- `universe`: Universe model

**Returns**: Vector of peak positions [ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅]

**Typical values**: [219, 438, 658, 877, 1096]

**Example**:
```rust
use andam::cmb::fluctuations::acoustic_peak_positions;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let peaks = acoustic_peak_positions(&universe);
println!("First peak at l = {}", peaks[0]);
```

---

## Structure

### Module: `andam::structure`

Large-scale structure and matter power spectrum.

#### Power Spectrum

##### `primordial_power_spectrum(k: f64, a_s: f64, n_s: f64) -> f64`

Calculate primordial power spectrum from inflation.

**Parameters**:
- `k`: Wavenumber in h/Mpc
- `a_s`: Scalar amplitude (typically 2.1×10^-9)
- `n_s`: Spectral index (typically 0.9665)

**Returns**: Primordial power P_R(k)

**Formula**: P_R(k) = A_s (k/k_pivot)^(n_s - 1)

**Example**:
```rust
use andam::structure::power_spectrum::primordial_power_spectrum;

let p_r = primordial_power_spectrum(0.05, 2.1e-9, 0.9665);
println!("Primordial power: {:.4e}", p_r);
```

##### `transfer_function_eh(k: f64, omega_m: f64, omega_b: f64, h: f64) -> f64`

Calculate Eisenstein-Hu transfer function.

**Parameters**:
- `k`: Wavenumber in h/Mpc
- `omega_m`: Total matter density parameter
- `omega_b`: Baryon density parameter
- `h`: Dimensionless Hubble constant (H₀/100)

**Returns**: Transfer function T(k)

**Reference**: Eisenstein & Hu (1998), ApJ 496:605

**Example**:
```rust
use andam::structure::power_spectrum::transfer_function_eh;

let t_k = transfer_function_eh(0.1, 0.3111, 0.049, 0.6766);
println!("Transfer function: {:.4}", t_k);
```

##### `matter_power_spectrum(k: f64, z: f64, omega_m: f64, omega_b: f64, h: f64, a_s: f64, n_s: f64) -> f64`

Calculate matter power spectrum P(k,z).

**Parameters**:
- `k`: Wavenumber in h/Mpc
- `z`: Redshift
- `omega_m`: Total matter density
- `omega_b`: Baryon density
- `h`: Dimensionless Hubble constant
- `a_s`: Scalar amplitude
- `n_s`: Spectral index

**Returns**: Matter power spectrum in (Mpc/h)³

**Formula**: P(k,z) = P_R(k) T²(k) D²(z)

**Example**:
```rust
use andam::structure::power_spectrum::matter_power_spectrum;

let p_k = matter_power_spectrum(
 0.1, // k in h/Mpc
 0.0, // z
 0.3111, // omega_m
 0.049, // omega_b
 0.6766, // h
 2.1e-9, // A_s
 0.9665 // n_s
);
println!("P(k=0.1) = {:.2e} (Mpc/h)³", p_k);
```

##### `dimensionless_power(k: f64, z: f64, omega_m: f64, omega_b: f64, h: f64, a_s: f64, n_s: f64) -> f64`

Calculate dimensionless power Δ²(k) = k³P(k)/(2π²).

**Parameters**: Same as matter_power_spectrum

**Returns**: Dimensionless power Δ²(k)

**Used for**: Standard plotting convention

**Example**:
```rust
use andam::structure::power_spectrum::dimensionless_power;

let delta_sq = dimensionless_power(
 0.1, 0.0, 0.3111, 0.049, 0.6766, 2.1e-9, 0.9665
);
println!("Δ²(k=0.1) = {:.4}", delta_sq);
```

---

## Perturbations

### Module: `andam::perturbations`

Linear perturbation theory and structure growth.

#### Growth Factor

##### `growth_factor(a: f64, universe: &Universe) -> f64`

Calculate linear growth factor D(a), normalized to D(1) = 1.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)
- `universe`: Universe model

**Returns**: Linear growth factor

**Relation**: δ(a,k) = D(a) δ_init(k) on large scales

**Example**:
```rust
use andam::perturbations::growth::growth_factor;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let d = growth_factor(0.5, &universe);
println!("Growth factor at a=0.5: {:.4}", d);
```

##### `growth_rate(a: f64, universe: &Universe) -> f64`

Calculate logarithmic growth rate f = d ln D / d ln a.

**Parameters**:
- `a`: Scale factor (0 < a ≤ 1)
- `universe`: Universe model

**Returns**: Growth rate f(a)

**Approximation**: f ≈ Ω_m(a)^0.55 (Peebles)

**Example**:
```rust
use andam::perturbations::growth::growth_rate;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let f = growth_rate(1.0, &universe);
println!("Growth rate today: {:.3}", f);
```

#### Boltzmann Solver

##### Struct: `BoltzmannSolver`

Simplified Boltzmann equation solver for photon perturbations.

```rust
pub struct BoltzmannSolver {
 pub universe: Universe,
 pub k: f64,
 pub mode: PerturbationMode,
 pub l_max: usize,
}
```

##### `BoltzmannSolver::new(universe: Universe, k: f64, mode: PerturbationMode, l_max: usize) -> Self`

Create new Boltzmann solver.

**Parameters**:
- `universe`: Universe model
- `k`: Wavenumber in h/Mpc
- `mode`: Scalar or Vector
- `l_max`: Maximum multipole

**Returns**: New solver instance

**Example**:
```rust
use andam::perturbations::*;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let solver = BoltzmannSolver::new(
 universe,
 0.1,
 PerturbationMode::Scalar,
 10
);
```

##### `evolve(&mut self, a_initial: f64, a_final: f64) -> Vec<PerturbationState>`

Evolve perturbations from a_initial to a_final.

**Parameters**:
- `a_initial`: Starting scale factor
- `a_final`: Ending scale factor

**Returns**: Vector of states at each timestep

**Method**: Euler integration

**Example**:
```rust
let mut solver = BoltzmannSolver::new(
 universe, 0.1, PerturbationMode::Scalar, 10
);
let results = solver.evolve(1e-5, 1e-2);
println!("Evolved {} steps", results.len());
```

---

## Advanced

### Module: `andam::advanced`

Advanced features including weak lensing.

#### Weak Lensing

##### Struct: `ConvergenceField`

2D convergence field κ(θ₁, θ₂).

```rust
pub struct ConvergenceField {
 pub kappa: Array2<f64>,
 pub size: usize,
}
```

##### `ConvergenceField::new(size: usize) -> Self`

Create zero convergence field.

**Parameters**:
- `size`: Grid size (size × size)

**Returns**: Zero-initialized field

**Example**:
```rust
use andam::advanced::weak_lensing::ConvergenceField;

let field = ConvergenceField::new(128);
println!("Created {}x{} field", field.size, field.size);
```

##### `ConvergenceField::from_power_spectrum(universe: &Universe, source_z: f64, size: usize, field_size_deg: f64) -> Self`

Generate convergence field from power spectrum.

**Parameters**:
- `universe`: Universe model
- `source_z`: Source redshift
- `size`: Grid size
- `field_size_deg`: Field size in degrees

**Returns**: Convergence field

**Note**: Current implementation uses simplified model

**Example**:
```rust
use andam::advanced::weak_lensing::ConvergenceField;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let field = ConvergenceField::from_power_spectrum(
 &universe,
 1.0, // source at z=1
 128, // 128x128 grid
 1.0 // 1 degree field
);
```

##### Struct: `Shear`

Shear components (γ₁, γ₂).

```rust
pub struct Shear {
 pub gamma1: f64,
 pub gamma2: f64,
}
```

**Methods**:

- `magnitude(&self) -> f64`: |γ| = √(γ₁² + γ₂²)
- `angle(&self) -> f64`: φ = arctan(γ₂/γ₁) / 2

**Example**:
```rust
use andam::advanced::weak_lensing::Shear;

let shear = Shear { gamma1: 0.03, gamma2: 0.04 };
println!("Shear: |γ| = {:.4}", shear.magnitude());
println!("Angle: φ = {:.2} rad", shear.angle());
```

##### `convergence_to_shear(kappa: &Array2<f64>) -> (Array2<f64>, Array2<f64>)`

Convert convergence to shear using Kaiser-Squires.

**Parameters**:
- `kappa`: Convergence field

**Returns**: Tuple of (γ₁ field, γ₂ field)

**Method**: Finite differences for second derivatives

**Example**:
```rust
use andam::advanced::weak_lensing::*;

let field = ConvergenceField::new(128);
let (gamma1, gamma2) = convergence_to_shear(&field.kappa);
```

##### `lensing_power_spectrum(l: usize, source_z: f64, universe: &Universe) -> f64`

Calculate lensing convergence power spectrum C_ℓ^κκ.

**Parameters**:
- `l`: Multipole
- `source_z`: Source redshift
- `universe`: Universe model

**Returns**: Lensing power at multipole ℓ

**Method**: Limber approximation

**Example**:
```rust
use andam::advanced::weak_lensing::lensing_power_spectrum;
use andam::dynamics::Universe;

let universe = Universe::benchmark();
let c_l = lensing_power_spectrum(1000, 1.0, &universe);
println!("C_l^κκ at l=1000: {:.4e}", c_l);
```

---

## Visualization

### Module: `andam::visualization`

Plotting and visualization utilities.

#### Static 2D Plots

Module: `andam::visualization::plots_2d`

##### `create_line_plot(data: &[(f64, f64)], filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>>`

Create static line plot as PNG.

**Parameters**:
- `data`: Vector of (x, y) points
- `filename`: Output filename (e.g., "plot.png")
- `title`: Plot title
- `xlabel`: X-axis label
- `ylabel`: Y-axis label

**Returns**: Result indicating success or error

**Example**:
```rust
use andam::visualization::plots_2d::create_line_plot;

let data = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
create_line_plot(
 &data,
 "plot.png",
 "Quadratic Function",
 "x",
 "y = x²"
)?;
```

##### `create_loglog_plot(data: &[(f64, f64)], filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>>`

Create log-log plot.

**Parameters**: Same as create_line_plot

**Note**: Both axes use logarithmic scale

**Example**:
```rust
use andam::visualization::plots_2d::create_loglog_plot;

let data = vec![(1.0, 1.0), (10.0, 100.0), (100.0, 10000.0)];
create_loglog_plot(&data, "loglog.png", "Power Law", "x", "y")?;
```

##### `create_multiline_plot(series: &[Vec<(f64, f64)>], labels: &[&str], filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>>`

Create plot with multiple series.

**Parameters**:
- `series`: Vector of data series
- `labels`: Labels for each series
- `filename`: Output filename
- `title`: Plot title
- `xlabel`: X-axis label
- `ylabel`: Y-axis label

**Returns**: Result indicating success or error

**Example**:
```rust
use andam::visualization::plots_2d::create_multiline_plot;

let s1 = vec![(0.0, 0.0), (1.0, 1.0)];
let s2 = vec![(0.0, 0.0), (1.0, 2.0)];
create_multiline_plot(
 &[s1, s2],
 &["Linear", "Double"],
 "multiline.png",
 "Comparison",
 "x",
 "y"
)?;
```

#### Interactive HTML Plots

Module: `andam::visualization::plotly_plots`

##### `create_interactive_plot(data: &[(f64, f64)], filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>>`

Create interactive HTML plot with Plotly.

**Parameters**: Same as create_line_plot

**Returns**: Result indicating success or error

**Features**: Zoom, pan, hover tooltips

**Example**:
```rust
use andam::visualization::plotly_plots::create_interactive_plot;

let data = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
create_interactive_plot(
 &data,
 "interactive.html",
 "Interactive Plot",
 "x",
 "y"
)?;
```

##### `create_loglog_interactive(data: &[(f64, f64)], filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>>`

Create interactive log-log plot.

**Parameters**: Same as create_line_plot

**Returns**: Result indicating success or error

**Example**:
```rust
use andam::visualization::plotly_plots::create_loglog_interactive;

let data = vec![(1.0, 1.0), (10.0, 100.0)];
create_loglog_interactive(
 &data,
 "interactive_log.html",
 "Log Plot",
 "x",
 "y"
)?;
```

---

## Type Aliases and Enums

### PerturbationMode

```rust
pub enum PerturbationMode {
 Scalar,
 Vector,
 Tensor,
}
```

Modes of cosmological perturbations:
- **Scalar**: Density perturbations (most important)
- **Vector**: Rotational modes (decay rapidly)
- **Tensor**: Gravitational waves

### PerturbationState

```rust
pub struct PerturbationState {
 pub a: f64,
 pub delta: f64,
 pub theta: Vec<f64>,
}
```

State of perturbations at given scale factor:
- `a`: Scale factor
- `delta`: Density perturbation
- `theta`: Photon multipoles [Θ₀, Θ₁, Θ₂, ...]

---

## Error Handling

Most functions return concrete types rather than Results. Invalid inputs may produce:

- `NaN`: For mathematically undefined operations
- `Inf`: For overflow conditions
- Panic: For out-of-bounds array access

Best practice: Validate inputs before calling functions.

---

## Performance Characteristics

### Computational Complexity

- `Universe::hubble()`: O(n) where n = number of components (typically 3)
- `Universe::age()`: O(N) where N = integration steps (default 1000)
- `comoving_distance()`: O(N) where N = integration steps (default 1000)
- `angular_power_spectrum()`: O(l_max)
- `matter_power_spectrum()`: O(1)
- `growth_factor()`: O(N) where N = integration steps (default 1000)

### Memory Usage

All calculations use stack-allocated variables except:
- `Vec` allocations for results
- `Array2` for convergence fields

Typical memory footprint: < 1 MB for all calculations.

---

## Version Compatibility

This API documentation is for version 0.1.0.

Future versions may introduce breaking changes with appropriate deprecation warnings.

---

## See Also

- [User Guide](USER_GUIDE.md) - Comprehensive usage guide
- [Project Overview](PROJECT_OVERVIEW.md) - Project structure and timeline
- [Examples](../examples/) - Complete example programs
- Rust documentation: `cargo doc --open`

---

## References

1. Ryden, B. (2016). Introduction to Cosmology (2nd ed.)
2. Dodelson, S. (2003). Modern Cosmology
3. Planck Collaboration (2018). Planck 2018 results
4. Eisenstein & Hu (1998). ApJ 496:605
