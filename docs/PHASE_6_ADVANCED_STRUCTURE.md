# Phase 6: Advanced Structure Formation (Weeks 18-20)

## Overview
Implement advanced structure formation including non-linear power spectrum (HALOFIT), halo mass function, galaxy bias, redshift-space distortions, and N-body simulation visualization. Features 3D cosmic web visualizations and high-resolution plots with equations.

---

## Prerequisites
✅ Phases 1-5 completed
✅ Linear power spectrum working
✅ Growth factor implemented

---

## Objectives
- [x] Non-linear power spectrum (HALOFIT)
- [x] Halo mass function (Press-Schechter, Sheth-Tormen)
- [x] Galaxy bias models
- [x] Redshift-space distortions
- [x] Correlation function ξ(r)
- [x] 3D cosmic web visualization
- [x] N-body simulation framework
- [x] Publication-quality plots with equations

---

## Task 1: Non-Linear Power Spectrum

### Step 1.1: Create `src/structure/nonlinear.rs`

```rust
//! Non-linear power spectrum calculations

use crate::structure::power_spectrum::*;
use std::f64::consts::PI;

/// HALOFIT non-linear power spectrum (Smith et al. 2003, Takahashi et al. 2012)
pub struct HalofitSpectrum {
    /// Linear power spectrum calculator
    pub omega_m: f64,
    pub omega_b: f64,
    pub h: f64,
    pub sigma_8: f64,
    pub n_s: f64,
}

impl HalofitSpectrum {
    /// Create new HALOFIT calculator
    pub fn new(omega_m: f64, omega_b: f64, h: f64, sigma_8: f64, n_s: f64) -> Self {
        HalofitSpectrum {
            omega_m,
            omega_b,
            h,
            sigma_8,
            n_s,
        }
    }
    
    /// Non-linear scale k_nl at redshift z
    pub fn k_nonlinear(&self, z: f64) -> f64 {
        // Find k where Δ²(k) = 1
        // Binary search
        let mut k_min = 0.001;
        let mut k_max = 10.0;
        
        for _ in 0..50 {
            let k_mid = (k_min + k_max) / 2.0;
            let delta2 = dimensionless_power(k_mid, z, self.omega_m, self.omega_b, 
                                           self.h, 2.1e-9, self.n_s);
            
            if delta2 > 1.0 {
                k_max = k_mid;
            } else {
                k_min = k_mid;
            }
        }
        
        (k_min + k_max) / 2.0
    }
    
    /// Effective spectral index at scale k
    pub fn n_eff(&self, k: f64, z: f64) -> f64 {
        let dk = k * 0.01;
        let p1 = self.linear_power(k - dk, z);
        let p2 = self.linear_power(k + dk, z);
        
        (p2 / p1).ln() / (2.0 * dk / k).ln()
    }
    
    /// Curvature of power spectrum
    pub fn n_curv(&self, k: f64, z: f64) -> f64 {
        let dk = k * 0.01;
        let n1 = self.n_eff(k - dk, z);
        let n2 = self.n_eff(k + dk, z);
        
        (n2 - n1) / (2.0 * dk / k)
    }
    
    /// Linear power spectrum
    fn linear_power(&self, k: f64, z: f64) -> f64 {
        matter_power_spectrum(k, z, self.omega_m, self.omega_b, 
                            self.h, 2.1e-9, self.n_s)
    }
    
    /// HALOFIT non-linear power spectrum
    pub fn nonlinear_power(&self, k: f64, z: f64) -> f64 {
        let omega_m_z = self.omega_m * (1.0 + z).powi(3);
        let omega_de_z = 1.0 - omega_m_z;
        
        let p_lin = self.linear_power(k, z);
        let k_nl = self.k_nonlinear(z);
        
        let n_eff = self.n_eff(k, z);
        let n_curv = self.n_curv(k, z);
        
        // HALOFIT parameters (Takahashi et al. 2012)
        let a_n = 10.0_f64.powf(1.5222 + 2.8553 * n_eff + 2.3706 * n_eff.powi(2)
                              + 0.9903 * n_eff.powi(3) + 0.2250 * n_eff.powi(4)
                              - 0.6038 * n_curv);
        
        let b_n = 10.0_f64.powf(-0.5642 + 0.5864 * n_eff + 0.5716 * n_eff.powi(2)
                              - 1.5474 * n_curv);
        
        let c_n = 10.0_f64.powf(0.3698 + 2.0404 * n_eff + 0.8161 * n_eff.powi(2)
                              + 0.5869 * n_curv);
        
        let gamma_n = 0.1971 - 0.0843 * n_eff + 0.8460 * n_curv;
        let alpha_n = (3.5634 + 2.3313 * n_eff + 0.9847 * n_eff.powi(2)
                     + 0.2599 * n_eff.powi(3) - 0.2729 * n_curv).abs();
        let beta_n = (0.8291 + 0.9854 * n_eff + 0.3400 * n_eff.powi(2)).abs();
        let mu_n = 0.0;
        let nu_n = 10.0_f64.powf(5.2105 + 3.6902 * n_eff);
        
        let f1 = omega_m_z.powf(-0.0307);
        let f2 = omega_m_z.powf(-0.0585);
        let f3 = omega_m_z.powf(0.0743);
        
        let y = k / k_nl;
        
        // Quasi-linear contribution
        let delta_q = p_lin * (1.0 + p_lin).powf(beta_n / (1.0 + alpha_n * p_lin))
                    * (-y / 4.0 - y.powi(2) / 8.0).exp();
        
        // Halo term
        let delta_h = a_n * y.powf(3.0 * f1) 
                    / (1.0 + b_n * y.powf(f2) + (c_n * f3 * y).powf(3.0 - gamma_n));
        
        let delta_h = delta_h / (1.0 + mu_n / y + nu_n / y.powi(2));
        
        // Total non-linear power
        p_lin * (1.0 + delta_q).powi(2) * (1.0 + delta_h)
    }
    
    /// Ratio of non-linear to linear power
    pub fn boost_factor(&self, k: f64, z: f64) -> f64 {
        let p_nl = self.nonlinear_power(k, z);
        let p_lin = self.linear_power(k, z);
        (p_nl / p_lin).sqrt()
    }
}

/// Dimensionless non-linear power spectrum
pub fn nonlinear_dimensionless_power(
    k: f64,
    z: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
    sigma_8: f64,
    n_s: f64,
) -> f64 {
    let halofit = HalofitSpectrum::new(omega_m, omega_b, h, sigma_8, n_s);
    let p_nl = halofit.nonlinear_power(k, z);
    k.powi(3) * p_nl / (2.0 * PI.powi(2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_halofit_k_nonlinear() {
        let halofit = HalofitSpectrum::new(0.3, 0.05, 0.7, 0.8, 0.96);
        let k_nl = halofit.k_nonlinear(0.0);
        
        // Should be around 0.1-1 h/Mpc
        assert!(k_nl > 0.05 && k_nl < 5.0);
    }
    
    #[test]
    fn test_boost_factor_increases_at_small_scales() {
        let halofit = HalofitSpectrum::new(0.3, 0.05, 0.7, 0.8, 0.96);
        
        let boost_large = halofit.boost_factor(0.01, 0.0);
        let boost_small = halofit.boost_factor(1.0, 0.0);
        
        // Non-linear effects stronger at small scales
        assert!(boost_small > boost_large);
    }
}
```

---

## Task 2: Halo Mass Function

### Step 2.1: Create `src/structure/halos.rs`

```rust
//! Halo mass function and bias

use crate::constants::*;
use crate::structure::power_spectrum::*;
use std::f64::consts::PI;

/// Halo mass function formalism
#[derive(Debug, Clone, Copy)]
pub enum MassFunctionType {
    PressSchechter,
    ShethTormen,
    Tinker,
}

/// Halo mass function calculator
pub struct HaloMassFunction {
    pub omega_m: f64,
    pub sigma_8: f64,
    pub mf_type: MassFunctionType,
}

impl HaloMassFunction {
    /// Create new mass function calculator
    pub fn new(omega_m: f64, sigma_8: f64, mf_type: MassFunctionType) -> Self {
        HaloMassFunction {
            omega_m,
            sigma_8,
            mf_type,
        }
    }
    
    /// RMS fluctuation σ(M) at mass scale M
    pub fn sigma_mass(&self, mass_msun: f64, z: f64) -> f64 {
        // R = (3M / 4πρ̄)^(1/3)
        let rho_crit = critical_density(70.0); // kg/m³
        let rho_m = self.omega_m * rho_crit;
        let rho_m_msun_mpc3 = rho_m * (PARSEC * 1e6).powi(3) / M_SUN;
        
        let r_mpc = ((3.0 * mass_msun) / (4.0 * PI * rho_m_msun_mpc3)).powf(1.0/3.0);
        
        // σ²(R) = ∫ dk/k P(k) W²(kR) k³/(2π²)
        // Using approximate scaling: σ(M) ∝ M^(-0.5)
        let m_pivot = 1e12; // M_sun
        let sigma_pivot = 1.0;
        
        sigma_pivot * (mass_msun / m_pivot).powf(-0.5) / (1.0 + z).powf(0.5)
    }
    
    /// Peak height ν = δ_c / σ(M)
    pub fn peak_height(&self, mass_msun: f64, z: f64) -> f64 {
        let delta_c = 1.686; // Critical overdensity for collapse
        let sigma = self.sigma_mass(mass_msun, z);
        delta_c / sigma
    }
    
    /// Multiplicity function f(ν)
    pub fn multiplicity_function(&self, nu: f64) -> f64 {
        match self.mf_type {
            MassFunctionType::PressSchechter => {
                // f(ν) = √(2/π) ν exp(-ν²/2)
                (2.0 / PI).sqrt() * nu * (-nu.powi(2) / 2.0).exp()
            },
            MassFunctionType::ShethTormen => {
                // Sheth-Tormen (1999)
                let a = 0.707;
                let p = 0.3;
                let a_nu2 = a * nu.powi(2);
                
                let a_factor = (2.0 * a / PI).sqrt();
                a_factor * (1.0 + a_nu2.powf(-p)) * (a_nu2).sqrt() 
                    * (-a_nu2 / 2.0).exp()
            },
            MassFunctionType::Tinker => {
                // Tinker et al. (2008) - simplified
                let alpha = 0.368;
                let beta = 0.589;
                let gamma = 0.864;
                let phi = -0.729;
                
                alpha * (1.0 + (beta / nu).powf(2.0 * phi)) 
                    * nu.powf(2.0 * beta) * (-gamma * nu.powi(2) / 2.0).exp()
            },
        }
    }
    
    /// Halo mass function dn/dM [h⁴ Mpc⁻³ M_sun⁻¹]
    pub fn dn_dm(&self, mass_msun: f64, z: f64) -> f64 {
        let rho_crit = critical_density(70.0);
        let rho_m = self.omega_m * rho_crit;
        let rho_m_msun_mpc3 = rho_m * (PARSEC * 1e6).powi(3) / M_SUN;
        
        let nu = self.peak_height(mass_msun, z);
        let f_nu = self.multiplicity_function(nu);
        
        let sigma = self.sigma_mass(mass_msun, z);
        let dlnσ_dlnM = 0.5; // Approximate
        
        (rho_m_msun_mpc3 / mass_msun) * f_nu * (dlnσ_dlnM / sigma)
    }
    
    /// Number density of halos above mass M_min
    pub fn number_density(&self, m_min: f64, m_max: f64, z: f64, n_bins: usize) -> f64 {
        let log_m_min = m_min.ln();
        let log_m_max = m_max.ln();
        let dlogm = (log_m_max - log_m_min) / n_bins as f64;
        
        let mut n_total = 0.0;
        
        for i in 0..n_bins {
            let log_m = log_m_min + (i as f64 + 0.5) * dlogm;
            let m = log_m.exp();
            let dn_dm = self.dn_dm(m, z);
            n_total += dn_dm * m * dlogm;
        }
        
        n_total
    }
    
    /// Halo bias b(M)
    pub fn halo_bias(&self, mass_msun: f64, z: f64) -> f64 {
        let nu = self.peak_height(mass_msun, z);
        let delta_c = 1.686;
        
        match self.mf_type {
            MassFunctionType::PressSchechter => {
                // b = 1 + (ν² - 1)/δ_c
                1.0 + (nu.powi(2) - 1.0) / delta_c
            },
            MassFunctionType::ShethTormen => {
                // Sheth-Tormen bias
                let a = 0.707;
                let b_param = 0.5;
                let c = 0.6;
                
                let a_nu2 = a * nu.powi(2);
                1.0 + (a_nu2 - 1.0) / delta_c 
                    + 2.0 * b_param / (delta_c * (1.0 + a_nu2.powf(c)))
            },
            MassFunctionType::Tinker => {
                // Tinker et al. (2010) bias
                let y = nu.powi(2).ln();
                let a = 1.0 + 0.24 * y * (-4.0_f64.powf(4.0));
                let b_t = 0.183;
                let c = 0.019 + 0.107 * y + 0.19 * (-4.0_f64.powf(4.0)).exp();
                
                1.0 - a * nu.powf(b_t) / (nu.powf(b_t) + delta_c.powf(b_t)) + c
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_halo_mass_function() {
        let hmf = HaloMassFunction::new(0.3, 0.8, MassFunctionType::ShethTormen);
        
        let m = 1e14; // M_sun
        let z = 0.0;
        
        let dn_dm = hmf.dn_dm(m, z);
        assert!(dn_dm > 0.0);
    }
    
    #[test]
    fn test_halo_bias_increases_with_mass() {
        let hmf = HaloMassFunction::new(0.3, 0.8, MassFunctionType::ShethTormen);
        
        let b1 = hmf.halo_bias(1e12, 0.0);
        let b2 = hmf.halo_bias(1e15, 0.0);
        
        assert!(b2 > b1);
    }
}
```

---

## Task 3: Correlation Function and RSD

### Step 3.1: Create `src/structure/correlation.rs`

```rust
//! Two-point correlation function and redshift-space distortions

use crate::structure::power_spectrum::*;
use std::f64::consts::PI;

/// Real-space correlation function ξ(r)
pub fn correlation_function(r_mpc: f64, z: f64, omega_m: f64) -> f64 {
    // ξ(r) = ∫ dk/(2π²) k² P(k) j₀(kr)
    // where j₀(x) = sin(x)/x is spherical Bessel function
    
    let n_k = 1000;
    let k_min = 0.001;
    let k_max = 10.0;
    let dk = (k_max.ln() - k_min.ln()) / n_k as f64;
    
    let mut xi = 0.0;
    
    for i in 0..n_k {
        let log_k = k_min.ln() + (i as f64 + 0.5) * dk;
        let k = log_k.exp();
        
        let p_k = matter_power_spectrum(k, z, omega_m, 0.05, 0.7, 2.1e-9, 0.96);
        let kr = k * r_mpc;
        let j0 = if kr < 1e-3 {
            1.0 - kr.powi(2) / 6.0
        } else {
            kr.sin() / kr
        };
        
        xi += k.powi(2) * p_k * j0 * k * dk;
    }
    
    xi / (2.0 * PI.powi(2))
}

/// Redshift-space distortion parameter β = f/b
pub fn rsd_parameter(z: f64, omega_m: f64, bias: f64) -> f64 {
    // Linear growth rate f ≈ Ω_m^0.55
    let f = omega_m.powf(0.55);
    f / bias
}

/// Redshift-space correlation function ξ_s(r_p, π)
pub fn correlation_function_rsd(
    r_parallel: f64,  // π (parallel to line of sight)
    r_perp: f64,      // r_p (perpendicular to line of sight)
    z: f64,
    omega_m: f64,
    beta: f64,
) -> f64 {
    let r = (r_parallel.powi(2) + r_perp.powi(2)).sqrt();
    let mu = r_parallel / r;
    
    let xi_real = correlation_function(r, z, omega_m);
    
    // Kaiser formula (linear theory)
    xi_real * (1.0 + 2.0 * beta * mu.powi(2) / 3.0 + beta.powi(2) * mu.powi(4) / 5.0)
}

/// Monopole of correlation function
pub fn correlation_monopole(r_mpc: f64, z: f64, omega_m: f64) -> f64 {
    correlation_function(r_mpc, z, omega_m)
}

/// Quadrupole of correlation function
pub fn correlation_quadrupole(r_mpc: f64, z: f64, omega_m: f64, beta: f64) -> f64 {
    // ξ₂(r) = (4β/3 + 4β²/7) ξ(r)
    let xi_real = correlation_function(r_mpc, z, omega_m);
    (4.0 * beta / 3.0 + 4.0 * beta.powi(2) / 7.0) * xi_real
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_correlation_function_positive_at_small_scales() {
        let xi = correlation_function(5.0, 0.0, 0.3);
        assert!(xi > 0.0);
    }
    
    #[test]
    fn test_correlation_decreases_with_distance() {
        let xi1 = correlation_function(5.0, 0.0, 0.3);
        let xi2 = correlation_function(50.0, 0.0, 0.3);
        
        assert!(xi1.abs() > xi2.abs());
    }
}
```

---

## Task 4: High-Resolution Visualizations

### Step 4.1: Create `examples/structure/power_spectrum_comparison.rs`

```rust
//! Compare linear and non-linear power spectra with equations

use andam::structure::nonlinear::*;
use andam::structure::power_spectrum::*;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Creating power spectrum comparison...");
    
    let filename = "power_spectrum_comparison_hires.png";
    let root = BitMapBackend::new(filename, (2400, 1600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let omega_m = 0.3;
    let omega_b = 0.05;
    let h = 0.7;
    let sigma_8 = 0.8;
    let n_s = 0.96;
    
    let redshifts = vec![0.0, 0.5, 1.0, 2.0];
    let colors = vec![BLUE, RED, GREEN, MAGENTA];
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Linear vs Non-Linear Matter Power Spectrum", 
                ("DejaVu Sans", 48))
        .margin(20)
        .x_label_area_size(80)
        .y_label_area_size(100)
        .build_cartesian_2d(
            (1e-3..10.0).log_scale(),
            (1e-4..100.0).log_scale(),
        )?;
    
    chart.configure_mesh()
        .x_desc("Wavenumber k [h/Mpc]")
        .y_desc("Δ²(k) = k³P(k)/(2π²)")
        .x_label_style(("DejaVu Sans", 36))
        .y_label_style(("DejaVu Sans", 36))
        .draw()?;
    
    let halofit = HalofitSpectrum::new(omega_m, omega_b, h, sigma_8, n_s);
    
    for (i, (&z, &color)) in redshifts.iter().zip(colors.iter()).enumerate() {
        // Linear power
        let mut data_linear = Vec::new();
        let mut data_nonlinear = Vec::new();
        
        for j in 0..200 {
            let log_k = -3.0 + (j as f64) * 0.025;
            let k = 10_f64.powf(log_k);
            
            let delta2_lin = dimensionless_power(k, z, omega_m, omega_b, h, 2.1e-9, n_s);
            let delta2_nl = nonlinear_dimensionless_power(k, z, omega_m, omega_b, h, sigma_8, n_s);
            
            data_linear.push((k, delta2_lin));
            data_nonlinear.push((k, delta2_nl));
        }
        
        // Draw linear (dashed)
        chart.draw_series(LineSeries::new(
            data_linear.iter().map(|&(k, d)| (k, d)),
            color.stroke_width(2).dash(),
        ))?;
        
        // Draw non-linear (solid)
        chart.draw_series(LineSeries::new(
            data_nonlinear.iter().map(|&(k, d)| (k, d)),
            color.stroke_width(3),
        ))?
        .label(format!("z = {:.1}", z))
        .legend(move |(x, y)| {
            PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(3))
        });
    }
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.9))
        .border_style(&BLACK)
        .label_font(("DejaVu Sans", 32))
        .draw()?;
    
    // Add equation annotations
    chart.draw_series(std::iter::once(Text::new(
        "Linear: Δ²(k,z) = D²(z) T²(k) A_s (k/k_pivot)^(n_s-1)",
        (0.002, 50.0),
        ("DejaVu Sans", 28).into_font(),
    )))?;
    
    chart.draw_series(std::iter::once(Text::new(
        "Non-linear: HALOFIT (Takahashi et al. 2012)",
        (0.002, 30.0),
        ("DejaVu Sans", 28).into_font(),
    )))?;
    
    // Add k_nl marker
    let k_nl = halofit.k_nonlinear(0.0);
    chart.draw_series(LineSeries::new(
        vec![(k_nl, 1e-4), (k_nl, 100.0)],
        BLACK.stroke_width(2).dash(),
    ))?;
    
    chart.draw_series(std::iter::once(Text::new(
        format!("k_nl = {:.3} h/Mpc", k_nl),
        (k_nl * 1.2, 1.0),
        ("DejaVu Sans", 24).into_font(),
    )))?;
    
    root.present()?;
    println!("Created: {}", filename);
    
    Ok(())
}
```

### Step 4.2: Create `examples/structure/halo_mass_function.rs`

```rust
//! Halo mass function comparison with equations

use andam::structure::halos::*;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Creating halo mass function plot...");
    
    let filename = "halo_mass_function_hires.png";
    let root = BitMapBackend::new(filename, (2400, 1600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let omega_m = 0.3;
    let sigma_8 = 0.8;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Halo Mass Function: Theory Comparison", 
                ("DejaVu Sans", 48))
        .margin(20)
        .x_label_area_size(80)
        .y_label_area_size(120)
        .build_cartesian_2d(
            (1e11..1e16).log_scale(),
            (1e-20..1e-8).log_scale(),
        )?;
    
    chart.configure_mesh()
        .x_desc("Halo Mass M [M_☉/h]")
        .y_desc("dn/dM [h⁴ Mpc⁻³ (M_☉/h)⁻¹]")
        .x_label_style(("DejaVu Sans", 36))
        .y_label_style(("DejaVu Sans", 36))
        .draw()?;
    
    let models = vec![
        (MassFunctionType::PressSchechter, "Press-Schechter", BLUE),
        (MassFunctionType::ShethTormen, "Sheth-Tormen", RED),
        (MassFunctionType::Tinker, "Tinker et al.", GREEN),
    ];
    
    let z = 0.0;
    
    for (mf_type, label, color) in models {
        let hmf = HaloMassFunction::new(omega_m, sigma_8, mf_type);
        
        let mut data = Vec::new();
        
        for i in 0..100 {
            let log_m = 11.0 + (i as f64) * 0.05;
            let m = 10_f64.powf(log_m);
            
            let dn_dm = hmf.dn_dm(m, z);
            data.push((m, dn_dm));
        }
        
        chart.draw_series(LineSeries::new(
            data.iter().map(|&(m, n)| (m, n)),
            color.stroke_width(4),
        ))?
        .label(label)
        .legend(move |(x, y)| {
            PathElement::new(vec![(x, y), (x + 30, y)], color.stroke_width(4))
        });
    }
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.9))
        .border_style(&BLACK)
        .label_font(("DejaVu Sans", 32))
        .draw()?;
    
    // Add equations
    let equations = vec![
        ("Press-Schechter:", "f(ν) = √(2/π) ν exp(-ν²/2)", 3e15, 1e-10),
        ("Sheth-Tormen:", "f(ν) = A(1+(aν²)^(-p))√(aν²)exp(-aν²/2)", 3e15, 5e-11),
        ("ν = δ_c/σ(M),  δ_c = 1.686", "", 3e15, 2e-11),
    ];
    
    for (label, formula, x, y) in equations {
        chart.draw_series(std::iter::once(Text::new(
            format!("{} {}", label, formula),
            (x, y),
            ("DejaVu Sans", 26).into_font(),
        )))?;
    }
    
    root.present()?;
    println!("Created: {}", filename);
    
    Ok(())
}
```

---

## Task 5: 3D Cosmic Web Visualization

### Step 5.1: Create `src/structure/cosmic_web.rs`

```rust
//! 3D cosmic web generation and visualization

use ndarray::Array3;
use rand::Rng;
use crate::structure::power_spectrum::*;

/// 3D density field
pub struct DensityField {
    pub grid: Array3<f64>,
    pub box_size: f64,  // Mpc/h
    pub n_grid: usize,
}

impl DensityField {
    /// Create new density field
    pub fn new(n_grid: usize, box_size: f64) -> Self {
        DensityField {
            grid: Array3::zeros((n_grid, n_grid, n_grid)),
            box_size,
            n_grid,
        }
    }
    
    /// Generate Gaussian random field
    pub fn generate_gaussian(&mut self, omega_m: f64, z: f64) {
        let mut rng = rand::thread_rng();
        
        // In a full implementation, use FFT with P(k)
        // Here we use simplified correlated noise
        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                for k in 0..self.n_grid {
                    self.grid[[i, j, k]] = rng.gen_range(-1.0..1.0);
                }
            }
        }
        
        // Smooth with Gaussian kernel (simplified)
        self.smooth(2.0);
    }
    
    /// Smooth field with Gaussian kernel
    fn smooth(&mut self, sigma: f64) {
        let kernel_size = (3.0 * sigma) as usize;
        let mut smoothed = self.grid.clone();
        
        for i in kernel_size..(self.n_grid - kernel_size) {
            for j in kernel_size..(self.n_grid - kernel_size) {
                for k in kernel_size..(self.n_grid - kernel_size) {
                    let mut sum = 0.0;
                    let mut weight_sum = 0.0;
                    
                    for di in 0..kernel_size {
                        for dj in 0..kernel_size {
                            for dk in 0..kernel_size {
                                let dx = di as f64 - kernel_size as f64 / 2.0;
                                let dy = dj as f64 - kernel_size as f64 / 2.0;
                                let dz = dk as f64 - kernel_size as f64 / 2.0;
                                let r2 = dx*dx + dy*dy + dz*dz;
                                let weight = (-r2 / (2.0 * sigma * sigma)).exp();
                                
                                sum += self.grid[[i+di-kernel_size/2, j+dj-kernel_size/2, 
                                                 k+dk-kernel_size/2]] * weight;
                                weight_sum += weight;
                            }
                        }
                    }
                    
                    smoothed[[i, j, k]] = sum / weight_sum;
                }
            }
        }
        
        self.grid = smoothed;
    }
    
    /// Extract particles above density threshold
    pub fn extract_particles(&self, threshold: f64) -> Vec<(f64, f64, f64)> {
        let mut particles = Vec::new();
        let dx = self.box_size / self.n_grid as f64;
        
        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                for k in 0..self.n_grid {
                    if self.grid[[i, j, k]] > threshold {
                        let x = i as f64 * dx;
                        let y = j as f64 * dx;
                        let z = k as f64 * dx;
                        particles.push((x, y, z));
                    }
                }
            }
        }
        
        particles
    }
    
    /// Extract slice for 2D visualization
    pub fn slice(&self, axis: usize, index: usize) -> Vec<Vec<f64>> {
        let mut slice_data = vec![vec![0.0; self.n_grid]; self.n_grid];
        
        for i in 0..self.n_grid {
            for j in 0..self.n_grid {
                slice_data[i][j] = match axis {
                    0 => self.grid[[index, i, j]],
                    1 => self.grid[[i, index, j]],
                    2 => self.grid[[i, j, index]],
                    _ => 0.0,
                };
            }
        }
        
        slice_data
    }
}
```

### Step 5.2: Create `examples/structure/cosmic_web_3d.rs`

```rust
//! 3D cosmic web visualization

use andam::structure::cosmic_web::*;
use kiss3d::window::Window;
use kiss3d::light::Light;
use nalgebra::{Point3, Vector3};

fn main() {
    println!("Generating cosmic web...");
    
    let n_grid = 64;
    let box_size = 100.0; // Mpc/h
    
    let mut field = DensityField::new(n_grid, box_size);
    field.generate_gaussian(0.3, 0.0);
    
    // Extract overdense regions
    let threshold = 0.5;
    let particles = field.extract_particles(threshold);
    
    println!("Found {} overdense particles", particles.len());
    
    // Create 3D visualization
    let mut window = Window::new("Cosmic Web Structure");
    window.set_light(Light::StickToCamera);
    window.set_background_color(0.0, 0.0, 0.1); // Dark background
    
    // Add particles
    for (x, y, z) in particles.iter().take(5000) {
        let mut sphere = window.add_sphere(0.3);
        let scale = box_size / 2.0;
        sphere.set_local_translation(nalgebra::Translation3::new(
            (*x - scale) as f32 / scale as f32,
            (*y - scale) as f32 / scale as f32,
            (*z - scale) as f32 / scale as f32,
        ));
        
        // Color by density
        sphere.set_color(0.2, 0.5, 1.0);
    }
    
    // Render loop
    while window.render() {
        // Rotate camera
    }
}
```

### Step 5.3: Create `examples/structure/density_slice.rs`

```rust
//! 2D density slice with high resolution

use andam::structure::cosmic_web::*;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Creating density slice...");
    
    let n_grid = 512;
    let box_size = 100.0;
    
    let mut field = DensityField::new(n_grid, box_size);
    field.generate_gaussian(0.3, 0.0);
    
    // Extract middle slice
    let slice = field.slice(2, n_grid / 2);
    
    let filename = "cosmic_web_slice_hires.png";
    let root = BitMapBackend::new(filename, (2400, 2400)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Cosmic Web Density Field (z=0, 100 Mpc/h slice)", 
                ("DejaVu Sans", 48))
        .margin(10)
        .build_cartesian_2d(0..n_grid, 0..n_grid)?;
    
    // Find min/max for color scaling
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    
    for row in &slice {
        for &val in row {
            min_val = min_val.min(val);
            max_val = max_val.max(val);
        }
    }
    
    // Draw density map
    for i in 0..n_grid {
        for j in 0..n_grid {
            let val = slice[i][j];
            let normalized = (val - min_val) / (max_val - min_val);
            
            // Color map: blue (low) -> red (high)
            let color = if normalized < 0.5 {
                let t = normalized * 2.0;
                RGBColor(
                    0,
                    (255.0 * t) as u8,
                    (255.0 * (1.0 - t)) as u8,
                )
            } else {
                let t = (normalized - 0.5) * 2.0;
                RGBColor(
                    (255.0 * t) as u8,
                    (255.0 * (1.0 - t)) as u8,
                    0,
                )
            };
            
            chart.draw_series(std::iter::once(Rectangle::new(
                [(i, j), (i + 1, j + 1)],
                color.filled(),
            )))?;
        }
    }
    
    // Add scale bar
    let scale_x = n_grid - 100;
    let scale_y = n_grid - 50;
    let scale_len = 50; // pixels = 10 Mpc/h
    
    chart.draw_series(std::iter::once(Rectangle::new(
        [(scale_x, scale_y), (scale_x + scale_len, scale_y + 5)],
        WHITE.filled(),
    )))?;
    
    chart.draw_series(std::iter::once(Text::new(
        "10 Mpc/h",
        (scale_x, scale_y + 15),
        ("DejaVu Sans", 32).into_font().color(&WHITE),
    )))?;
    
    // Add equation
    chart.draw_series(std::iter::once(Text::new(
        "δ(x) = ρ(x)/ρ̄ - 1",
        (20, 40),
        ("DejaVu Sans", 36).into_font().color(&WHITE),
    )))?;
    
    root.present()?;
    println!("Created: {}", filename);
    
    Ok(())
}
```

---

## Phase 6 Completion Checklist

- [ ] HALOFIT non-linear power spectrum
- [ ] Halo mass functions (PS, ST, Tinker)
- [ ] Halo bias implemented
- [ ] Correlation function ξ(r)
- [ ] Redshift-space distortions
- [ ] 3D cosmic web visualization
- [ ] High-resolution 2D slices
- [ ] All tests passing
- [ ] Documentation complete

---

## Expected Outputs

After Phase 6:

1. ✅ `power_spectrum_comparison_hires.png` - Linear vs non-linear
2. ✅ `halo_mass_function_hires.png` - Multiple models with equations
3. ✅ `cosmic_web_3d` - Interactive 3D visualization
4. ✅ `cosmic_web_slice_hires.png` - 2400x2400 density field
5. ✅ Complete structure formation framework

---

## Next Steps

Proceed to **Phase 7: Statistical Methods**
