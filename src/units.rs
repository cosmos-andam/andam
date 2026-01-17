//! Unit conversion utilities for cosmology
//!
//! Provides conversions between common astronomical and cosmological units

use crate::constants::*;

/// Length units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Length {
    Meter(f64),
    Kilometer(f64),
    AU(f64),
    Parsec(f64),
    Kiloparsec(f64),
    Megaparsec(f64),
    Gigaparsec(f64),
}

impl Length {
    /// Convert to meters
    pub fn to_meters(&self) -> f64 {
        match self {
            Length::Meter(x) => *x,
            Length::Kilometer(x) => x * 1e3,
            Length::AU(x) => x * AU,
            Length::Parsec(x) => x * PARSEC,
            Length::Kiloparsec(x) => x * PARSEC * 1e3,
            Length::Megaparsec(x) => x * PARSEC * 1e6,
            Length::Gigaparsec(x) => x * PARSEC * 1e9,
        }
    }

    /// Convert to megaparsecs
    pub fn to_mpc(&self) -> f64 {
        self.to_meters() / (PARSEC * 1e6)
    }

    /// Convert to parsecs
    pub fn to_pc(&self) -> f64 {
        self.to_meters() / PARSEC
    }
}

/// Mass units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mass {
    Kilogram(f64),
    SolarMass(f64),
    EarthMass(f64),
}

impl Mass {
    /// Convert to kilograms
    pub fn to_kg(&self) -> f64 {
        match self {
            Mass::Kilogram(x) => *x,
            Mass::SolarMass(x) => x * M_SUN,
            Mass::EarthMass(x) => x * 5.972e24,
        }
    }

    /// Convert to solar masses
    pub fn to_solar_masses(&self) -> f64 {
        self.to_kg() / M_SUN
    }
}

/// Time units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Time {
    Second(f64),
    Year(f64),
    Kiloyear(f64),
    Megayear(f64),
    Gigayear(f64),
}

impl Time {
    /// Convert to seconds
    pub fn to_seconds(&self) -> f64 {
        match self {
            Time::Second(x) => *x,
            Time::Year(x) => x * YEAR,
            Time::Kiloyear(x) => x * YEAR * 1e3,
            Time::Megayear(x) => x * YEAR * 1e6,
            Time::Gigayear(x) => x * YEAR * 1e9,
        }
    }

    /// Convert to gigayears
    pub fn to_gyr(&self) -> f64 {
        self.to_seconds() / GYR
    }
}

/// Energy units
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Energy {
    Joule(f64),
    ElectronVolt(f64),
    KiloElectronVolt(f64),
    MegaElectronVolt(f64),
    GeV(f64),
}

impl Energy {
    /// Convert to joules
    pub fn to_joules(&self) -> f64 {
        match self {
            Energy::Joule(x) => *x,
            Energy::ElectronVolt(x) => x * EV_TO_J,
            Energy::KiloElectronVolt(x) => x * EV_TO_J * 1e3,
            Energy::MegaElectronVolt(x) => x * EV_TO_J * 1e6,
            Energy::GeV(x) => x * EV_TO_J * 1e9,
        }
    }

    /// Convert to eV
    pub fn to_ev(&self) -> f64 {
        self.to_joules() * J_TO_EV
    }
}

/// Temperature conversion utilities
pub mod temperature {
    use super::*;

    /// Convert temperature to energy (E = kT)
    pub fn temp_to_energy_ev(temp_kelvin: f64) -> f64 {
        (K_B * temp_kelvin) * J_TO_EV
    }

    /// Convert energy to temperature (T = E/k)
    pub fn energy_to_temp(energy_ev: f64) -> f64 {
        (energy_ev * EV_TO_J) / K_B
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_length_conversions() {
        let d = Length::Megaparsec(1.0);
        assert_relative_eq!(d.to_mpc(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(d.to_pc(), 1e6, epsilon = 1e-4);
    }

    #[test]
    fn test_time_conversions() {
        let t = Time::Gigayear(1.0);
        assert_relative_eq!(t.to_gyr(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_temperature_energy() {
        let temp = 1e4; // 10,000 K
        let energy = temperature::temp_to_energy_ev(temp);
        let back = temperature::energy_to_temp(energy);
        assert_relative_eq!(temp, back, epsilon = 1e-6);
    }
}
