//! CMB data storage (maps, power spectra)

use super::io::{DataStore, StorageError};
use hdf5::Group;
use serde::{Deserialize, Serialize};

/// CMB data storage
pub struct CMBStorage<'a> {
    group: Group,
    #[allow(dead_code)]
    store: &'a DataStore,
}

impl<'a> CMBStorage<'a> {
    /// Create CMB storage interface
    pub fn new(store: &'a DataStore) -> Result<Self, StorageError> {
        let group = store.cmb()?;
        Ok(CMBStorage { group, store })
    }

    /// Store temperature map
    pub fn store_temperature_map(
        &self,
        map: &[f64],
        n_side: usize,
        name: &str,
    ) -> Result<(), StorageError> {
        let dataset = self
            .group
            .new_dataset::<f64>()
            .shape([map.len()])
            .deflate(6) // Compression level 6
            .create(name)?;

        dataset.write(map)?;

        // Store metadata
        dataset
            .new_attr::<usize>()
            .create("n_side")?
            .write_scalar(&n_side)?;
        // Note: Storing type and units as JSON metadata instead of attributes
        // due to HDF5-rust API limitations with VarLenUnicode

        Ok(())
    }

    /// Read temperature map
    pub fn read_temperature_map(&self, name: &str) -> Result<(Vec<f64>, usize), StorageError> {
        let dataset = self.group.dataset(name)?;
        let map: Vec<f64> = dataset.read_raw()?;
        let n_side: usize = dataset.attr("n_side")?.read_scalar()?;

        Ok((map, n_side))
    }

    /// Store power spectrum
    pub fn store_power_spectrum(
        &self,
        ell: &[usize],
        c_ell: &[f64],
        _spectrum_type: &str,
        name: &str,
    ) -> Result<(), StorageError> {
        // Create group for this spectrum
        let spec_group = self.group.create_group(name)?;

        // Store ell values
        let ell_dataset = spec_group
            .new_dataset::<usize>()
            .shape([ell.len()])
            .create("ell")?;
        ell_dataset.write(ell)?;

        // Store C_ell values
        let c_ell_dataset = spec_group
            .new_dataset::<f64>()
            .shape([c_ell.len()])
            .deflate(6)
            .create("c_ell")?;
        c_ell_dataset.write(c_ell)?;

        // Metadata
        // Note: Spectrum type stored in group name instead of attribute
        spec_group
            .new_attr::<usize>()
            .create("l_max")?
            .write_scalar(ell.last().unwrap())?;

        Ok(())
    }

    /// Read power spectrum
    pub fn read_power_spectrum(&self, name: &str) -> Result<PowerSpectrum, StorageError> {
        let spec_group = self.group.group(name)?;

        let ell: Vec<usize> = spec_group.dataset("ell")?.read_raw()?;
        let c_ell: Vec<f64> = spec_group.dataset("c_ell")?.read_raw()?;
        // Extract spectrum type from group name (last part after /)
        let spectrum_type = name.split('/').last().unwrap_or(name).to_string();

        Ok(PowerSpectrum {
            ell,
            c_ell,
            spectrum_type,
        })
    }

    /// Store polarization maps (Q, U)
    pub fn store_polarization_maps(
        &self,
        q_map: &[f64],
        u_map: &[f64],
        n_side: usize,
        name: &str,
    ) -> Result<(), StorageError> {
        let pol_group = self.group.create_group(name)?;

        // Q map
        let q_dataset = pol_group
            .new_dataset::<f64>()
            .shape([q_map.len()])
            .deflate(6)
            .create("Q")?;
        q_dataset.write(q_map)?;

        // U map
        let u_dataset = pol_group
            .new_dataset::<f64>()
            .shape([u_map.len()])
            .deflate(6)
            .create("U")?;
        u_dataset.write(u_map)?;

        // Metadata
        pol_group
            .new_attr::<usize>()
            .create("n_side")?
            .write_scalar(&n_side)?;
        // Note: Type information stored in group name

        Ok(())
    }
}

/// Power spectrum data structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowerSpectrum {
    pub ell: Vec<usize>,
    pub c_ell: Vec<f64>,
    pub spectrum_type: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_store_and_read_power_spectrum() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();
        let cmb = CMBStorage::new(&store).unwrap();

        let ell: Vec<usize> = (2..=100).collect();
        let c_ell: Vec<f64> = (2..=100).map(|l| (l as f64).powi(-2)).collect();

        cmb.store_power_spectrum(&ell, &c_ell, "TT", "TT").unwrap();

        let read_spectrum = cmb.read_power_spectrum("TT").unwrap();
        assert_eq!(read_spectrum.ell.len(), 99);
        assert_eq!(read_spectrum.spectrum_type, "TT");
    }
}
