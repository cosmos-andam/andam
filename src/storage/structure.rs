//! Structure formation data storage

use super::io::{DataStore, StorageError};
use hdf5::Group;
use ndarray::{Array2, Array3};

/// Structure storage
pub struct StructureStorage<'a> {
    group: Group,
    #[allow(dead_code)]
    store: &'a DataStore,
}

impl<'a> StructureStorage<'a> {
    pub fn new(store: &'a DataStore) -> Result<Self, StorageError> {
        let group = store.structure()?;
        Ok(StructureStorage { group, store })
    }

    /// Store matter power spectrum
    pub fn store_power_spectrum(
        &self,
        k: &[f64],
        p_k: &[f64],
        z: f64,
        linear: bool,
        name: &str,
    ) -> Result<(), StorageError> {
        let ps_group = self.group.create_group(name)?;

        // Store k values
        ps_group.new_dataset::<f64>()
            .shape([k.len()])
            .create("k")?.write(k)?;

        // Store P(k) values
        ps_group.new_dataset::<f64>()
            .shape([p_k.len()])
            .deflate(6)
            .create("power_spectrum")?.write(p_k)?;

        // Metadata
        ps_group.new_attr::<f64>().create("redshift")?.write_scalar(&z)?;
        ps_group.new_attr::<bool>().create("linear")?.write_scalar(&linear)?;
        // Note: Units information can be inferred from context

        Ok(())
    }

    /// Store 3D density field
    pub fn store_density_field(
        &self,
        field: &Array3<f64>,
        box_size: f64,
        z: f64,
        name: &str,
    ) -> Result<(), StorageError> {
        let field_dataset = self.group
            .new_dataset::<f64>()
            .shape(field.dim())
            .chunk([64.min(field.shape()[0]), 64.min(field.shape()[1]), 64.min(field.shape()[2])]) // Chunk for efficient slicing
            .deflate(6)
            .create(name)?;

        field_dataset.write(field)?;

        // Metadata
        field_dataset.new_attr::<f64>().create("box_size")?.write_scalar(&box_size)?;
        field_dataset.new_attr::<f64>().create("redshift")?.write_scalar(&z)?;
        field_dataset.new_attr::<usize>().create("n_grid")?.write_scalar(&field.shape()[0])?;

        Ok(())
    }

    /// Read density field (full 3D field)
    pub fn read_density_field(&self, name: &str) -> Result<Array3<f64>, StorageError> {
        let dataset = self.group.dataset(name)?;
        let field: Array3<f64> = dataset.read()?;
        Ok(field)
    }

    /// Read density field slice
    pub fn read_density_slice(
        &self,
        name: &str,
        axis: usize,
        index: usize,
    ) -> Result<Array2<f64>, StorageError> {
        let dataset = self.group.dataset(name)?;
        let shape = dataset.shape();

        if axis >= 3 {
            return Err(StorageError::InvalidDimensions);
        }

        if index >= shape[axis] {
            return Err(StorageError::InvalidDimensions);
        }

        // Read full field then extract slice (HDF5 partial reads are complex)
        let field: Array3<f64> = dataset.read()?;

        let slice = match axis {
            0 => field.index_axis(ndarray::Axis(0), index).to_owned(),
            1 => field.index_axis(ndarray::Axis(1), index).to_owned(),
            2 => field.index_axis(ndarray::Axis(2), index).to_owned(),
            _ => unreachable!(),
        };

        Ok(slice)
    }

    /// Store halo catalog
    pub fn store_halo_catalog(
        &self,
        halos: &HaloCatalog,
        name: &str,
    ) -> Result<(), StorageError> {
        let halo_group = self.group.create_group(name)?;

        // Convert positions to flat array
        let pos_flat: Vec<f64> = halos.positions.iter()
            .flat_map(|p| p.iter().copied())
            .collect();

        // Positions
        halo_group.new_dataset::<f64>()
            .shape([halos.positions.len(), 3])
            .create("positions")?.write_raw(&pos_flat)?;

        // Masses
        halo_group.new_dataset::<f64>()
            .shape([halos.masses.len()])
            .create("masses")?.write(&halos.masses)?;

        // Convert velocities to flat array
        let vel_flat: Vec<f64> = halos.velocities.iter()
            .flat_map(|v| v.iter().copied())
            .collect();

        // Velocities
        halo_group.new_dataset::<f64>()
            .shape([halos.velocities.len(), 3])
            .create("velocities")?.write_raw(&vel_flat)?;

        // Metadata
        halo_group.new_attr::<usize>().create("n_halos")?.write_scalar(&halos.masses.len())?;
        halo_group.new_attr::<f64>().create("redshift")?.write_scalar(&halos.redshift)?;

        Ok(())
    }

    /// Read halo catalog
    pub fn read_halo_catalog(&self, name: &str) -> Result<HaloCatalog, StorageError> {
        let halo_group = self.group.group(name)?;

        let pos_flat: Vec<f64> = halo_group.dataset("positions")?.read_raw()?;
        let positions: Vec<[f64; 3]> = pos_flat.chunks_exact(3)
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect();

        let masses: Vec<f64> = halo_group.dataset("masses")?.read_raw()?;

        let vel_flat: Vec<f64> = halo_group.dataset("velocities")?.read_raw()?;
        let velocities: Vec<[f64; 3]> = vel_flat.chunks_exact(3)
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect();

        let redshift: f64 = halo_group.attr("redshift")?.read_scalar()?;

        Ok(HaloCatalog {
            positions,
            masses,
            velocities,
            redshift,
        })
    }
}

/// Halo catalog
#[derive(Debug, Clone)]
pub struct HaloCatalog {
    pub positions: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
    pub velocities: Vec<[f64; 3]>,
    pub redshift: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_store_power_spectrum() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();
        let structure = StructureStorage::new(&store).unwrap();

        let k: Vec<f64> = (0..100).map(|i| 0.01 * (i as f64 + 1.0)).collect();
        let p_k: Vec<f64> = k.iter().map(|k_val| 100.0 / k_val.powi(2)).collect();

        structure.store_power_spectrum(&k, &p_k, 0.0, true, "test_ps").unwrap();

        // Verify by reading the group
        let ps_group = structure.group.group("test_ps").unwrap();
        let read_k: Vec<f64> = ps_group.dataset("k").unwrap().read_raw().unwrap();
        assert_eq!(read_k.len(), 100);
    }

    #[test]
    fn test_store_density_field() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();
        let structure = StructureStorage::new(&store).unwrap();

        let field = Array3::from_shape_fn((32, 32, 32), |(i, j, k)| (i + j + k) as f64);

        structure.store_density_field(&field, 100.0, 0.0, "test_field").unwrap();

        let read_field = structure.read_density_field("test_field").unwrap();
        assert_eq!(read_field.shape(), &[32, 32, 32]);
    }
}
