//! Core HDF5 I/O utilities

use hdf5::{File, Group, Result as H5Result};
use std::path::Path;
use thiserror::Error;

/// Storage errors
#[derive(Error, Debug)]
pub enum StorageError {
    #[error("HDF5 error: {0}")]
    HDF5Error(#[from] hdf5::Error),

    #[error("Dataset not found: {0}")]
    DatasetNotFound(String),

    #[error("Invalid dimensions")]
    InvalidDimensions,

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("JSON error: {0}")]
    JsonError(#[from] serde_json::Error),
}

/// Main data store
pub struct DataStore {
    pub(crate) file: File,
}

impl DataStore {
    /// Create new HDF5 file
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self, StorageError> {
        let file = File::create(path)?;

        // Create standard group structure
        file.create_group("cmb")?;
        file.create_group("mcmc")?;
        file.create_group("structure")?;
        file.create_group("parameters")?;
        file.create_group("metadata")?;

        Ok(DataStore { file })
    }

    /// Open existing HDF5 file
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, StorageError> {
        let file = File::open(path)?;
        Ok(DataStore { file })
    }

    /// Open with read-write access
    pub fn open_rw<P: AsRef<Path>>(path: P) -> Result<Self, StorageError> {
        let file = File::open_rw(path)?;
        Ok(DataStore { file })
    }

    /// Get CMB group
    pub fn cmb(&self) -> H5Result<Group> {
        self.file.group("cmb")
    }

    /// Get MCMC group
    pub fn mcmc(&self) -> H5Result<Group> {
        self.file.group("mcmc")
    }

    /// Get structure group
    pub fn structure(&self) -> H5Result<Group> {
        self.file.group("structure")
    }

    /// Get parameters group
    pub fn parameters(&self) -> H5Result<Group> {
        self.file.group("parameters")
    }

    /// Store metadata as JSON
    pub fn store_metadata(&self, key: &str, value: &serde_json::Value) -> Result<(), StorageError> {
        let metadata = self.file.group("metadata")?;
        let json_str = serde_json::to_string_pretty(value)?;

        // Store as dataset instead of attribute for better compatibility
        let dataset = metadata
            .new_dataset::<u8>()
            .shape([json_str.len()])
            .create(key)?;
        dataset.write_raw(json_str.as_bytes())?;

        Ok(())
    }

    /// Read metadata
    pub fn read_metadata(&self, key: &str) -> Result<serde_json::Value, StorageError> {
        let metadata = self.file.group("metadata")?;
        let dataset = metadata.dataset(key)?;
        let bytes: Vec<u8> = dataset.read_raw()?;
        let json_str = String::from_utf8(bytes).map_err(|e| {
            StorageError::IoError(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        })?;
        let value = serde_json::from_str(&json_str)?;
        Ok(value)
    }

    /// List all datasets in a group
    pub fn list_datasets(&self, group_name: &str) -> Result<Vec<String>, StorageError> {
        let group = self.file.group(group_name)?;
        let mut datasets = Vec::new();

        for name in group.member_names()? {
            if group.dataset(&name).is_ok() {
                datasets.push(name);
            }
        }

        Ok(datasets)
    }

    /// Get file size in bytes
    pub fn file_size(&self) -> Result<u64, StorageError> {
        let path = self.file.filename();
        let metadata = std::fs::metadata(&path)?;
        Ok(metadata.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_create_datastore() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();

        assert!(store.cmb().is_ok());
        assert!(store.mcmc().is_ok());
    }

    #[test]
    fn test_metadata_storage() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();

        let meta = serde_json::json!({
            "omega_m": 0.3,
            "sigma_8": 0.8,
        });

        store.store_metadata("cosmology", &meta).unwrap();
        let read_meta = store.read_metadata("cosmology").unwrap();

        assert_eq!(meta, read_meta);
    }
}
