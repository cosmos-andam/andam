//! MCMC chain storage and retrieval

use super::io::{DataStore, StorageError};
use hdf5::Group;
use ndarray::Array2;

/// MCMC storage
pub struct MCMCStorage<'a> {
    group: Group,
    #[allow(dead_code)]
    store: &'a DataStore,
}

impl<'a> MCMCStorage<'a> {
    pub fn new(store: &'a DataStore) -> Result<Self, StorageError> {
        let group = store.mcmc()?;
        Ok(MCMCStorage { group, store })
    }

    /// Store MCMC chain
    pub fn store_chain(
        &self,
        chain_name: &str,
        samples: &Array2<f64>,
        parameter_names: &[String],
        log_probs: &[f64],
    ) -> Result<(), StorageError> {
        let chain_group = self.group.create_group(chain_name)?;

        // Store samples (n_samples x n_params)
        let samples_dataset = chain_group
            .new_dataset::<f64>()
            .shape(samples.dim())
            .chunk([1000.min(samples.nrows()), samples.ncols()]) // Chunk for efficient access
            .deflate(6)
            .create("samples")?;
        samples_dataset.write(samples)?;

        // Store log probabilities
        let logprob_dataset = chain_group
            .new_dataset::<f64>()
            .shape([log_probs.len()])
            .deflate(6)
            .create("log_probability")?;
        logprob_dataset.write(log_probs)?;

        // Store parameter names as JSON dataset
        let names_json = serde_json::to_string(&parameter_names)?;
        let names_dataset = chain_group.new_dataset::<u8>()
            .shape([names_json.len()])
            .create("parameter_names")?;
        names_dataset.write_raw(names_json.as_bytes())?;

        // Store chain info
        chain_group.new_attr::<usize>()
            .create("n_samples")?.write_scalar(&samples.nrows())?;
        chain_group.new_attr::<usize>()
            .create("n_params")?.write_scalar(&samples.ncols())?;

        Ok(())
    }

    /// Read MCMC chain
    pub fn read_chain(&self, chain_name: &str) -> Result<MCMCChain, StorageError> {
        let chain_group = self.group.group(chain_name)?;

        let samples_dataset = chain_group.dataset("samples")?;
        let samples: Array2<f64> = samples_dataset.read()?;

        let log_probs: Vec<f64> = chain_group.dataset("log_probability")?.read_raw()?;

        // Read parameter names from JSON dataset
        let names_bytes: Vec<u8> = chain_group.dataset("parameter_names")?.read_raw()?;
        let names_json = String::from_utf8(names_bytes).map_err(|e| {
            StorageError::IoError(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        })?;
        let parameter_names: Vec<String> = serde_json::from_str(&names_json)?;

        Ok(MCMCChain {
            samples,
            log_probs,
            parameter_names,
        })
    }

    /// Store chain statistics
    pub fn store_statistics(
        &self,
        chain_name: &str,
        stats: &ChainStatistics,
    ) -> Result<(), StorageError> {
        let stats_group = self.group
            .group(chain_name)?
            .create_group("statistics")?;

        // Store means
        stats_group.new_dataset::<f64>()
            .shape([stats.means.len()])
            .create("means")?
            .write(&stats.means)?;

        // Store std devs
        stats_group.new_dataset::<f64>()
            .shape([stats.std_devs.len()])
            .create("std_devs")?
            .write(&stats.std_devs)?;

        // Store covariance matrix
        stats_group.new_dataset::<f64>()
            .shape(stats.covariance.dim())
            .create("covariance")?
            .write(&stats.covariance)?;

        // Convergence diagnostics
        stats_group.new_attr::<f64>()
            .create("gelman_rubin")?.write_scalar(&stats.gelman_rubin)?;
        stats_group.new_attr::<f64>()
            .create("acceptance_rate")?.write_scalar(&stats.acceptance_rate)?;

        Ok(())
    }

    /// Read chain statistics
    pub fn read_statistics(&self, chain_name: &str) -> Result<ChainStatistics, StorageError> {
        let stats_group = self.group.group(chain_name)?.group("statistics")?;

        let means: Vec<f64> = stats_group.dataset("means")?.read_raw()?;
        let std_devs: Vec<f64> = stats_group.dataset("std_devs")?.read_raw()?;
        let covariance: Array2<f64> = stats_group.dataset("covariance")?.read()?;
        let gelman_rubin: f64 = stats_group.attr("gelman_rubin")?.read_scalar()?;
        let acceptance_rate: f64 = stats_group.attr("acceptance_rate")?.read_scalar()?;

        Ok(ChainStatistics {
            means,
            std_devs,
            covariance,
            gelman_rubin,
            acceptance_rate,
        })
    }

    /// Store thinned chain (every nth sample)
    pub fn store_thinned_chain(
        &self,
        chain_name: &str,
        thin_factor: usize,
    ) -> Result<(), StorageError> {
        let chain = self.read_chain(chain_name)?;

        let n_samples = chain.samples.nrows();
        let n_thinned = n_samples / thin_factor;

        let mut thinned_samples = Array2::zeros((n_thinned, chain.samples.ncols()));
        let mut thinned_logprobs = Vec::with_capacity(n_thinned);

        for i in 0..n_thinned {
            let idx = i * thin_factor;
            thinned_samples.row_mut(i).assign(&chain.samples.row(idx));
            thinned_logprobs.push(chain.log_probs[idx]);
        }

        let thinned_name = format!("{}_thinned_{}", chain_name, thin_factor);
        self.store_chain(
            &thinned_name,
            &thinned_samples,
            &chain.parameter_names,
            &thinned_logprobs,
        )?;

        Ok(())
    }
}

/// MCMC chain data
#[derive(Debug, Clone)]
pub struct MCMCChain {
    pub samples: Array2<f64>,
    pub log_probs: Vec<f64>,
    pub parameter_names: Vec<String>,
}

/// Chain statistics
#[derive(Debug, Clone)]
pub struct ChainStatistics {
    pub means: Vec<f64>,
    pub std_devs: Vec<f64>,
    pub covariance: Array2<f64>,
    pub gelman_rubin: f64,
    pub acceptance_rate: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_store_and_read_chain() {
        let temp_file = NamedTempFile::new().unwrap();
        let store = DataStore::create(temp_file.path()).unwrap();
        let mcmc = MCMCStorage::new(&store).unwrap();

        let samples = Array2::from_shape_fn((100, 2), |(i, j)| (i + j) as f64);
        let log_probs: Vec<f64> = (0..100).map(|i| -0.5 * (i as f64)).collect();
        let param_names = vec!["param1".to_string(), "param2".to_string()];

        mcmc.store_chain("test_chain", &samples, &param_names, &log_probs).unwrap();

        let read_chain = mcmc.read_chain("test_chain").unwrap();
        assert_eq!(read_chain.samples.nrows(), 100);
        assert_eq!(read_chain.samples.ncols(), 2);
        assert_eq!(read_chain.parameter_names.len(), 2);
    }
}
