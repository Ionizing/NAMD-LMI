use std::path::PathBuf;

use rayon;

use clap::{Parser, ValueEnum};
use shared::{
    log,
    Result,
};
use crate::OptProcess;
use crate::version::Version;
use crate::waveslice;
use crate::core::NamdConfig;


#[derive(Debug, Parser)]
/// Slicing WAVECARs and PROCARs from pre-calculated AIMD trajectory.
#[command(arg_required_else_help(true))]
pub struct WavesliceCommand {
    #[arg(short='n', long, default_value_t=0)]
    /// Number of threads for parallel calculation.
    ///
    /// If 0 is set, it will fall back to the number of logical CPU cores of your machine.
    nthreads: usize,

    #[arg(short='c', long, default_value="waveslice_config.toml", aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: PathBuf,

    #[arg(long)]
    /// Reconstruct waveslice from existing waveslice and new config file.
    ///
    /// Only 'phasecorrection', 'unitary_transform' and 'rearrangement' changes are allowed in
    /// new configuration file.
    reconstruct: Option<PathBuf>,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation.
    ///
    /// The calculation will not run if this flag is set.
    ///
    /// Alias: "gen"
    generate: Option<TemplateGenerator>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, Ord, PartialOrd, ValueEnum)]
enum TemplateGenerator {
    #[value(aliases=["config", "cfg", "conf"])]
    /// Generate config template for Waveslice calculation. Aliases: "config", "cfg", "conf".
    ConfigTemplate,
}


impl OptProcess for WavesliceCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator::*;

        if let Some(g) = self.generate {
            return match g {
                ConfigTemplate => {
                    log::info!("Writing `00_waveslice_config_template.toml` ...");
                    waveslice::WavesliceConfig::default().to_file("00_waveslice_config_template.toml")
                },
            }
        }

        log::info!("\n{}", Version::new());

        log::info!("Running with {} threads.", self.nthreads);
        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let cfg = waveslice::WavesliceConfig::from_file(&self.config)?;
        log::info!("Got WaveSlice config:\n{}", &cfg);
        if cfg.get_waveslicefname().is_file() {
            log::info!("Found existing NAC: {:?}, exiting.", cfg.get_waveslicefname());
            return Ok(());
        }

        let ws = waveslice::Waveslice::from_config(&cfg)?;

        log::info!("Saving to {:?}.", cfg.get_waveslicefname());
        ws.save_to_h5(cfg.get_waveslicefname())
    }
}
