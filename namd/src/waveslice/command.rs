use std::path::PathBuf;

use rayon;

use clap::{Parser, ValueEnum, Subcommand};
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

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation.
    ///
    /// The calculation will not run if this flag is set.
    ///
    /// Alias: "gen"
    generate: Option<TemplateGenerator>,


    #[command(subcommand)]
    /// Perform transforms for given waveslice.
    transform: Option<Transform>,
}


#[derive(Subcommand, Debug, Clone)]
#[command(arg_required_else_help(true))]
enum Transform {
    Transform {
        #[arg(short='i', long, default_value="waveslice.h5")]
        /// Input waveslice filename.
        ifname: PathBuf,

        #[arg(short='o', long)]
        /// Output waveslice filename.
        ofname: Option<PathBuf>,

        #[arg(long, required(true), num_args(1..))]
        /// Operations to be applied.
        ///
        /// You can perform multiple operations in way.
        operations: Vec<Transformations>,

        #[arg(long, default_value_t=1.0)]
        /// The gausian width applied during the rearrangement.
        ///
        /// This factor works when calculating
        /// c_ij = |<phi_j(t') | phi_i(t)>|^2 * exp(-|Ej(t') - Ei(t)|^2)
        sigma: f64,

        #[arg(long)]
        /// Use the order file to rearrange the bands.
        ///
        /// This file can be obtained by performing rearrangement once.
        order_file: Option<PathBuf>,
    }
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, Ord, PartialOrd, ValueEnum)]
enum Transformations {
    #[value(aliases=["rearr", "reordering", "reorder"])]
    /// Re-arrange the band order to solve the band crossing problem.
    Rearrangement,

    #[value(aliases=["phase", "phasecorr"])]
    /// Apply phase correction of each band.
    PhaseCorrection,

    //UnitaryTransform,
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

        if let Some(transform) = self.transform.as_ref() {
            match transform {
                Transform::Transform { ifname, ofname, operations, order_file, sigma } => {
                    let rearrangement = operations.contains(&Transformations::Rearrangement);
                    let phase_correction = operations.contains(&Transformations::PhaseCorrection);

                    let mut suffix = String::new();

                    if rearrangement {
                        log::info!("Found rearrangement operation");
                        suffix.push_str("_rearr");
                    }

                    if phase_correction {
                        log::info!("Found rearrangement operation");
                        suffix.push_str("_phcor");
                    }

                    let ofname = if ofname.is_none() {
                        let mut fname = ifname.with_extension("").to_str().unwrap().to_string();
                        fname.push_str(&suffix);
                        PathBuf::from(fname).with_extension("h5")
                    } else {
                        ofname.as_ref().unwrap().to_owned()
                    };

                    log::info!("Output file {:?} will be written.", &ofname);
                    log::info!("Loading {:?} ...", ifname);
                    let mut ws = waveslice::Waveslice::from_h5(ifname)?;
                    
                    if rearrangement {
                        if let Some(order_fname) = order_file.as_ref() {
                            ws.apply_rearrangement_with_given_order(order_fname)?;
                        } else {
                            log::info!("Performing rearrangement ...");
                            ws.apply_rearrangement(*sigma);
                        }
                    }

                    if phase_correction {
                        log::info!("Performing phase correction ...");
                        ws.apply_phase_correction();
                    }

                    log::info!("Saving to {:?} ...", ofname);
                    ws.save_to_h5(ofname)?;
                }
            }

            return Ok(())
        }

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
