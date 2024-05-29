use std::path::PathBuf;

use rayon::prelude::*;

use clap::{Parser, ValueEnum};
use shared::Result;
use crate::OptProcess;

use crate::core::{
    SurfaceHopping,
    NamdConfig,
};
use crate::surfhop;


#[derive(Debug, Parser)]
/// Perform the surface-hopping process with given Hamiltonian file and config file.
pub struct SurfhopCommand {
    #[arg(short='n', long, default_value_t=0)]
    /// Number of threads for parallel calculation.
    /// 
    /// If 0 is set, it will fall back to the number of logic CPU cores of you machine.
    nthreads: usize,

    #[arg(short='c', long, default_value="surfhop_config.toml", aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: PathBuf,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation and analysis.
    ///
    /// The surface-hopping will not run if this flag is set.
    ///
    /// Alias: "gen".
    generate: Option<TemplateGenerator>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum TemplateGenerator {
    #[value(aliases=["config", "cfg", "conf"])]
    /// Generate config template for Hamiltonian generation.
    /// Aliases: "config", "cfg" and "conf".
    ConfigTemplate,

    #[value(aliases=["post-process", "postprocess", "pp"])]
    /// Generate post-process scripts for surface-hopping analysis.
    /// Aliases: "post-process", "postprocess", "pp".
    PostprocessTemplate,
}


impl OptProcess for SurfhopCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator::*;

        if let Some(g) = self.generate {
            return match g {
                ConfigTemplate => surfhop::SurfhopConfig::default().to_file("surfhop_config_template.toml"),
                PostprocessTemplate => todo!(),
            }
        }

        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let cfg = surfhop::SurfhopConfig::from_file(&self.config)?;
        crate::logging::logger_redirect(&cfg.get_outdir())?;
        let sh = surfhop::Surfhop::from_config(&cfg)?;

        sh.into_par_iter()
            .map(|mut v| v.run())
            .collect::<Result<Vec<()>>>()?;

        Ok(())
    }
}
