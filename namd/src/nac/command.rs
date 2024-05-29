use std::path::PathBuf;

use rayon;

use clap::{Parser, ValueEnum};
use shared::{
    log,
    Result,
};
use crate::OptProcess;

use crate::core::{
    Couplings,
    NamdConfig,
};

use crate::nac;

#[derive(Debug, Parser)]
/// Calculate non-adiabatic coupling (NAC) including `<j| d/dt |k>` and
/// momentum matrix `<i| p |j>`.
pub struct NacCommand {
    #[arg(short='n', long, default_value_t=0)]
    /// Number of threads for parallel calculation.
    /// 
    /// If 0 is set, it will fall back to the number of logic CPU cores of you machine.
    nthreads: usize,

    #[arg(short='c', long, default_value="nac_config.toml", aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: PathBuf,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation and analysis.
    ///
    /// The calculation will not run if this flag is set.
    ///
    /// Alias: "gen"
    generate: Option<TemplateGenerator>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum TemplateGenerator {
    #[value(aliases=["config", "cfg", "conf"])]
    /// Generate config template for NAC calculation. Aliases: "config", "cfg", "conf".
    ConfigTemplate,

    #[value(aliases=["post-process", "pp", "postprocess"])]
    /// Generate post-process scripts for NAC analysis. Aliases: "post-process", "postprocess", "pp".
    PostprocessTemplate,
}


impl OptProcess for NacCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator::*;

        if let Some(g) = self.generate {
            return match g {
                ConfigTemplate => nac::NacConfig::default().to_file("nac_config_template.toml"),
                PostprocessTemplate => todo!(),
            }
        }

        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let cfg  = nac::NacConfig::from_file(&self.config)?;
        if cfg.get_nacfname().is_file() {
            log::info!("Found existing NAC: {:?}, exiting.", cfg.get_nacfname());
            return Ok(())
        }

        let coup = nac::Nac::from_config(&cfg)?;
        coup.save_to_h5(cfg.get_nacfname())
    }
}
