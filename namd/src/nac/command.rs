use std::path::PathBuf;
use std::fs;

use rayon;

use clap::{Parser, ValueEnum};
use shared::{
    log,
    Result,
};
use crate::OptProcess;
use crate::cli::write_script;

use crate::core::{
    Couplings,
    NamdConfig,
};

use crate::nac;

#[derive(Debug, Parser)]
/// Calculate non-adiabatic coupling (NAC) including `<j| d/dt |k>` and
/// momentum matrix `<i| p |j>`.
#[command(arg_required_else_help(true))]
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
                ConfigTemplate => {
                    log::info!("Writing `01_nac_config_template.toml ...`");
                    nac::NacConfig::default().to_file("01_nac_config_template.toml")
                },
                PostprocessTemplate => {
                    log::info!("Writing `nac_plot.py` ...");
                    write_script("nac_plot.py", include_str!("./nac_plot.py"), true)
                },
            }
        }

        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let cfg  = nac::NacConfig::from_file(&self.config)?;
        log::info!("Got NAC config:\n{}", &cfg);
        if cfg.get_nacfname().is_file() {
            log::info!("Found existing NAC: {:?}, exiting.", cfg.get_nacfname());
            return Ok(())
        }

        let coup = nac::Nac::from_config(&cfg)?;
        coup.save_to_h5(cfg.get_nacfname())
    }
}
