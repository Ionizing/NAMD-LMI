use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use shared::Result;
use crate::OptProcess;

use crate::core::{
    Hamiltonian,
    NamdConfig,
};
use crate::hamil;


#[derive(Debug, Parser)]
/// Generate the Hamiltonian from NAC according to config file.
pub struct HamilCommand {
    #[arg(short='c', long, default_value="hamil_config.toml", aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: PathBuf,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation and analysis.
    ///
    /// The generation of Hamiltonian will not run if this flag is set.
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

    #[value(aliases=["efield", "ef"])]
    /// Generate script template for external electric field.
    /// Aliases: "efield", "ef".
    EfieldTemplate,

    #[value(aliases=["post-process", "postprocess", "pp"])]
    /// Generate post-process scripts for Hamiltonian analysis.
    /// Aliases: "post-process", "postprocess", "pp".
    PostprocessTemplate,
}


impl OptProcess for HamilCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator::*;

        if let Some(g) = self.generate {
            return match g {
                ConfigTemplate => hamil::HamilConfig::default().to_file("hamil_config_template.toml"),
                EfieldTemplate => hamil::Efield::template_to_file("efield_template.rhai"),
                PostprocessTemplate => todo!(),
            }
        }

        let cfg = hamil::HamilConfig::from_file(&self.config)?;
        let ham = hamil::SPHamiltonian::from_config(&cfg)?;
        ham.save_to_h5(cfg.get_hamil_fname())
    }
}
