use std::path::PathBuf;
use std::fs;

use clap::{Parser, ValueEnum};
use shared::{
    log,
    Result,
};
use crate::OptProcess;
use crate::cli::write_script;
use crate::version::Version;

use crate::core::{
    Hamiltonian,
    NamdConfig,
};
use crate::hamil;


#[derive(Debug, Parser)]
/// Generate the Hamiltonian from NAC according to config file.
#[command(arg_required_else_help(true))]
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
                ConfigTemplate => {
                    log::info!("Writing `02_hamil_config_template.toml` ...");
                    hamil::HamilConfig::default().to_file("02_hamil_config_template.toml")
                },
                EfieldTemplate => {
                    log::info!("Writing `efield_template.rhai` ...");
                    hamil::Efield::template_to_file("efield_template.rhai")
                },
                PostprocessTemplate => {
                    log::info!("Writing `hamil_plot.py` ...");
                    write_script("hamil_plot.py", include_str!("./hamil_plot.py"), true)
                },
            }
        }

        log::info!("\n{}", Version::new());

        let cfg = hamil::HamilConfig::from_file(&self.config)?;
        log::info!("Got Hamiltonnian config:\n{}", &cfg);
        if let Some(e) = cfg.get_efield_fname() {
            let src = fs::read_to_string(e)?;
            log::info!("Got electric field from file {:?} with content of:\n{}", e, src);
        }
        let ham = hamil::SPHamiltonian::from_config(&cfg)?;
        ham.save_to_h5(cfg.get_hamil_fname())
    }
}
