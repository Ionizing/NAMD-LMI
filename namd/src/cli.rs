use std::sync::OnceLock;

use clap::{
    Parser,
    builder::styling::{
        AnsiColor,
        Effects,
        Styles,
    },
};

use shared::{log, Result};
use crate::version::Version;
use crate::logging::logger_init;


pub fn get_style() -> Styles {
    static INSTANCE: OnceLock<Styles> = OnceLock::new();
    INSTANCE.get_or_init(|| {
        Styles::styled()
            .header(AnsiColor::Yellow.on_default() | Effects::BOLD)
            .usage(AnsiColor::Green.on_default()   | Effects::BOLD)
            .literal(AnsiColor::Green.on_default() | Effects::BOLD)
            .placeholder(AnsiColor::BrightBlue.on_default())
            .error(AnsiColor::BrightRed.on_default())
            .valid(AnsiColor::BrightYellow.on_default())
    }).to_owned()
}


pub trait OptProcess : Parser {
    fn process(&self) -> Result<()>;
}


#[derive(Debug, Parser)]
#[command(name = "NAMD-LUMI",
          about = Version::new().to_string(),
          long_about = format!("{:#}", Version::new()),
          version,
          author = "@Ionizing github.com/Ionizing/NAMD-lumi",
          styles = get_style())]
enum Opt {
    Nac(crate::nac::NacCommand),
    Hamil(crate::hamil::HamilCommand),
    Surfhop(crate::surfhop::SurfhopCommand),
}


impl OptProcess for Opt {
    fn process(&self) -> Result<()> {
        use Opt::*;
        
        logger_init();
        log::info!("Global logger initialized with targets being stderr and \"./globalrun.log\"");

        match self {
            Nac(cmd) => cmd.process(),
            Hamil(cmd) => cmd.process(),
            Surfhop(cmd) => cmd.process(),
        }
    }
}


pub fn run() -> Result<()> {
    Opt::parse().process()
}
