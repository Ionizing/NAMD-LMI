use std::sync::OnceLock;

use clap::{
    Parser,
    builder::styling::{
        AnsiColor,
        Effects,
        Styles,
    },
};

use shared::Result;

use crate::nac::command::NacCommand;


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
          about = "Non-Adiabatic Molecular Dynamics with external fields.",
          version,
          author = "@Ionizing github.com/Ionizing/NAMD-lumi",
          styles = get_style())]
enum Opt {
    Nac(NacCommand),
    Hamil,
    Surfhop,
}


impl OptProcess for Opt {
    fn process(&self) -> Result<()> {
        use Opt::*;
        
        match self {
            Nac(cmd) => cmd.process(),
            _ => todo!(),
        }
    }
}


pub fn run() -> Result<()> {
    Opt::parse().process()
}
