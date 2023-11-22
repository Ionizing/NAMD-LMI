use clap::Parser;

pub use shared::Result;

pub mod generate;
pub mod run;
pub mod scripts;
pub mod tutorial;


pub trait OptProcess {
    fn process(&self) -> Result<()>;
}


#[derive(Parser)]
#[command(author, version, about, long_about=None)]
#[command(propagate_version = true)]
pub enum Command {
    //#[command(subcommand)]
    Run(run::Run),
}
