use clap::Parser;

use shared::Result;
use crate::version::Version;

pub mod generate;
pub mod run;
pub mod scripts;
pub mod tutorial;


pub trait OptProcess {
    fn process(&self) -> Result<()>;
}


#[derive(Parser)]
#[command(author,
          version,
          propagate_version = true,
          about = Version::new().to_string(),
          long_about= format!("{:#}", Version::new()) )]
pub enum Command {
    Run(run::Run),
}
