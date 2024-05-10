use clap::Args;
use crate::commands::OptProcess;
use crate::commands::Result;


#[derive(Args)]
pub struct Generate {
    #[arg(default_value_t = String::from("input.toml") )]
    /// Where to write the example input file.
    filename: String,
}


impl OptProcess for Generate {
    fn process(&self) -> Result<()> {
        todo!()
    }
}
