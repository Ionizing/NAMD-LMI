use clap::Parser;
use shared::Result;
use crate::OptProcess;

#[derive(Debug, Parser)]
pub struct NacCommand {

}

impl OptProcess for NacCommand {
    fn process(&self) -> Result<()> {
        todo!()
    }
}
