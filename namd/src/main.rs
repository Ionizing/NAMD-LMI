use shared::{
    Result,
    //Context,
    //bail,
    info,
    //debug,
};
use clap::Parser;
//use namd::commands::OptProcess;
//use namd::commands::Command;

//use namd::logging::logger_init;

fn main() -> Result<()> {
    let now = std::time::Instant::now();
    //let cml = Command::parse();

    //logger_init();
    info!("Global logger initialized with targets being stderr and \"./globalrun.log\"");

    //match cml {
        //Command::Run(run) => run.process()?,
    //}

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
