use shared::{
    Result,
    //Context,
    //bail,
    info,
    //debug,
};

use namd::logging::logger_init;
use namd::cli::run;

fn main() -> Result<()> {
    let now = std::time::Instant::now();

    logger_init();
    info!("Global logger initialized with targets being stderr and \"./globalrun.log\"");

    run()?;

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
