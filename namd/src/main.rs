use shared::{
    Result,
    //Context,
    //bail,
    info,
    //debug,
};
/*
 *use mpi::{
 *    self,
 *    topology::Communicator,
 *};
 */

//use namd::version::Version;
use namd::logging::logger_init;

fn main() -> Result<()> {
/*
 *    let universe = mpi::initialize()
 *        .with_context(|| "MPI initialization failed!")?;
 *    let world    = universe.world();
 *    let nrank    = world.size();
 *    let irank    = world.rank();
 *    let root_rank = world.process_at_rank(0);
 *
 *
 *    if nrank != 1 {
 *        if 0 == irank {
 *            bail!("MPI parallelism not implemented yet, exiting ...");
 *        } else {
 *            return Ok(());
 *        }
 *    }
 *
 *    if 0 == irank {
 *        println!("MPI initialized with {nrank} ranks.");
 *    }
 */

    let now = std::time::Instant::now();
    logger_init();
    info!("Global logger initialized with targets being stderr and \"./globalrun.log\"");

    {
        use clap::Parser;
        use namd::commands::OptProcess;
        use namd::commands::Command;
        match Command::parse() {
            Command::Run(run) => run.process()?,
        }
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
