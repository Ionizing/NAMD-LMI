use shared::{
    //ndarray,
    Result,
    Context,
    bail,
    info,
};
use mpi::{
    self,
    topology::Communicator,
};
use env_logger::init_from_env;

use namd::{
    version::Version,
    input::Input,
    nac::Nac,
    hamiltonian::Hamiltonian,
    surface_hopping::SurfaceHopping,
};

fn main() -> Result<()> {
    let universe = mpi::initialize()
        .with_context(|| "MPI initialization failed!")?;
    let world    = universe.world();
    let nrank     = world.size();
    let irank    = world.rank();
    //let root_rank = world.process_at_rank(0);

    if 0 == irank {
        println!("MPI initialized with {nrank} ranks.");

        let version = Version::new();
        println!("{:#}", version);

        if nrank != 1 {
            bail!("MPI parallelism not implemented yet, exiting ...");
        }
    }

    let now = std::time::Instant::now();
    init_from_env(env_logger::Env::new().filter_or("RSGRAD_LOG", "info"));

    {
        use clap::Parser;
        use namd::commands;
        use namd::commands::OptProcess;
        use namd::commands::Command::*;
        match commands::Command::parse() {
            Run(run) => run.process()?,
        }
    }

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
