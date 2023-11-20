use shared::{
    //ndarray,
    Result,
    Context,
    bail,
};
use mpi::{
    self,
    topology::Communicator,
};

use namd::version::Version;

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


    Ok(())
}
