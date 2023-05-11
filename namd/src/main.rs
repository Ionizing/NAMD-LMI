use shared::{
    //ndarray,
    Result,
    Context,
};
use mpi::{
    self,
    topology::Communicator,
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
    }

    Ok(())
}
