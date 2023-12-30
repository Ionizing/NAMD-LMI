use std::io::IsTerminal;

use shared::{
    //ndarray,
    Result,
    Context,
    bail,
    info,
    debug,
    tracing::{self, Level},
};
use mpi::{
    self,
    topology::Communicator,
};

use tracing_subscriber::{
    fmt,
    layer::SubscriberExt,
    FmtSubscriber,
};
use tracing_appender::{
    self,
    non_blocking::WorkerGuard,
};

use namd::version::Version;

fn main() -> Result<()> {
    let universe = mpi::initialize()
        .with_context(|| "MPI initialization failed!")?;
    let world    = universe.world();
    let nrank    = world.size();
    let irank    = world.rank();
    //let root_rank = world.process_at_rank(0);

    let _guard = init_tracing();

    if nrank != 1 {
        if 0 == irank {
            bail!("MPI parallelism not implemented yet, exiting ...");
        } else {
            return Ok(());
        }
    }

    if 0 == irank {
        println!("MPI initialized with {nrank} ranks.");

        let version = Version::new();
        info!("\n{:#}", version);
    }

    let now = std::time::Instant::now();

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


fn init_tracing() -> WorkerGuard {
    let file_appender = tracing_appender::rolling::never("./", "runlog.txt");
    let (file_writer, guard) = tracing_appender::non_blocking(file_appender);

    let isatty = std::io::stderr().is_terminal();

    let subscriber = FmtSubscriber::builder()
        .with_max_level(Level::INFO)
        .with_writer(std::io::stderr)
        .with_ansi(isatty)
        .finish()
        .with(fmt::Layer::default().with_ansi(false).with_writer(file_writer));

    tracing::subscriber::set_global_default(subscriber)
        .expect("Setting global tracing subscriber failed.");

    debug!("Tracing initialized.");

    return guard;
}
