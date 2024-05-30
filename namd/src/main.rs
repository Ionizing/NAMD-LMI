use shared::{
    Result,
    info,
};

use namd::cli::run;

fn main() -> Result<()> {
    let now = std::time::Instant::now();

    run()?;

    info!("Time used: {:?}", now.elapsed());
    Ok(())
}
