use clap::Args;

//use shared::bail;
use rayon::{
    ThreadPoolBuilder,
    prelude::*
};
use crate::{
    input::Input,
    nac::Nac,
    hamiltonian::Hamiltonian,
    surface_hopping::SurfaceHopping,
};
use crate::commands::{
    OptProcess,
    Result,
};

#[derive(Args)]
pub struct Run {
    #[arg(default_value_t = String::from("./input.toml") )]
    /// Specify the input file name for namd_lumi.
    input: String,

    #[arg(long, default_value_t = 0)]
    /// Set the number of threads for parallel calculation.
    nthreads: usize,
}


impl OptProcess for Run {
    fn process(&self) -> Result<()> {
        ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let input = Input::from_file(&self.input)?;
        let nac = Nac::from_inp(&input)?;
        let ninibands = input.inibands.len();
        Hamiltonian::from_input(&nac, &input, 0).save_to_h5("HAMIL.h5")?;
        for iniband_idx in 0 .. ninibands {
            let hamil = Hamiltonian::from_input(&nac, &input, iniband_idx);
            let mut sh = SurfaceHopping::from_input(hamil, &input);
            sh.run();
        }
        Ok(())
    }
}
