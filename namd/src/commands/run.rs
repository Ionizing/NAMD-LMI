use std::time::Instant;

use clap::Args;

//use shared::bail;
use rayon::{
    ThreadPoolBuilder,
    prelude::*
};
use shared::info;
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

        {
            let hamil = Hamiltonian::from_input(&nac, &input, 0);
            hamil.save_to_h5("HAMIL.h5")?;

            if let Some(e) = hamil.efield.as_ref() {
                info!("Writing TDEFIELD.txt ...");
                e.print_efield_to_file("TDEFIELD.txt")?;

                info!("Writing TDAFIELD.txt ...");
                e.print_afield_to_file("TDAFIELD.txt")?;
            }
        }

        (0 .. ninibands).into_par_iter()
            .for_each(move |iniband_idx| {
                let now = Instant::now();

                let hamil = Hamiltonian::from_input(&nac, &input, iniband_idx);
                let mut sh = SurfaceHopping::from_input(hamil, &input);
                sh.run();

                info!("Time used for {}(th/st) inicon: {:?}", iniband_idx + 1, now.elapsed());
            });
        Ok(())
    }
}
