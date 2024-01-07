use std::time::Instant;

use clap::Args;
use rayon::{
    ThreadPoolBuilder,
    prelude::*
};
use shared::{
    info,
    copy_file_to,
    link_file_to,
};
use crate::{
    commands::scripts,
    input::Input,
    nac::Nac,
    hamiltonian::Hamiltonian,
    surface_hopping::SurfaceHopping,
    logging::logger_redirect,
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

    #[arg(short, long, default_value_t = 0)]
    /// Set the number of threads for parallel calculation.
    nthreads: usize,

    #[arg(short, long)]
    write_scripts: bool,
}


impl OptProcess for Run {
    fn process(&self) -> Result<()> {
        ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let input = Input::from_file(&self.input)?;
        logger_redirect(&input.outdir)?;
        input.print_to_log();

        let nac = Nac::from_inp(&input)?;
        let ninibands = input.inibands.len();

        {
            let hamil = Hamiltonian::from_input(&nac, &input, 0);
            hamil.save_to_h5(&input.outdir.join("HAMIL.h5"))?;

            if let Some(e) = hamil.efield.as_ref() {
                info!("Writing TDEFIELD.txt ...");
                e.print_efield_to_file(&input.outdir.join("TDEFIELD.txt"))?;

                info!("Writing TDAFIELD.txt ...");
                e.print_afield_to_file(&input.outdir.join("TDAFIELD.txt"))?;
            }
        }

        (0 .. ninibands).into_par_iter()
            .for_each(|iniband_idx| {
                let now = Instant::now();

                let hamil = Hamiltonian::from_input(&nac, &input, iniband_idx);
                let mut sh = SurfaceHopping::from_input(hamil, &input);
                sh.run();

                info!("Time used for {}(th/st) inicon: {:?}", iniband_idx + 1, now.elapsed());
            });

        info!("Copying & linking essential files for post-processing ...");
        copy_file_to(&self.input, &input.outdir)?;
        if let Some(efield) = input.efield.as_ref() {
            copy_file_to(&efield.0, &input.outdir)?;
        }
        link_file_to(&input.nacfname, &input.outdir)?;

        if self.write_scripts {
            info!("Writing post-processing Python scripts to result folder ...");

            use scripts::Scripts::PostProcess;
            scripts::write_script(&input.outdir, PostProcess)?;
        }

        info!("The results are written to {:?}", input.outdir);

        Ok(())
    }
}
