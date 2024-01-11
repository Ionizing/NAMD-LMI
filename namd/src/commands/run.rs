use std::time::Instant;

use clap::Args;
use rayon::{
    ThreadPoolBuilder,
    prelude::*
};
use hdf5::File as H5File;
use shared::{
    anyhow::Context,
    ndarray::{
        Array1,
        Array2,
    },
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
    version::Version,
};
use crate::commands::{
    OptProcess,
    Result,
};


#[derive(Args, Debug)]
pub struct Run {
    #[arg(default_value_t = String::from("./input.toml") )]
    /// Specify the input file name for namd_lumi.
    input: String,

    #[arg(short, long, default_value_t = 0)]
    /// Set the number of threads for parallel calculation.
    nthreads: usize,

    #[arg(short, long)]
    /// Write post-processing Python scripts or not.
    write_scripts: bool,
}


impl OptProcess for Run {
    fn process(&self) -> Result<()> {
        ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let input = Input::from_file(&self.input)?;
        logger_redirect(&input.outdir)?;
        info!("\n{:#}", Version::new());
        input.print_to_log();

        info!("Copying & linking essential files for post-processing ...");
        copy_file_to(&self.input, &input.outdir)?;
        if let Some(efield) = input.efield.as_ref() {
            copy_file_to(&efield.0, &input.outdir)?;
        }
        link_file_to(&input.nacfname, &input.outdir)?;

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

        if self.write_scripts {
            info!("Writing post-processing Python scripts to result folder ...");

            use scripts::Scripts::PostProcess;
            scripts::write_script(&input.outdir, PostProcess)?;
        }

        let collected_fname = "results_avg.h5";
        info!("Collecting results to {:?} ...", input.outdir.join(collected_fname));
        collect_results(&input, collected_fname)?;
        info!("All the results are written to {:?}", input.outdir);

        Ok(())
    }
}


fn collect_results(inp: &Input, collected_fname: &str) -> Result<()> {
    type ResultType = (Array1<f64>,     // time         [namdtime]
                       Array1<f64>,     // prop_energy  [namdtime]
                       Array2<f64>,     // psi_t        [namdtime, nbasis]
                       Array1<f64>,     // sh_energy    [namdtime]
                       Array2<f64>,);   // sh_pops      [namdtime, nbasis]

    let ndigit = inp.ndigit;
    let results = inp.inisteps.par_iter()
        .map(|i| -> Result<ResultType> {
            let fname = inp.outdir.join(&format!("result_{:0ndigit$}.h5", i));
            let f = H5File::open(fname)?;

            let time:        Array1<f64> = f.dataset("time")?.read()?;
            let prop_energy: Array1<f64> = f.dataset("prop_energy")?.read()?;
            let psi_t_r:     Array2<f64> = f.dataset("psi_t_r")?.read()?;
            let psi_t_i:     Array2<f64> = f.dataset("psi_t_i")?.read()?;
            let psi_t = psi_t_r.mapv(|x| x*x) + psi_t_i.mapv(|x| x * x);
            let sh_energy:   Array1<f64> = f.dataset("sh_energy")?.read()?;
            let sh_pops:     Array2<f64> = f.dataset("shpops")?.read()?;

            return Ok((time,
                       prop_energy,
                       psi_t,
                       sh_energy,
                       sh_pops));
        })
        .reduce_with(|acc, e| -> Result<ResultType> {
            let acc = acc?;
            let (time1, prop_energy1, psi_t1, sh_energy1, sh_pops1) = e?;
            Ok((time1, acc.1+prop_energy1, acc.2+psi_t1, acc.3+sh_energy1, acc.4+sh_pops1))
        }).context("No results collected")??;

    let nresults = inp.inisteps.len() as f64;
    let (time, mut prop_energy, mut psi_t, mut sh_energy, mut shpops) = results;
    prop_energy /= nresults;
    psi_t       /= nresults;
    sh_energy   /= nresults;
    shpops      /= nresults;

    let f = H5File::create(inp.outdir.join(collected_fname))?;

    f.new_dataset_builder().with_data(&time).create("time")?;
    f.new_dataset_builder().with_data(&prop_energy).create("prop_energy")?;
    f.new_dataset_builder().with_data(&psi_t).create("psi_t")?;
    f.new_dataset_builder().with_data(&sh_energy).create("sh_energy")?;
    f.new_dataset_builder().with_data(&shpops).create("shpops")?;

    Ok(())
}
