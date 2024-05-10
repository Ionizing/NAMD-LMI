use std::thread;
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
        Array3,
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

        let nac = Nac::from_inp(&input)?;
        link_file_to(&input.nacfname, &input.outdir)?;

        let ninibands = input.inibands.len();

        let thread_join_handle = {
            let hamil = Hamiltonian::from_input(&nac, &input, 0);
            let outdir = input.outdir.clone();

            thread::spawn(move || -> Result<()> {
                hamil.save_to_h5(&outdir.join("HAMIL.h5"))?;

                if let Some(e) = hamil.efield.as_ref() {
                    info!("Writing TDEFIELD.txt ...");
                    e.print_efield_to_file(&outdir.join("TDEFIELD.txt"))?;

                    info!("Writing TDAFIELD.txt ...");
                    e.print_afield_to_file(&outdir.join("TDAFIELD.txt"))?;
                }

                Ok(())
            })
        };

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

        thread_join_handle.join().unwrap()?;

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
            Ok((time1,
                acc.1+prop_energy1,
                acc.2+psi_t1,
                acc.3+sh_energy1,
                acc.4+sh_pops1))
        })
        .context("No results collected")??;

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

    if inp.ldbg_hamil_t {
        type DebugInfoType = (Array3<f64>, Array3<f64>, Array2<f64>);

        let dbg_info = inp.inisteps.par_iter()
            .map(|i| -> Result<DebugInfoType> {
                let fname = inp.outdir.join(&format!("result_{:0ndigit$}.h5", i));
                let f = H5File::open(fname)?;

                let ham_t_r: Array3<f64> = f.dataset("ham_t_r")?.read()?;
                let ham_t_i: Array3<f64> = f.dataset("ham_t_i")?.read()?;
                let ham_t = ham_t_r.mapv(|x| x*x) + ham_t_i.mapv(|x| x*x);

                let lmi_t_r: Array3<f64> = f.dataset("lmi_t_r")?.read()?;
                let lmi_t_i: Array3<f64> = f.dataset("lmi_t_i")?.read()?;
                let lmi_t = lmi_t_r.mapv(|x| x*x) + lmi_t_i.mapv(|x| x*x);

                let lvl_t: Array2<f64> = f.dataset("lvl_t")?.read()?;

                return Ok((ham_t,
                           lmi_t,
                           lvl_t))
            })
            .reduce_with(|acc, e| -> Result<DebugInfoType> {
                let acc = acc?;
                let (ham_t, lmi_t, lvl_t) = e?;
                Ok((acc.0 + ham_t,
                    acc.1 + lmi_t,
                    acc.2 + lvl_t))
            })
            .context("No debug info collected")??;

        let (mut ham_t, mut lmi_t, mut lvl_t) = dbg_info;
        ham_t /= nresults;
        lmi_t /= nresults;
        lvl_t /= nresults;

        f.new_dataset_builder().with_data(&ham_t).create("ham_t")?;
        f.new_dataset_builder().with_data(&lmi_t).create("lmi_t")?;
        f.new_dataset_builder().with_data(&lvl_t).create("lvl_t")?;
    }

    Ok(())
}
