use std::path::Path;
use std::path::PathBuf;
use std::fs::create_dir_all;

use rayon::prelude::*;

use itertools::iproduct;
use hdf5::File as H5File;
use clap::{Parser, ValueEnum};
use shared::{
    ndarray as nd,
    ndrustfft as fft,
    c64,
    Context,
    Result,
    bail,
    log,
    copy_file_to,
    rfft_freq_1d,
    auto_correlation,
    lower_triangle_matrix_index,
};
use crate::OptProcess;
use crate::cli::write_script;
use crate::version::Version;

use crate::core::{
    Hamiltonian,
    SurfaceHopping,
    NamdConfig,
};
use crate::surfhop;
use crate::hamil;


#[derive(Debug, Parser)]
/// Perform the surface-hopping process with given Hamiltonian file and config file.
#[command(arg_required_else_help(true))]
pub struct SurfhopCommand {
    #[arg(short='n', long, default_value_t=0)]
    /// Number of threads for parallel calculation.
    /// 
    /// If 0 is set, it will fall back to the number of logic CPU cores of you machine.
    nthreads: usize,

    #[arg(short='c', long, aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: Option<PathBuf>,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation and analysis.
    ///
    /// The surface-hopping will not run if this flag is set.
    ///
    /// Alias: "gen".
    generate: Option<TemplateGenerator>,

    #[arg(long, aliases=["collect", "cr"])]
    /// Collect results produced by the surface-hopping.
    ///
    /// The surface-hopping will not run if this flag is set.
    ///
    /// Aliases: "collect", "cr"
    collect_results: Option<PathBuf>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum TemplateGenerator {
    #[value(aliases=["config", "cfg", "conf"])]
    /// Generate config template for Hamiltonian generation.
    /// Aliases: "config", "cfg" and "conf".
    ConfigTemplate,

    #[value(name="inistep.py", aliases=["inistep", "ini"])]
    /// Generate Python script to help append `inistep` field.
    /// Aliases: "inistep" and "ini"
    Inistep,

    #[value(aliases=["post-process", "postprocess", "pp"])]
    /// Generate post-process scripts for surface-hopping analysis.
    /// Aliases: "post-process", "postprocess", "pp".
    PostprocessTemplate,
}


impl OptProcess for SurfhopCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator as TG;

        if let Some(g) = self.generate {
            return match g {
                TG::ConfigTemplate => {
                    {
                        log::info!("writing `03_surfhop_config_template.toml` ...");
                        surfhop::SurfhopConfig::default().to_file("03_surfhop_config_template.toml")
                    }.and_then(|_| {
                        log::info!("Writing `inisteps.py` ...");
                        surfhop::SurfhopConfig::write_inistep_py("inisteps.py")
                    })
                },
                TG::Inistep => {
                    log::info!("Writing `inisteps.py` ...");
                    surfhop::SurfhopConfig::write_inistep_py("inisteps.py")
                },
                TG::PostprocessTemplate => {
                    log::info!("Writing `surfhop_plot.py` ...");
                    write_script("surfhop_plot.py", include_str!("./surfhop_plot.py"), true)
                },
            }
        }

        let config_fname = self.config.clone()
            .context("A connfig file is required via `-c surfhop_config.toml`.")?;

        if let Some(outdir) = self.collect_results.to_owned() {
            log::info!("Collecting results from existing surface hopping artifact with {} threads ...", self.nthreads);
            rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();
            let cfg = surfhop::SurfhopConfig::from_file(&config_fname)?;
            collect_results(&cfg, outdir, "averaged_results.h5")?;
            return Ok(())
        }


        // Start running surface-hopping method.
        log::info!("\n{}", Version::new());
        log::info!("Prepare to run surface-hopping in {} threads ...", self.nthreads);
        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let mut cfg = surfhop::SurfhopConfig::from_file(&config_fname)?;
        create_outputdir(cfg.get_outdir_mut())?;
        let cfg = cfg;      // cancel mutability

        crate::logging::logger_redirect(&cfg.get_outdir())?;
        copy_file_to(&config_fname, cfg.get_outdir())?;

        log::info!("Got Surface Hopping config:\n{}", &cfg);
        let sh = surfhop::Surfhop::from_config(&cfg)?;
        let _ = sh.into_par_iter()
            .map(|mut v| v.run())
            .collect::<Result<Vec<()>>>()?;

        log::info!("Collecting results ...");
        collect_results(&cfg, cfg.get_outdir(), "averaged_results.h5")?;
        Ok(())
    }
}


fn create_outputdir(dir: &mut PathBuf) -> Result<()> {
    if dir.is_file() {
        bail!("The output dir {:?} exists as a regular file, please change.", dir);
    }

    if dir.file_name().is_none() {
        bail!("The output dir {:?} cannot be current working dir, please change.", dir);
    }

    if dir.is_dir() {
        let parent = dir.parent().unwrap();
        let subdir = dir.file_name().unwrap().to_str().unwrap();
        let mut newdir: Option<PathBuf> = None;
        let mut tmpdir = PathBuf::new();

        for i in 1 ..= 99 {
            let dirstr = format!("{}_{:02}", &subdir, i);
            tmpdir = parent.join(&dirstr);
            if !tmpdir.is_file() && !tmpdir.is_dir() {
                newdir = Some(tmpdir.clone());
                break;
            }
        }

        if let Some(newdir) = newdir {
            log::warn!("The outdir {:?} is already exists and will be switched to {:?} for this run.", dir, newdir);
            *dir = newdir;
        } else {
            bail!("Existed outdir reached maximum homonymy outdirs: {:?}", tmpdir);
        }
    }

    log::info!("Log and output files will be stored in {:?} .", dir);
    create_dir_all(dir)?;

    Ok(())
}


fn collect_results<P1: AsRef<Path>, P2: AsRef<Path>>(
    cfg: &surfhop::SurfhopConfig,
    outdir: P1,
    collected_fname: P2) -> Result<()>
{
    struct ResultType {
        time:           nd::Array1<f64>, // [namdtime]
        prop_energy:    nd::Array1<f64>, // [namdtime]
        psi_t:          nd::Array2<f64>, // [namdtime, nbasis]
        sh_energy:      nd::Array1<f64>, // [namdtime]
        sh_pops:        nd::Array2<f64>, // [namdtime, nbasis]
        eigs_t:         nd::Array2<f64>, // [namdtime, nbasis]
        proj_t:         nd::Array4<f64>, // [namdtime, nbasis, nions, nproj]
        sh_phonons_t:   nd::Array3<f64>, // [namdtime, nbasis, nbasis]
        sh_photons_t:   nd::Array3<f64>, // [namdtime, nbasis, nbasis]
        photons_emit_t: nd::Array2<f64>, // [namdtime, npoints], May it's better to plot theses things in python
        photons_absp_t: nd::Array2<f64>, // [namdtime, npoints]
        phonons_emit_t: nd::Array2<f64>, // [namdtime, nfreqs]
        phonons_absp_t: nd::Array2<f64>, // [namdtime, nfreqs]
    }

    let outdir = outdir.as_ref();

    let smearing_method = cfg.get_smearing_method();
    let smearing_sigma  = cfg.get_smearing_sigma();

    let hamil_fname = cfg.get_hamil_fname();
    let hamil    = hamil::SPHamiltonian::from_h5(hamil_fname)?;
    let ndigit   = hamil.get_ndigit();
    let potim    = hamil.get_potim();
    let namdtime = cfg.get_namdtime();

    let eigs_t = hamil.get_eigs_t();
    let proj_t = hamil.get_proj();
    let (nsw, nbasis, nions, nproj) = proj_t.dim();

    let (delta_et, frequencies, spectra) = get_transition_phonons(hamil.get_eigs_t(), potim);
    let xvals   = generate_xvals_linspace(delta_et.view(), cfg.get_npoints_per_ev());
    let npoints = xvals.len();
    let nfreqs  = frequencies.len();

    let result_sum = cfg.get_inisteps().par_iter()
        .map(|&namdinit| -> Result<ResultType> {
            let fname = outdir.join(format!("result_{:0ndigit$}.h5", namdinit));
            log::info!("Processing {:?} ...", fname);
            let f = H5File::open(fname)?;

            let time:        nd::Array1<f64> = f.dataset("time")?.read()?;
            let prop_energy: nd::Array1<f64> = f.dataset("prop_energy")?.read()?;
            let psi_t_r:     nd::Array2<f64> = f.dataset("psi_t_r")?.read()?;
            let psi_t_i:     nd::Array2<f64> = f.dataset("psi_t_i")?.read()?;
            let psi_t = psi_t_r.map(|x| x*x) + psi_t_i.map(|x| x*x);
            let sh_energy: nd::Array1<f64> = f.dataset("sh_energy")?.read()?;
            let sh_pops:   nd::Array2<f64> = f.dataset("sh_pops")?.read()?;

            let mut eigs = nd::Array2::<f64>::zeros((namdtime, nbasis));
            let mut proj = nd::Array4::<f64>::zeros((namdtime, nbasis, nions, nproj));

            let tdphotons: nd::Array3<f64> = f.dataset("sh_photons_t")?.read()?;
            let mut photons_emit_t = nd::Array2::<f64>::zeros((namdtime, npoints));
            let mut photons_absp_t = nd::Array2::<f64>::zeros((namdtime, npoints));

            let tdphonons: nd::Array3<f64> = f.dataset("sh_phonons_t")?.read()?;
            let mut phonons_emit_t = nd::Array2::<f64>::zeros((namdtime, nfreqs));
            let mut phonons_absp_t = nd::Array2::<f64>::zeros((namdtime, nfreqs));

            for iion in 0 .. namdtime {
                let [rtime, _xtime] = hamil::SPHamiltonian::get_rtime_xtime(iion, nsw, namdinit);
                eigs.slice_mut(nd::s![iion, ..]).assign(&eigs_t.slice(nd::s![rtime, ..]));
                proj.slice_mut(nd::s![iion, .., .., ..]).assign(&proj_t.slice(nd::s![rtime, .., .., ..]));

                let mut photon_emit_de  = vec![0.0f64; 0];
                let mut photon_emit_num = vec![0.0f64; 0];
                let mut photon_absp_de  = vec![0.0f64; 0];
                let mut photon_absp_num = vec![0.0f64; 0];
                for (iband, jband) in iproduct!(0 .. nbasis, 0 .. nbasis) {
                    if iband == jband { continue; }
                    let de = (eigs[(iion, iband)] - eigs[(iion, jband)]).abs();
                    let oc = tdphotons[(iion, iband, jband)];

                    if oc < 0.0 {
                        photon_absp_de.push(de);        // absorption part
                        photon_absp_num.push(oc);
                    } else if oc == 0.0 {
                        continue;
                    } else {
                        photon_emit_de.push(de);        // emission part
                        photon_emit_num.push(oc);
                    }
                }

                let photons_emit: nd::Array1<f64> = smearing_method.apply_smearing(
                        &xvals.as_slice().unwrap(),
                        &photon_emit_de,
                        smearing_sigma,
                        Some(&photon_emit_num)
                    );
                let photons_absp: nd::Array1<f64> = smearing_method.apply_smearing(
                        &xvals.as_slice().unwrap(),
                        &photon_absp_de,
                        smearing_sigma,
                        Some(&photon_absp_num)
                    );

                photons_emit_t.row_mut(iion).assign(&photons_emit);
                photons_absp_t.row_mut(iion).assign(&photons_absp);


                let mut phonons_emit = nd::Array1::<f64>::zeros(nfreqs);
                let mut phonons_absp = nd::Array1::<f64>::zeros(nfreqs);
                for (iband, jband) in iproduct!(0 .. nbasis, 0 .. nbasis) {
                    if iband == jband { continue; }
                    let oc = tdphonons[(iion, iband, jband)];

                    let idx = lower_triangle_matrix_index(iband, jband);
                    if oc < 0.0 {
                        phonons_absp += &(spectra.row(idx).to_owned() * oc);
                    } else if oc == 0.0 {
                        continue;
                    } else {
                        phonons_emit += &(spectra.row(idx).to_owned() * oc);
                    }
                }

                phonons_emit_t.row_mut(iion).assign(&phonons_emit);
                phonons_absp_t.row_mut(iion).assign(&phonons_absp);
            }

            Ok(ResultType {
                time,
                prop_energy,
                psi_t,
                sh_energy,
                sh_pops,
                eigs_t: eigs,
                proj_t: proj,
                sh_phonons_t: tdphonons,
                sh_photons_t: tdphotons,
                photons_emit_t,
                photons_absp_t,
                phonons_emit_t,
                phonons_absp_t,
            })
        })
        .reduce_with(|acc, e| -> Result<ResultType> {
            let acc = acc?;
            let e = e?;
            Ok(ResultType {
                time:        e.time,

                prop_energy: acc.prop_energy + e.prop_energy,
                psi_t:       acc.psi_t     + e.psi_t,
                sh_energy:   acc.sh_energy + e.sh_energy,
                sh_pops:     acc.sh_pops   + e.sh_pops,

                eigs_t:      acc.eigs_t + e.eigs_t,
                proj_t:      acc.proj_t + e.proj_t,

                sh_phonons_t: acc.sh_phonons_t + e.sh_phonons_t,
                sh_photons_t: acc.sh_photons_t + e.sh_photons_t,

                photons_emit_t: acc.photons_emit_t + e.photons_emit_t,
                photons_absp_t: acc.photons_absp_t + e.photons_absp_t,

                phonons_emit_t: acc.phonons_emit_t + e.phonons_emit_t,
                phonons_absp_t: acc.phonons_absp_t + e.phonons_absp_t,
            })
        })
        .context("No results collected")??;

    let time    = result_sum.time;
    let nsample = cfg.get_inisteps().len() as f64;

    let prop_energy = result_sum.prop_energy / nsample;
    let psi_t     = result_sum.psi_t / nsample;
    let sh_energy = result_sum.sh_energy / nsample;
    let sh_pops   = result_sum.sh_pops / nsample;
    let eigs_t    = result_sum.eigs_t / nsample;
    let proj_t    = result_sum.proj_t / nsample;
    let sh_phonons_t = result_sum.sh_phonons_t / nsample;
    let sh_photons_t = result_sum.sh_photons_t / nsample;
    let photons_emit_t = result_sum.photons_emit_t / nsample;
    let photons_absp_t = result_sum.photons_absp_t / nsample;
    let phonons_emit_t = result_sum.phonons_emit_t / nsample;
    let phonons_absp_t = result_sum.phonons_absp_t / nsample;

    // Writing results.
    let ret_fname = outdir.join(collected_fname);
    log::info!("Collecting done. Writing to {:?} ...", &ret_fname);

    let (basis_list, basis_labels_src) = {
        let namdinit = cfg.get_inisteps()[0];
        let fname = outdir.join(format!("result_{:0ndigit$}.h5", namdinit));
        let f = H5File::open(fname)?;

        let basis_list: Vec<i32> = f.dataset("basis_list")?.read_raw::<i32>()?;
        let basis_labels = if f.dataset("basis_labels").is_err() {
            None
        } else {
            let raw: Vec<u8> = f.dataset("basis_labels")?.read_raw()?;
            Some(raw)
        };
        (basis_list, basis_labels)
    };

    let f = H5File::create(ret_fname)?;
    f.new_dataset::<usize>().create("ndigit")?.write_scalar(&ndigit)?;
    f.new_dataset::<f64>().create("potim")?.write_scalar(&potim)?;

    f.new_dataset_builder().with_data(&time).create("time")?;
    f.new_dataset_builder().with_data(&basis_list).create("basis_list")?;
    if let Some(src) = basis_labels_src {
        f.new_dataset_builder().with_data(&src).create("basis_labels")?;
    }
    
    f.new_dataset_builder().with_data(&prop_energy).create("prop_energy")?;
    f.new_dataset_builder().with_data(&psi_t).create("psi_t")?;
    f.new_dataset_builder().with_data(&sh_energy).create("sh_energy")?;
    f.new_dataset_builder().with_data(&sh_pops).create("sh_pops")?;

    f.new_dataset_builder().with_data(&eigs_t).create("eigs_t")?;
    f.new_dataset_builder().with_data(&proj_t).create("proj_t")?;

    f.new_dataset_builder().with_data(&sh_phonons_t).create("sh_phonons_t")?;
    f.new_dataset_builder().with_data(&sh_photons_t).create("sh_photons_t")?;

    f.new_dataset_builder().with_data(&phonons_emit_t).create("phonons_emit_t")?;
    f.new_dataset_builder().with_data(&phonons_absp_t).create("phonons_absp_t")?;

    f.new_dataset_builder().with_data(&photons_emit_t).create("photons_emit_t")?;
    f.new_dataset_builder().with_data(&photons_absp_t).create("photons_absp_t")?;

    f.new_dataset_builder().with_data(&frequencies).create("phonon_spectra_frequencies")?;
    f.new_dataset_builder().with_data(&spectra).create("phonons_spectra")?;
    f.new_dataset_builder().with_data(&xvals).create("photon_spectra_xvals")?;
    f.new_dataset_builder().with_data(&delta_et).create("delta_et")?;

    Ok(())
}


// Calculate the phonon contribution by (Ej(t) - Ei(t)) |> ACF |> rFFT
//
// where
//     - Ei(t) and Ej(t) are time evolution of band energies,
//     - ACF means auto correlation function, defined in `shared::numeric_methods`,
//     - rFFT is the real-to-complex Fast Fourier Transform.
//
// Returns:
//     - delta_et: [nbasis*(nbasis+1)/2, nsw]
//     - frequencies: [nsw/2 + 1]
//     - spectra: [nbasis*(nbasis+1)/2, nfreqs], (nfreqs = nsw / 2 + 1)
//
// The spectra will have the shape of [nfreq, nbasis*(nbasis+1)/2]
//
// O . . .
// O O . .
// O O O .
// O O O O
//
fn get_transition_phonons(eigs: nd::ArrayView2<f64>, potim: f64)
    -> (nd::Array2<f64>, nd::Array1<f64>, nd::Array2<f64>)   // (Eit - Ejt, frequencies, spectra)
{
    #[allow(non_upper_case_globals)]
    const Fs2InvCm: f64 = 33_356.420_646_319_475;

    let (nsw, nbasis) = eigs.dim();
    let frequencies = rfft_freq_1d::<nd::Array1<_>>(nsw, potim) * Fs2InvCm;

    let mut acf_total = nd::Array2::<f64>::zeros((nbasis * (nbasis + 1) / 2, nsw));
    let mut delta_et  = nd::Array2::<f64>::zeros((nbasis * (nbasis + 1) / 2, nsw));
    let mut spectra   = nd::Array2::<c64>::zeros((nbasis * (nbasis + 1) / 2, nsw / 2 + 1));

    for iband in 0 .. nbasis {
        for jband in 0 .. iband {
            let idx = lower_triangle_matrix_index(iband, jband);
            let de  = eigs.column(iband).to_owned() - eigs.column(jband);
            auto_correlation(
                de.as_slice().unwrap(),
                acf_total.row_mut(idx).as_slice_mut().unwrap(),
            );
            delta_et.row_mut(idx).assign(&de);
        }

        let idx = lower_triangle_matrix_index(iband, iband);
        let e = eigs.column(iband).to_owned();
        auto_correlation(
            e.as_slice().unwrap(),
            acf_total.row_mut(idx).as_slice_mut().unwrap(),
        );
        delta_et.row_mut(idx).assign(&e);
    }

    let mut handler = fft::R2cFftHandler::<f64>::new(nsw);
    fft::ndfft_r2c_par(&acf_total, &mut spectra, &mut handler, 1);
    let spectra = spectra.mapv(|x| x.norm());

    (delta_et, frequencies, spectra)
}


fn generate_xvals_linspace<D>(delta_et: nd::ArrayView<f64, D>, npoints_per_ev: usize) -> nd::Array1<f64>
where D: nd::Dimension
{
    let de_min  = *delta_et.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap() - 0.5;
    let de_max  = *delta_et.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap() + 0.5;
    let npoints = ((de_max - de_min) * npoints_per_ev as f64).ceil() as usize + 1;
    nd::Array1::<f64>::linspace(de_min, de_max, npoints)
}
