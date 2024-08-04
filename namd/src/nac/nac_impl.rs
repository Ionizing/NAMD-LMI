use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};

use hdf5::File as H5File;
use rayon::prelude::*;

#[cfg(not(test))]
use shared::{info, anyhow::ensure};
use shared::{
    c64,
    ndarray as nd,
    Context,
    Result,
    //ndarray_linalg as nl
};
//use nl::Norm;
#[cfg(test)]
use std::{println as info, assert as ensure};

use vasp_parsers::{procar::Procar, Wavecar, WavecarType};

use crate::core::Couplings;
use crate::nac::config::NacConfig;

/// All the indices are counted from 1, and use closed inverval.
#[derive(Clone, PartialEq, Debug)]
pub struct Nac {
    ikpoint: usize,         // In this struct, ikpoint starts from 1.
    nspin: usize,
    nbands: usize,
    ndigit: usize,

    /// stores `(brange[0] ..= brange[1])` where `brange[1]` is included
    brange: [usize; 2],
    nbrange: usize,
    nsw: usize,
    efermi: f64,
    potim: f64,
    temperature: f64,
    phasecorrection: bool,

    olaps: nd::Array4<c64>, // istep, ispin, iband, iband
    eigs: nd::Array3<f64>,  // istep, ispin, iband
    pij: nd::Array5<c64>,   // istep, ispin, ixyz, iband, iband

    // projections in PROCAR
    proj: nd::Array5<f64>, // istep ispin, iband, iion, iorbit
}

struct CoupIjRet {
    c_ij: nd::Array3<c64>,
    e_ij: nd::Array2<f64>,
    p_ij: nd::Array4<c64>,
    proj: nd::Array4<f64>,
    efermi: f64,
}

struct CoupTotRet {
    olaps: nd::Array4<c64>,
    eigs: nd::Array3<f64>,
    pij: nd::Array5<c64>,
    proj: nd::Array5<f64>,
    efermi: f64,
}

impl Couplings for Nac {
    type TdCoupType<'a> = nd::ArrayView4<'a, c64>;
    type TdPijType<'a>  = nd::ArrayView5<'a, c64>;
    type TdProjType<'a> = nd::ArrayView5<'a, f64>;
    type TdEigsType<'a> = nd::ArrayView3<'a, f64>;
    type ConfigType = NacConfig;

    fn get_nspin(&self) -> usize { self.nspin }
    fn get_nbands(&self) -> usize { self.nbands }
    fn get_ikpoint(&self) -> usize { self.ikpoint }
    fn get_brange(&self) -> [usize; 2] { self.brange }
    fn get_nsw(&self) -> usize { self.nsw }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_temperature(&self) -> f64 { self.temperature }
    fn get_efermi(&self) -> f64 { self.efermi }
    fn get_tdcoup<'a>(&'a self) -> Self::TdCoupType<'a> { self.olaps.view() }
    fn get_tdpij<'a>(&'a self) -> Self::TdPijType<'a> { self.pij.view() }
    fn get_tdrij<'a>(&'a self) -> Self::TdPijType<'a> { self.pij.view() }
    fn get_tdproj<'a>(&'a self) -> Self::TdProjType<'a> { self.proj.view() }
    fn get_tdeigs<'a>(&'a self) -> Self::TdEigsType<'a> { self.eigs.view() }

    fn from_config(cfg: &NacConfig) -> Result<Self> {
        if cfg.get_nacfname().is_file() {
            info!("Found pre-calculated NAC available in {:?}, reading NAC from it ...", cfg.get_nacfname());
            let nac = Self::from_h5(cfg.get_nacfname())?;
            ensure!(cfg.get_ikpoint() == nac.ikpoint,     "Inconsistent ikpoint from config and NAC file.");
            ensure!(cfg.get_brange()  == nac.brange,      "Inconsistent brange from config and NAC file.");
            ensure!(cfg.get_nsw()     == nac.nsw,         "Inconsistent nsw from config and NAC file.");
            ensure!(cfg.get_potim()   == nac.potim,       "Inconsistent potim from config and NAC file.");
            ensure!(cfg.get_phasecorrection() == nac.phasecorrection, "Inconsistent phasecorrection from config and NAC file.");
            return Ok(nac);
        }

        Self::calculate_from_scratch(&cfg)
    }

    fn from_h5<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        let f = H5File::open(fname)?;

        // ikpoint counts from 1
        let ikpoint = f.dataset("ikpoint")?.read_scalar::<usize>()?;
        let nspin   = f.dataset("nspin")?.read_scalar::<usize>()?;
        let nbands  = f.dataset("nbands")?.read_scalar::<usize>()?;
        let ndigit  = f.dataset("ndigit")?.read_scalar::<usize>()?;
        let brange  = f.dataset("brange")?.read_scalar::<[usize;2]>()?;
        let nbrange = f.dataset("nbrange")?.read_scalar::<usize>()?;
        let nsw     = f.dataset("nsw")?.read_scalar::<usize>()?;
        let efermi  = f.dataset("efermi")?.read_scalar::<f64>()?;
        let potim   = f.dataset("potim")?.read_scalar::<f64>()?;
        let temperature = f.dataset("temperature")?.read_scalar::<f64>()?;
        let phasecorrection = f.dataset("phasecorrection")?.read_scalar::<bool>()?;

        let olaps = {
            let olaps_r: nd::Array4<f64> = f.dataset("olaps_r")?.read()?;
            let olaps_i: nd::Array4<f64> = f.dataset("olaps_i")?.read()?;
            olaps_r.mapv(|v| c64::new(v, 0.0)) + olaps_i.mapv(|v| c64::new(0.0, v))
        };

        let eigs: nd::Array3<f64> = f.dataset("eigs")?.read()?;

        let pij = {
            let pij_r: nd::Array5<f64> = f.dataset("pij_r")?.read()?;
            let pij_i: nd::Array5<f64> = f.dataset("pij_i")?.read()?;
            pij_r.mapv(|v| c64::new(v, 0.0)) + pij_i.mapv(|v| c64::new(0.0, v))
        };

        let proj: nd::Array5<f64> = f.dataset("proj")?.read()?;

        Ok(Self {
            ikpoint,
            nspin,
            nbands,
            ndigit,
            brange,
            nbrange,
            nsw,
            efermi,
            potim,
            temperature,
            phasecorrection,
            olaps,
            eigs,
            pij,
            proj,
        })
    }

    fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        let f = H5File::create(fname)?;

        // In NAC.h5, ikpoint counts from 1, but in struct, ikpoint starts from 0.
        f.new_dataset::<usize>().create("ikpoint")?.write_scalar(&self.ikpoint)?;
        f.new_dataset::<usize>().create("nspin")?.write_scalar(&self.nspin)?;
        f.new_dataset::<usize>().create("nbands")?.write_scalar(&self.nbands)?;
        f.new_dataset::<usize>().create("ndigit")?.write_scalar(&self.ndigit)?;
        f.new_dataset::<[usize;2]>().create("brange")?.write_scalar(&self.brange)?;
        f.new_dataset::<usize>().create("nbrange")?.write_scalar(&self.nbrange)?;
        f.new_dataset::<usize>().create("nsw")?.write_scalar(&self.nsw)?;
        f.new_dataset::<f64>().create("efermi")?.write_scalar(&self.efermi)?;
        f.new_dataset::<f64>().create("potim")?.write_scalar(&self.potim)?;
        f.new_dataset::<f64>().create("temperature")?.write_scalar(&self.temperature)?;
        f.new_dataset::<bool>().create("phasecorrection")?.write_scalar(&self.phasecorrection)?;

        f.new_dataset_builder().with_data(&self.olaps.mapv(|v| v.re)).create("olaps_r")?;
        f.new_dataset_builder().with_data(&self.olaps.mapv(|v| v.im)).create("olaps_i")?;

        f.new_dataset_builder().with_data(&self.eigs).create("eigs")?;

        f.new_dataset_builder().with_data(&self.pij.mapv(|v| v.re)).create("pij_r")?;
        f.new_dataset_builder().with_data(&self.pij.mapv(|v| v.im)).create("pij_i")?;

        f.new_dataset_builder().with_data(&self.proj).create("proj")?;

        Ok(())
    }
}


impl Nac {
    fn calculate_from_scratch(cfg: &NacConfig) -> Result<Self> {
        info!("No pre-calculated NAC available, start calculating from scratch in {:?}/.../WAVECARs ...", cfg.get_rundir());
        let rundir  = Path::new(cfg.get_rundir());
        let nsw     = cfg.get_nsw();
        let ikpoint = cfg.get_ikpoint() - 1;        // count from 0
        let brange  = Range { start: cfg.get_brange()[0] - 1, end: cfg.get_brange()[1] };
        let nbrange = brange.len();
        let ndigit  = cfg.get_ndigit();
        let potim   = cfg.get_potim();
        let temperature = cfg.get_temperature();
        let phasecorrection = cfg.get_phasecorrection();

        let path_1 = rundir.join(format!("{:0ndigit$}", 1)).join("WAVECAR");
        let w1 = Wavecar::from_file(&path_1)
            .with_context(|| format!("Failed to parse {:?} as WAVECAR.", &path_1))
            .unwrap();

        let nspin   = w1.nspin as usize;
        let nbands  = w1.nbands as usize;
        let nspinor = match w1.wavecar_type {
            WavecarType::NonCollinear => 2,
            _ => 1usize,
        };
        let nplw = w1.nplws[ikpoint] as usize;

        let phi_1s = {
            let mut phi = nd::Array3::<c64>::zeros((nspin, nbrange, nplw));

            for ispin in 0 .. nspin {
                for iband in brange.clone().into_iter() {
                    phi.slice_mut(nd::s![ispin, iband - brange.start, ..]).assign(
                        &w1._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                            .into_shape((nplw,))
                            .with_context(|| format!("Wavefunction reshape failed."))?
                    );
                }
            }

            phi
        };

        let procar_1 = rundir.join(format!("{:0ndigit$}", 1)).join("PROCAR");
        let p1 = Procar::from_file(&procar_1)
            .with_context(|| format!("Failed to parse {:?} as PROCAR.", &procar_1))
            .unwrap();
        let nions = p1.pdos.nions as usize;
        let nproj = p1.pdos.projected.shape()[4];

        let gvecs = nd::arr2(&w1.generate_fft_grid_cart(ikpoint as u64))
            .rows()
            .into_iter()
            .map(|g| [
                c64::new(g[0], 0.0),
                c64::new(g[1], 0.0),
                c64::new(g[2], 0.0),
            ])
            .cycle()
            .take(nplw)
            .flatten()
            .collect::<nd::Array1<c64>>()
            .into_shape((nplw, 3))
            .unwrap();

        let lncl = match w1.wavecar_type {
            WavecarType::NonCollinear => true,
            _ => false,
        };

        let CoupTotRet {olaps, eigs, pij, proj, efermi} = Self::from_wavecars(
            &phi_1s, &rundir, nsw, ikpoint, brange, ndigit, nspin, lncl, nions, nproj, &gvecs, phasecorrection
        )?;

        Ok(Self {
            ikpoint: ikpoint + 1,
            nspin,
            nbands,
            ndigit,
            brange: cfg.get_brange(),
            nbrange,
            nsw,
            efermi,
            potim,
            temperature,
            phasecorrection,

            olaps,
            eigs,
            pij,

            proj,
        })
    }

    // This function calculates non-adiabatic coupling (NAC), and transition dipole moment (TDM)
    fn from_wavecars(
        phi_1s: &nd::Array3<c64>,
        rundir: &Path, nsw: usize, ikpoint: usize, brange: Range<usize>, ndigit: usize,
        nspin: usize, lncl: bool, nions: usize, nproj: usize, gvecs: &nd::Array2<c64>, phasecorrection: bool,
        ) -> Result<CoupTotRet>
    {
        let nbrange = brange.clone().count();

        let ret_c_ij = Arc::new(Mutex::new(
                nd::Array4::<c64>::zeros((nsw-1, nspin, nbrange, nbrange))
                ));
        let ret_e_ij = Arc::new(Mutex::new(
                nd::Array3::<f64>::zeros((nsw-1, nspin, nbrange))
                ));
        let ret_p_ij = Arc::new(Mutex::new(
                nd::Array5::<c64>::zeros((nsw-1, nspin, 3, nbrange, nbrange))
                ));

        let nspinors = if lncl { 4 } else { nspin };
        let ret_proj = Arc::new(Mutex::new(
                nd::Array5::<f64>::zeros((nsw-1, nspinors, nbrange, nions, nproj))
                ));
        let efermi_sum = Arc::new(Mutex::new(0.0f64));

        let remain_count = Arc::new(Mutex::new(nsw - 1));

        (0 .. nsw-1).into_par_iter().for_each(|isw| {
            let path_i = rundir.join(format!("{:0ndigit$}", isw + 1));
            let path_j = rundir.join(format!("{:0ndigit$}", isw + 2));

            {
                let mut remain_now = remain_count.lock().unwrap();
                let remains: usize = *remain_now;
                info!(" Calculating couplings between {:?} and {:?} ..., remains: {:6}", &path_i, &path_j, &remains);
                *remain_now -= 1;
            }

            let CoupIjRet {c_ij, e_ij, p_ij, proj, efermi} = Self::coupling_ij(
                phi_1s, path_i.as_path(), path_j.as_path(), ikpoint, brange.clone(), gvecs, phasecorrection
                )
                .with_context(|| format!("Failed to calculate couplings between {:?} and {:?}.", &path_i, &path_j))
                .unwrap();

            ret_c_ij.lock().unwrap()
                .slice_mut(nd::s![isw, .., .., ..]).assign(&c_ij);

            ret_e_ij.lock().unwrap()
                .slice_mut(nd::s![isw, .., ..]).assign(&e_ij);

            ret_p_ij.lock().unwrap()
                .slice_mut(nd::s![isw, .., .., .., ..]).assign(&p_ij);

            ret_proj.lock().unwrap()
                .slice_mut(nd::s![isw, .., .., .., ..]).assign(&proj);

            let mut sum = efermi_sum.lock().unwrap();
            *sum += efermi;
        });

        Ok( CoupTotRet {
            olaps: Arc::try_unwrap(ret_c_ij).unwrap().into_inner()?,
            eigs: Arc::try_unwrap(ret_e_ij).unwrap().into_inner()?,
            pij: Arc::try_unwrap(ret_p_ij).unwrap().into_inner()?,
            proj: Arc::try_unwrap(ret_proj).unwrap().into_inner()?,
            efermi: Arc::try_unwrap(efermi_sum).unwrap().into_inner()? / (nsw - 1) as f64,
        })

    }


    fn coupling_ij(phi_1s: &nd::Array3<c64>,
        path_i: &Path, path_j: &Path,
        ikpoint: usize, brange: Range<usize>,
        gvecs: &nd::Array2<c64>, phasecorrection: bool)
        -> Result<CoupIjRet>
    {
        let wi = Wavecar::from_file(&path_i.join("WAVECAR"))?;
        let wj = Wavecar::from_file(&path_j.join("WAVECAR"))?;

        let nspinor  = match wi.wavecar_type {
            WavecarType::NonCollinear => 2,
            _ => 1usize,
        };

        let nplw     = wi.nplws[ikpoint] as usize;
        let nbrange  = brange.clone().count();
        let nspin    = wi.nspin as usize;

        assert_eq!(nplw, gvecs.shape()[0]);

        let mut c_ij = nd::Array3::<c64>::zeros((nspin, nbrange, nbrange));
        let mut e_ij = nd::Array2::<f64>::zeros((nspin, nbrange));
        let mut p_ij = nd::Array4::<c64>::zeros((nspin, 3, nbrange, nbrange));

        let mut phi_i = nd::Array2::<c64>::zeros((nbrange, nplw));
        let mut phi_j = phi_i.clone();

        let mut phase_i = nd::Array1::<c64>::zeros(nbrange);
        let mut phase_j = phase_i.clone();

        let eigs_i = wi.band_eigs.slice(nd::s![.., ikpoint, brange.clone()]).to_owned();
        let eigs_j = wj.band_eigs.slice(nd::s![.., ikpoint, brange.clone()]).to_owned();

        for ispin in 0 .. nspin {
            for iband in brange.clone().into_iter() {
                phi_i.slice_mut(nd::s![iband - brange.start, ..]).assign(
                    &wi._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .with_context(|| format!("Wavefunction reshape failed."))?
                );
                phi_j.slice_mut(nd::s![iband - brange.start, ..]).assign(
                    &wj._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .with_context(|| format!("Wavefunction reshape failed."))?
                );
            }

            // phase correction:
            //             < phi_0 | phi_i >
            // phase_0i = -------------------
            //            |< phi_0 | phi_i >|
            // phi_i[:] *= conj(phase_0i)
            if phasecorrection {
                phase_i.assign(&(phi_1s.slice(nd::s![ispin, .., ..])
                                       .mapv(|x| x.conj()) * &phi_i)
                                       .sum_axis(nd::Axis(1)));
                phase_i.mapv_inplace(|x| x.conj() / x.norm());
                phi_i.axis_iter_mut(nd::Axis(0)).zip(phase_i.iter())
                     .par_bridge()
                     .for_each(|(mut row, phase)| row.mapv_inplace(|x| x * phase));

                phase_j.assign(&(phi_1s.slice(nd::s![ispin, .., ..])
                                       .mapv(|x| x.conj()) * &phi_j)
                                       .sum_axis(nd::Axis(1)));
                phase_j.mapv_inplace(|x| x.conj() / x.norm());
                phi_j.axis_iter_mut(nd::Axis(0)).zip(phase_j.iter())
                     .par_bridge()
                     .for_each(|(mut row, phase)| row.mapv_inplace(|x| x * phase));
            }


            let phi_ij = phi_i.mapv(|v| v.conj()).dot(&phi_j.t());
            let phi_ji = phi_ij.t().mapv(|v| v.conj());  // phi_ji = phi_ij^H
            c_ij.slice_mut(nd::s![ispin, .., ..]).assign(&(phi_ij - phi_ji));

            for idirect in 0 .. 3 {
                let phi_x_gvecs: nd::Array2<_> = phi_j.clone() * gvecs.slice(nd::s![nd::NewAxis, .., idirect]);

                // <i | p | j>, in eV*fs/Angstrom
                let p_ij_tmp = match wi.wavecar_type {
                    WavecarType::GammaHalf(_) => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t())
                                               - phi_x_gvecs.mapv(|v| v.conj()).dot(&phi_j.t()),
                    _ => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t()),
                };
                p_ij.slice_mut(nd::s![ispin, idirect, .., ..]).assign(&p_ij_tmp);
            }
        }

        // read PROCAR
        let proj = {
            let fpath = path_i.join("PROCAR");
            let proj_i = Procar::from_file(&fpath)
                .with_context(|| format!("Failed to parse {:?}.", &fpath))
                .unwrap();
            proj_i.pdos.projected.slice(nd::s![.., ikpoint, brange.clone(), .., ..]).to_owned()
        };

        e_ij.slice_mut(nd::s![.., ..]).assign(&(
            ( eigs_i + eigs_j ) / 2.0
        ));

        Ok( CoupIjRet {
            c_ij,
            e_ij, 
            p_ij, 
            proj,
            efermi: wi.efermi,
        })
    }

    pub fn get_ndigit(&self) -> usize { self.ndigit }
}
