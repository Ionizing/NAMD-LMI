//! This module calculates the couplings involved in NAMD.
//!
//!

use std::path::Path;
use std::ops::Range;
use std::sync::{
    Arc,
    Mutex,
};

use rayon::prelude::*;
//use mpi::traits::*;
use hdf5::File as H5File;
use pathfinding::prelude::{
    kuhn_munkres,
    Matrix as pfMatrix
};
use ordered_float::OrderedFloat;

use shared::{
    Context,
    Result,
    c64,
    ndarray::{
        s,
        arr2,
        ArrayBase,
        Array1,
        Array2,
        Array3,
        Array4,
        Array5,
        NewAxis,
        Axis,
        Data,
        //Dimension,
        Ix2,
    },
};

#[cfg(not(test))]
use shared::info;
#[cfg(test)]
use std::{println as info, println as warn};

use vasp_parsers::{
    Wavecar,
    WavecarType,
    procar::Procar,
};

use crate::input::Input;

// All the indices are counted from 1, and use closed inverval
#[derive(Clone, PartialEq, Debug)]
pub struct Nac {
    pub ikpoint: usize,
    pub nspin:   usize,
    pub nbands:  usize,

    /// stores `(brange[0] ..= brange[1])` where `brange[1]` is included
    pub brange:  [usize; 2],
    pub nbrange: usize,
    pub nsw:     usize,
    pub efermi:  f64,
    pub dt:      f64,
    pub lreal:   bool,

    pub olaps:   Array4<c64>,   // istep, ispin, iband, iband
    pub eigs:    Array3<f64>,   // istep, ispin, iband
    pub pij:     Array5<c64>,   // istep, ispin, ixyz, iband, iband

    // projections in PROCAR
    pub proj:    Array5<f64>,   // istep ispin, iband, iion, iorbit
    pub order:   Array3<usize>,
}


struct CoupIjRet {
    c_ij: Array3<c64>,
    e_ij: Array2<f64>,
    p_ij: Array4<c64>,
    proj: Array4<f64>,
    efermi: f64,
    order: Array2<usize>,
}


struct CoupTotRet {
    olaps: Array4<c64>,
    eigs: Array3<f64>,
    pij: Array5<c64>,
    proj: Array5<f64>,
    efermi: f64,
    order: Array3<usize>,
}


impl Nac {
    pub fn from_h5<P>(fname: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let f = H5File::open(fname)?;

        let ikpoint = f.dataset("ikpoint")?.read_scalar::<usize>()?;
        let nspin   = f.dataset("nspin")?.read_scalar::<usize>()?;
        let nbands  = f.dataset("nbands")?.read_scalar::<usize>()?;
        let brange  = f.dataset("brange")?.read_scalar::<[usize;2]>()?;
        let nbrange = f.dataset("nbrange")?.read_scalar::<usize>()?;
        let nsw     = f.dataset("nsw")?.read_scalar::<usize>()?;
        let efermi  = f.dataset("efermi")?.read_scalar::<f64>()?;
        let dt      = f.dataset("dt")?.read_scalar::<f64>()?;
        let lreal   = f.dataset("lreal")?.read_scalar::<bool>()?;

        let olaps = {
            let olaps_r: Array4<f64> = f.dataset("olaps_r")?.read()?;
            let olaps_i: Array4<f64> = f.dataset("olaps_i")?.read()?;
            olaps_r.mapv(|v| c64::new(v, 0.0)) + olaps_i.mapv(|v| c64::new(0.0, v))
        };

        let eigs: Array3<f64> = f.dataset("eigs")?.read()?;

        let pij = {
            let pij_r: Array5<f64> = f.dataset("pij_r")?.read()?;
            let pij_i: Array5<f64> = f.dataset("pij_i")?.read()?;
            pij_r.mapv(|v| c64::new(v, 0.0)) + pij_i.mapv(|v| c64::new(0.0, v))
        };

        let proj: Array5<f64> = f.dataset("proj")?.read()?;
        let order: Array3<usize> = f.dataset("order")?.read()?;

        Ok(Self {
            ikpoint,
            nspin,
            nbands,
            brange,
            nbrange,
            nsw,
            efermi,
            dt,
            lreal,
            olaps,
            eigs,
            pij,
            proj,
            order,
        })
    }


    pub fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let f = H5File::create(fname)?;

        f.new_dataset::<usize>().create("ikpoint")?.write_scalar(&self.ikpoint)?;
        f.new_dataset::<usize>().create("nspin")?.write_scalar(&self.nspin)?;
        f.new_dataset::<usize>().create("nbands")?.write_scalar(&self.nbands)?;
        f.new_dataset::<[usize;2]>().create("brange")?.write_scalar(&self.brange)?;
        f.new_dataset::<usize>().create("nbrange")?.write_scalar(&self.nbrange)?;
        f.new_dataset::<usize>().create("nsw")?.write_scalar(&self.nsw)?;
        f.new_dataset::<f64>().create("efermi")?.write_scalar(&self.efermi)?;
        f.new_dataset::<f64>().create("dt")?.write_scalar(&self.dt)?;
        f.new_dataset::<bool>().create("lreal")?.write_scalar(&self.lreal)?;

        f.new_dataset_builder().with_data(&self.olaps.mapv(|v| v.re)).create("olaps_r")?;
        f.new_dataset_builder().with_data(&self.olaps.mapv(|v| v.im)).create("olaps_i")?;

        f.new_dataset_builder().with_data(&self.eigs).create("eigs")?;

        f.new_dataset_builder().with_data(&self.pij.mapv(|v| v.re)).create("pij_r")?;
        f.new_dataset_builder().with_data(&self.pij.mapv(|v| v.im)).create("pij_i")?;

        f.new_dataset_builder().with_data(&self.proj).create("proj")?;
        f.new_dataset_builder().with_data(&self.order).create("order")?;

        Ok(())
    }


    pub fn from_inp(inp: &Input) -> Result<Self> {
        if inp.nacfname.is_file() {
            info!("Found pre-calculated NAC available in {:?}, reading NAC from it ...", inp.nacfname);
            let nac = Self::from_h5(&inp.nacfname)?;
            assert_eq!(inp.ikpoint, nac.ikpoint + 1);
            assert_eq!(inp.brange, nac.brange);
            assert_eq!(inp.nsw, nac.nsw);
            assert_eq!(inp.dt, nac.dt);
            assert_eq!(inp.lreal, nac.lreal);
            return Ok(nac);
        }

        info!("No pre-calculated NAC available, start calculating from scratch in {:?}/.../WAVECARs ...", inp.rundir);
        let rundir  = Path::new(&inp.rundir);
        let nsw     = inp.nsw;
        let ikpoint = inp.ikpoint - 1;
        let brange  = Range { start: inp.brange[0] - 1, end: inp.brange[1] };
        let nbrange = brange.len();
        let ndigit  = inp.ndigit;
        let dt      = inp.dt;
        let lreal   = inp.lreal;

        // get nspin and nbands
        let path_1  = rundir.join(format!("{:0ndigit$}", 1)).join("WAVECAR");
        let w1      = Wavecar::from_file(&path_1)
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
            let mut phi = Array3::<c64>::zeros((nspin, nbrange, nplw));

            for ispin in 0 .. nspin {
                for iband in brange.clone().into_iter() {
                    phi.slice_mut(s![ispin, iband - brange.start, ..]).assign(
                        &w1._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                            .into_shape((nplw,))
                            .with_context(|| format!("Wavefunction reshape failed."))?
                    );
                }
            }

            phi
        };

        let procar_1 = rundir.join(format!("{:0ndigit$}", 1)).join("PROCAR");
        let p1      = Procar::from_file(&procar_1)
            .with_context(|| format!("Failed to parse {:?} as PROCAR.", &procar_1))
            .unwrap();
        let nions   = p1.pdos.nions as usize;
        let nproj   = p1.pdos.projected.shape()[4];

        let gvecs   = arr2(&w1.generate_fft_grid_cart(ikpoint as u64))
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
            .collect::<Array1<c64>>()
            .into_shape((nplw, 3))
            .unwrap();

        let lncl = match w1.wavecar_type {
            WavecarType::NonCollinear => true,
            _ => false,
        };

        let CoupTotRet {olaps, eigs, pij, proj, efermi, order} = Self::from_wavecars(
            &phi_1s, &rundir, nsw, ikpoint, brange, ndigit, nspin, lncl, nions, nproj, &gvecs
        )?;

        let ret = Self {
            ikpoint,
            nspin,
            nbands,
            brange: inp.brange,
            nbrange,
            nsw,
            efermi,
            dt,
            lreal,

            olaps,
            eigs,
            pij,

            proj,
            order,
        };

        ret.save_to_h5(&inp.nacfname)?;
        Ok(ret)
    }


    // This function calculates non-adiabatic coupling (NAC), and transition dipole moment (TDM)
    fn from_wavecars(
        phi_1s: &Array3<c64>,
        rundir: &Path, nsw: usize, ikpoint: usize, brange: Range<usize>, ndigit: usize,
        nspin: usize, lncl: bool, nions: usize, nproj: usize, gvecs: &Array2<c64>
        ) -> Result<CoupTotRet>
    {
        let nbrange = brange.clone().count();

        let ret_c_ij = Arc::new(Mutex::new(
                Array4::<c64>::zeros((nsw-1, nspin, nbrange, nbrange))
                ));
        let ret_e_ij = Arc::new(Mutex::new(
                Array3::<f64>::zeros((nsw-1, nspin, nbrange))
                ));
        let ret_p_ij = Arc::new(Mutex::new(
                Array5::<c64>::zeros((nsw-1, nspin, 3, nbrange, nbrange))
                ));

        let nspinors = if lncl { 4 } else { nspin };
        let ret_proj = Arc::new(Mutex::new(
                Array5::<f64>::zeros((nsw-1, nspinors, nbrange, nions, nproj))
                ));
        let efermi_sum = Arc::new(Mutex::new(0.0f64));

        let ret_order = Arc::new(Mutex::new(
                Array3::<usize>::zeros((nsw-1, nspin, nbrange))
                ));

        (0 .. nsw-1).into_par_iter().for_each(|isw| {
            let path_i = rundir.join(format!("{:0ndigit$}", isw + 1));
            let path_j = rundir.join(format!("{:0ndigit$}", isw + 2));

            info!(" Calculating couplings between {:?} and {:?} ...", &path_i, &path_j);

            let CoupIjRet {c_ij, e_ij, p_ij, proj, efermi, order} = Self::coupling_ij(
                phi_1s, path_i.as_path(), path_j.as_path(), ikpoint, brange.clone(), gvecs
                )
                .with_context(|| format!("Failed to calculate couplings between {:?} and {:?}.", &path_i, &path_j))
                .unwrap();

            ret_c_ij.lock().unwrap()
                .slice_mut(s![isw, .., .., ..]).assign(&c_ij);

            ret_e_ij.lock().unwrap()
                .slice_mut(s![isw, .., ..]).assign(&e_ij);

            ret_p_ij.lock().unwrap()
                .slice_mut(s![isw, .., .., .., ..]).assign(&p_ij);

            ret_proj.lock().unwrap()
                .slice_mut(s![isw, .., .., .., ..]).assign(&proj);

            let mut sum = efermi_sum.lock().unwrap();
            *sum += efermi;

            ret_order.lock().unwrap()
                .slice_mut(s![isw, .., ..]).assign(&order);
        });

        Ok( CoupTotRet {
            olaps: Arc::try_unwrap(ret_c_ij).unwrap().into_inner()?,
            eigs: Arc::try_unwrap(ret_e_ij).unwrap().into_inner()?,
            pij: Arc::try_unwrap(ret_p_ij).unwrap().into_inner()?,
            proj: Arc::try_unwrap(ret_proj).unwrap().into_inner()?,
            efermi: Arc::try_unwrap(efermi_sum).unwrap().into_inner()? / (nsw - 1) as f64,
            order: Arc::try_unwrap(ret_order).unwrap().into_inner()?,
        })

    }


    fn coupling_ij(phi_1s: &Array3<c64>, path_i: &Path, path_j: &Path, ikpoint: usize, brange: Range<usize>, gvecs: &Array2<c64>)
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

        let mut c_ij = Array3::<c64>::zeros((nspin, nbrange, nbrange));
        let mut e_ij = Array2::<f64>::zeros((nspin, nbrange));
        let mut p_ij = Array4::<c64>::zeros((nspin, 3, nbrange, nbrange));

        let mut phi_i = Array2::<c64>::zeros((nbrange, nplw));
        let mut phi_j = phi_i.clone();

        let mut phase_i = Array1::<c64>::zeros(nbrange);
        let mut phase_j = phase_i.clone();

        let mut order = Array2::<usize>::zeros((nspin, nbrange));

        for ispin in 0 .. nspin {
            for iband in brange.clone().into_iter() {
                phi_i.slice_mut(s![iband - brange.start, ..]).assign(
                    &wi._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .with_context(|| format!("Wavefunction reshape failed."))?
                );
                phi_j.slice_mut(s![iband - brange.start, ..]).assign(
                    &wj._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .with_context(|| format!("Wavefunction reshape failed."))?
                );
            }

            // TODO: impl reordering
            //
            // Reorder band indices to solve band crossing issues.
            //
            // This is a linear assignment problem, use kuhn munkres or hungarian algorithm to
            // solve it
            //     use pathfinding::kuhn_munkres::kuhn_munkres (or hungarian algorithm)

            // Reordering j, u1i = <phi_1 | phi_j>, the order of first step must be ordered,
            // i.e. (0, 1, ...), so we only store the order of j-th step.
            let order_j: Array1<usize> = Self::find_order(&phi_1s.slice(s![ispin, .., ..]), &phi_j).into();
            order.slice_mut(s![ispin, ..]).assign(&order_j);

            //info!("path = {:?}, order_i = {:?}", path_i, order_i);
            //info!("path = {:?}, order_j = {:?}", path_j, order_j);

            // phase correction:
            //             < phi_0 | phi_i >
            // phase_0i = -------------------
            //            |< phi_0 | phi_i >|
            // phi_i[:] *= conj(phase_0i)

            phase_i.assign( &(phi_1s.slice(s![ispin, .., ..]).mapv(|x| x.conj()) * &phi_i).sum_axis(Axis(1)) );
            phase_i.mapv_inplace(|x| x.conj() / x.norm());
            phi_i.axis_iter_mut(Axis(0)).zip(phase_i.iter())
                .par_bridge()
                .for_each(|(mut row, phase)| row.mapv_inplace(|x| x * phase));

            phase_j.assign( &(phi_1s.slice(s![ispin, .., ..]).mapv(|x| x.conj()) * &phi_j).sum_axis(Axis(1)) );
            phase_j.mapv_inplace(|x| x.conj() / x.norm());
            phi_j.axis_iter_mut(Axis(0)).zip(phase_j.iter())
                .par_bridge()
                .for_each(|(mut row, phase)| row.mapv_inplace(|x| x * phase));


            let phi_ij = phi_i.mapv(|v| v.conj()).dot(&phi_j.t());
            let phi_ji = phi_ij.t().mapv(|v| v.conj());  // phi_ji = phi_ij^H
            c_ij.slice_mut(s![ispin, .., ..]).assign(&(phi_ij - phi_ji));

            for idirect in 0 .. 3 {
                let phi_x_gvecs: Array2<_> = phi_j.clone() * gvecs.slice(s![NewAxis, .., idirect]);

                // <i | p | j>, in eV*fs/Angstrom
                let p_ij_tmp = match wi.wavecar_type {
                    WavecarType::GammaHalf(_) => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t())
                                               - phi_x_gvecs.mapv(|v| v.conj()).dot(&phi_j.t()),
                    _ => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t()),
                };
                p_ij.slice_mut(s![ispin, idirect, .., ..]).assign(&p_ij_tmp);
            }
        }

        // read PROCAR
        let proj = {
            let fpath = path_i.join("PROCAR");
            let proj_i = Procar::from_file(&fpath)
                .with_context(|| format!("Failed to parse {:?}.", &fpath))
                .unwrap();
            proj_i.pdos.projected.slice(s![.., ikpoint, brange.clone(), .., ..]).to_owned()
        };

        e_ij.slice_mut(s![.., ..]).assign(&(
            ( wi.band_eigs.slice(s![.., ikpoint, brange.clone()]).to_owned() + 
              wj.band_eigs.slice(s![.., ikpoint, brange.clone()]) ) / 2.0
        ));

        Ok( CoupIjRet{
            c_ij,
            e_ij, 
            p_ij, 
            proj,
            efermi: wi.efermi,
            order
        })
    }


    fn find_order<R1, R2>(phi_i: &ArrayBase<R1, Ix2>, phi_j: &ArrayBase<R2, Ix2>) -> Vec<usize>
    where
        R1: Data<Elem = c64>,
        R2: Data<Elem = c64>,
    {

        let uij = phi_i.mapv(|x| x.conj())
            .dot(&phi_j.t())
            .mapv(|x| OrderedFloat(x.norm_sqr()));

        let weights = pfMatrix::square_from_vec(uij.into_raw_vec()).unwrap();
        let (_maxcoup, order) = kuhn_munkres(&weights);

        return order;
    }
}



#[cfg(test)]
mod tests {
    use std::time::Instant;
    use super::*;
    use tempfile::tempdir;

    #[test]
    #[ignore]
    fn test_from_wavecars() {
        let ndigit  = 5;
        let rundir  = Path::new("/data2/chenlj/Work/Work2022/NAMD_lumi/2022-11-14_GaAs/2x2x2/aimd/static_0.1fs/run_3000fs");
        let path_1  = rundir.join(format!("{:0ndigit$}", 1)).join("WAVECAR");
        let w1      = Wavecar::from_file(&path_1).unwrap();
        let nsw     = 300;
        let ikpoint = 0;
        let brange  = 0 .. 150;
        let nbrange = brange.len();
        let nplw    = w1.nplws[ikpoint] as usize;
        let nions   = 4usize;
        let nproj   = 9usize;
        let gvecs   = arr2(&w1.generate_fft_grid_cart(ikpoint as u64)).mapv(|v| c64::new(v, 0.0));
        let phi_1s  = Array3::<c64>::zeros((1, nbrange, nplw));

        let now = Instant::now();
        let (_c_ij, _e_ij, _p_ij, _proj, efermi) = Nac::from_wavecars(&phi_1s, &rundir, nsw, ikpoint, brange, ndigit,
            1, false, nions, nproj, &gvecs).unwrap();
        println!("efermi = {}, time used: {:?}", efermi, now.elapsed());
    }

    #[test]
    fn test_h5() {
        let ikpoint = 0usize;
        let nspin   = 1usize;
        let nbands  = 114usize;
        let brange  = [10usize, 30];
        let nbrange = 20usize;
        let nsw     = 3000usize;
        let efermi  = 1.14514f64;
        let dt      = 1.0f64;
        let lreal   = true;
        let proj    = Array5::zeros((nsw-1, nspin, nbands, 3, 3));

        let nac = Nac {
            ikpoint,
            nspin,
            nbands,
            brange,
            nbrange,
            nsw,
            efermi,
            dt,
            lreal,

            olaps: Array4::<c64>::zeros((nsw-1, nspin, nbrange, nbrange)),
            eigs: Array3::<f64>::zeros((nsw-1, nspin, nbrange)),
            pij: Array5::<c64>::zeros((nsw-1, nspin, 3, nbrange, nbrange)),

            proj,
        };

        let dir = tempdir().unwrap();
        let fname = dir.path().join("test_nac.h5");
        nac.save_to_h5(&fname).unwrap();
        let nac_read = Nac::from_h5(&fname).unwrap();
        assert_eq!(nac, nac_read);
    }
}
