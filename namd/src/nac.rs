use std::path::Path;
use std::ops::RangeInclusive;

use mpi::traits::*;

use shared::{
    Result,
    c64,
    ndarray::{
        s,
        Array2,
        Array3,
        Array4,
        Array5,
        NewAxis,
    },
};
use vasp_parsers::{
    Wavecar,
    WavecarType,
};

pub struct Nac {
    pub ikpoint: usize,
    pub nspin:   usize,
    pub nbands:  usize,
    pub brange:  [usize; 2],
    pub nbrange: usize,
    pub nsw:     usize,
    pub efermi:  f64,
    pub dt:      f64,
    pub lreal:   bool,

    pub olaps:   Array4<c64>,
    pub eigs:    Array3<f64>,
    pub pij:     Array5<c64>,
}

impl Nac {
    pub fn from_h5<P>(_fname: &P) -> Self
    where
        P: AsRef<Path> + ?Sized,
    {
        todo!()
    }


    pub fn from_inp() -> Self {
        todo!()
    }


    fn from_wavecars(&mut self) -> Result<()> {
        todo!()
    }


    fn nac_ij(wi: &Wavecar, wj: &Wavecar, ikpoint: usize, brange: RangeInclusive<usize>, gvecs: &Array2<c64>)
        -> Result<(Array3<c64>, Array2<f64>, Array4<c64>)>
    {
        let nplw     = wi.nplws[ikpoint] as usize;
        let nbrange  = brange.clone().count();
        let nspin    = wi.nspin as usize;
        let nspinor  = match wi.wavecar_type {
            WavecarType::NonCollinear => 2,
            _ => 1usize,
        };

        assert_eq!(nplw, gvecs.shape()[0]);

        let mut c_ij = Array3::<c64>::zeros((nspin, nbrange, nbrange));
        let mut e_ij = Array2::<f64>::zeros((nspin, nbrange));
        let mut p_ij = Array4::<c64>::zeros((nspin, 3, nbrange, nbrange));

        let mut phi_i = Array2::<c64>::zeros((nbrange, nplw));
        let mut phi_j = phi_i.clone();

        for ispin in 0 .. (wi.nspin as usize) {
            for iband in brange.clone().into_iter() {
                phi_i.slice_mut(s![iband, ..]).assign(
                    &wi._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw)
                );
                phi_j.slice_mut(s![iband, ..]).assign(
                    &wj._wav_kspace(ispin as u64, ikpoint as u64, iband as u64, nplw)
                        .into_shape((nspinor * nplw,)).unwrap()
                );
            }

            let phi_ij = phi_i.mapv(|v| v.conj()).dot(&phi_j.t());
            let phi_ji = phi_i.t().mapv(|v| v.conj());  // phi_ji = phi_ij^H
            c_ij.slice_mut(s![ispin, .., ..]).assign(&(phi_ij - phi_ji));

            for idirect in 0 .. 3 {
                let phi_x_gvecs: Array2<_> = phi_j.clone() * gvecs.slice(s![NewAxis, .., idirect]);

                // <phi_i | p | phi_j>
                let p_ij_tmp = match wi.wavecar_type {
                    WavecarType::GammaHalf(_) => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t())
                                               - phi_x_gvecs.mapv(|v| v.conj()).dot(&phi_j.t()),
                    _ => phi_j.mapv(|v| v.conj()).dot(&phi_x_gvecs.t()),
                };
                p_ij.slice_mut(s![ispin, idirect, .., ..]).assign(&p_ij_tmp);
            }

            e_ij.slice_mut(s![ispin, ..]).assign(&(
                ( wi.band_eigs.slice(s![ispin, ikpoint, brange.clone()]).to_owned() + 
                  wj.band_eigs.slice(s![ispin, ikpoint, brange.clone()]) ) / 2.0
            ));
        }

        Ok((c_ij, e_ij, p_ij))
    }
}
