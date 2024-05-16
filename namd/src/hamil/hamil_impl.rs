use std::path::Path;
use std::fmt;

use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::{c64, Result};

use crate::core::Hamiltonian;
use crate::nac::nac_impl::Nac;
use crate::hamil::config::HamilConfig;
use crate::hamil::efield::Efield;


/// Sing-Particle Hamiltonian
pub struct SPHamiltonian<'a> {
    ikpoint:     usize,
    basis_up:    [usize; 2],
    basis_dn:    [usize; 2],
    nbasis:      usize,
    potim:       f64,
    nsw:         usize,
    temperature: f64,

    eig_t:     nd::Array2<f64>,     // [nsw-1, nbasis]
    nac_t:     nd::Array3<c64>,     // [nsw-1, nbasis, nbasis]
    pij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    rij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    proj_t:    nd::Array4<f64>,     // [nsw-1, nbasis, nions, nproj]
    efield:    Option<Efield<'a>>,
}


impl<'a> Hamiltonian<'a> for SPHamiltonian<'a> {
    type ConfigType   = HamilConfig;
    type CouplingType = Nac;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nsw(&self) -> usize { self.nsw }
    fn get_temperature(&self) -> f64 { self.temperature }

    fn get_hamil(&self, iion: usize) -> nd::ArrayView2<c64> {
        todo!()
    } // [nbasis, nbasis]

    fn from_config(fname: &Self::ConfigType, coup: &Self::CouplingType) -> Result<Self> {
        todo!()
    }

    fn from_h5<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        let f = H5File::open(fname)?;

        let ikpoint  = f.dataset("ikpoint")?.read_scalar::<usize>()?;
        let basis_up = f.dataset("basis_up")?.read_scalar::<[usize;2]>()?;
        let basis_dn = f.dataset("basis_dn")?.read_scalar::<[usize;2]>()?;
        let nbasis   = f.dataset("nbasis")?.read_scalar::<usize>()?;
        let potim    = f.dataset("potim")?.read_scalar::<f64>()?;
        let nsw      = f.dataset("nsw")?.read_scalar::<usize>()?;
        let temperature = f.dataset("temperature")?.read_scalar::<f64>()?;

        let eig_t: nd::Array2<f64> = f.dataset("eig_t")?.read()?;
        let nac_t = {
            let nac_t_r: nd::Array3<f64> = f.dataset("nac_t_r")?.read()?;
            let nac_t_i: nd::Array3<f64> = f.dataset("nac_t_i")?.read()?;
            nac_t_r.mapv(|v| c64::new(v, 0.0)) + nac_t_i.mapv(|v| c64::new(0.0, v))
        };
        let pij_t = {
            let pij_t_r: nd::Array4<f64> = f.dataset("pij_t_r")?.read()?;
            let pij_t_i: nd::Array4<f64> = f.dataset("pij_t_i")?.read()?;
            pij_t_r.mapv(|v| c64::new(v, 0.0)) + pij_t_i.mapv(|v| c64::new(0.0, v))
        };
        let rij_t = {
            let rij_t_r: nd::Array4<f64> = f.dataset("rij_t_r")?.read()?;
            let rij_t_i: nd::Array4<f64> = f.dataset("rij_t_i")?.read()?;
            rij_t_r.mapv(|v| c64::new(v, 0.0)) + rij_t_i.mapv(|v| c64::new(0.0, v))
        };
        let proj_t: nd::Array4<f64> = f.dataset("proj_t")?.read()?;

        let efield = {
            if f.dataset("efield").is_err() {
                None
            } else {
                let raw: Vec<u8> = f.dataset("efield")?.read_raw()?;
                let efield_src = String::from_utf8(raw)?;
                Some(Efield::from_str(&efield_src)?)
            }
        };

        Ok(SPHamiltonian {
            ikpoint,
            basis_up,
            basis_dn,
            nbasis,
            potim,
            nsw,
            temperature,

            eig_t,
            nac_t,
            pij_t,
            rij_t,
            proj_t,
            efield
        })
    }

    fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        let f = H5File::create(fname)?;

        f.new_dataset::<usize>().create("ikpoint")?.write_scalar(&self.ikpoint)?;
        f.new_dataset::<[usize;2]>().create("basis_up")?.write_scalar(&self.basis_up)?;
        f.new_dataset::<[usize;2]>().create("basis_dn")?.write_scalar(&self.basis_dn)?;
        f.new_dataset::<usize>().create("nbasis")?.write_scalar(&self.nbasis)?;
        f.new_dataset::<f64>().create("potim")?.write_scalar(&self.potim)?;
        f.new_dataset::<usize>().create("nsw")?.write_scalar(&self.nsw)?;
        f.new_dataset::<f64>().create("temperature")?.write_scalar(&self.temperature)?;

        f.new_dataset_builder().with_data(&self.eig_t).create("eig_t")?;

        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.re)).create("nac_t_r")?;
        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.im)).create("nac_t_i")?;

        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.re)).create("pij_t_r")?;
        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_i")?;

        f.new_dataset_builder().with_data(&self.rij_t.mapv(|v| v.re)).create("rij_t_r")?;
        f.new_dataset_builder().with_data(&self.rij_t.mapv(|v| v.im)).create("rij_t_i")?;

        f.new_dataset_builder().with_data(&self.proj_t).create("proj_t")?;

        if let Some(efield) = self.efield.as_ref() {
            f.new_dataset_builder().with_data(efield.get_src()).create("efield");
        }

        Ok(())
    }
} 
