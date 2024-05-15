use std::path::Path;
use std::fmt;

use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::{c64, Result};

use crate::core::Hamiltonian;
use crate::nac::nac_impl::Nac;
use crate::hamil::config::{
    HamilConfig,
    PropagateMethod,
};
use crate::hamil::efield::Efield;


/// Sing-Particle Hamiltonian
pub struct SPHamiltonian {
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
    efield:    Option<Efield>,

    hamil:     nd::Array2<c64>,     // [nbasis, nbasis]
}


impl<'a> Hamiltonian<'a> for SPHamiltonian {
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
        todo!()
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
        }
        
        todo!()
    }
} 
