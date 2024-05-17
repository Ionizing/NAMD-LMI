use std::path::Path;
//use std::fmt;

use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::{
    c64,
    Result,
    anyhow::ensure,
};

use crate::core::{
    Couplings,
    Hamiltonian,
};
use crate::nac::nac_impl::Nac;
use crate::hamil::config::HamilConfig;
use crate::hamil::config::PropagateMethod;
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
    propmethod:  PropagateMethod,

    eig_t:     nd::Array2<f64>,     // [nsw-1, nbasis]
    nac_t:     nd::Array3<c64>,     // [nsw-1, nbasis, nbasis]
    pij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    rij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    proj_t:    nd::Array4<f64>,     // [nsw-1, nbasis, nions, nproj]
    efield:    Option<Efield<'a>>,
}


impl<'a> Hamiltonian for SPHamiltonian<'a> {
    type ConfigType   = HamilConfig;
    type CouplingType = Nac;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nsw(&self) -> usize { self.nsw }
    fn get_temperature(&self) -> f64 { self.temperature }

    fn get_hamil(&self, iion: usize) -> nd::ArrayView2<c64> {
        todo!()
    } // [nbasis, nbasis]

    fn from_config(cfg: &Self::ConfigType) -> Result<Self> {
        let coup = Nac::from_h5(cfg.get_nac_fname())?;
        Self::with_config_and_coupling(cfg, &coup)
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
        let propmethod = {
            let raw: Vec<u8> = f.dataset("propmethod")?.read_raw()?;
            let src = String::from_utf8(raw)?;
            PropagateMethod::from_str(&src)?
        };

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
            propmethod,

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

        f.new_dataset_builder().with_data(&self.propmethod.to_string()).create("propmethod")?;
        f.new_dataset_builder().with_data(&self.eig_t).create("eig_t")?;

        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.re)).create("nac_t_r")?;
        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.im)).create("nac_t_i")?;

        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.re)).create("pij_t_r")?;
        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_i")?;

        f.new_dataset_builder().with_data(&self.rij_t.mapv(|v| v.re)).create("rij_t_r")?;
        f.new_dataset_builder().with_data(&self.rij_t.mapv(|v| v.im)).create("rij_t_i")?;

        f.new_dataset_builder().with_data(&self.proj_t).create("proj_t")?;

        if let Some(efield) = self.efield.as_ref() {
            f.new_dataset_builder().with_data(efield.get_src()).create("efield")?;
        }

        Ok(())
    }
} 


impl<'a> SPHamiltonian<'a> {
    fn with_config_and_coupling(cfg: &HamilConfig, coup: &Nac) -> Result<Self> {
        cfg.check_config()?;
        ensure!(cfg.get_ikpoint() == coup.get_ikpoint());

        let ikpoint: usize       = cfg.get_ikpoint();
        let basis_up: [usize; 2] = cfg.get_basis_up();
        let basis_dn: [usize; 2] = cfg.get_basis_dn();
        let potim: f64           = coup.get_potim();
        let nsw: usize           = coup.get_nsw();
        let temperature: f64     = coup.get_temperature();
        let propmethod           = cfg.get_propmethod();
        let brange               = coup.get_brange();

        let mut nb = [0usize; 2];
        nb[0] = if basis_up.contains(&0) {
            0
        } else {
            basis_up[1] - basis_up[0] + 1
        };
        nb[1] = if basis_dn.contains(&0) {
            0
        } else {
            basis_dn[1] - basis_dn[0] + 1
        };

        let nbasis = nb[0] + nb[1];
        let bup_nac = if nb[0] == 0 {
            0 .. 0
        } else {
            idx_convert(&brange, basis_up[0], 0)?
                ..
            idx_convert(&brange, basis_up[1], 1)?
        };
        let bup_hamil = 0 .. nb[0];

        let bdn_nac = if nb[1] == 0 {
            0 .. 0
        } else {
            idx_convert(&brange, basis_dn[0], 0)?
                ..
            idx_convert(&brange, basis_dn[1], 1)?
        };
        let bdn_hamil = nb[0] .. nbasis;

        let nions = coup.get_tdproj().shape()[3];
        let nproj = coup.get_tdproj().shape()[4];

        let mut eig_t = nd::Array2::<f64>::zeros((nsw-1, nbasis));
        let mut nac_t = nd::Array3::<c64>::zeros((nsw-1, nbasis, nbasis));
        let mut pij_t = nd::Array4::<c64>::zeros((nsw-1, 3, nbasis, nbasis));
        let mut rij_t = nd::Array4::<c64>::zeros((nsw-1, 3, nbasis, nbasis));
        let mut proj_t = nd::Array4::<f64>::zeros((nsw-1, nbasis, nions, nproj));

        // Damn, Range doesn't impl Copy trait. Fxxk up.
        eig_t.slice_mut(nd::s![.., bup_hamil.clone()])
            .assign(&coup.get_tdeigs().slice(nd::s![.., 0, bup_nac.clone()]));
        eig_t.slice_mut(nd::s![.., bdn_hamil.clone()])
            .assign(&coup.get_tdeigs().slice(nd::s![.., 1, bdn_nac.clone()]));

        nac_t.slice_mut(nd::s![.., bup_hamil.clone(), bup_hamil.clone()])
            .assign(&coup.get_tdcoup().slice(nd::s![.., 0, bup_nac.clone(), bup_nac.clone()]));
        nac_t.slice_mut(nd::s![.., bdn_hamil.clone(), bdn_hamil.clone()])
            .assign(&coup.get_tdcoup().slice(nd::s![.., 1, bdn_nac.clone(), bdn_nac.clone()]));

        pij_t.slice_mut(nd::s![.., .., bup_hamil.clone(), bup_hamil.clone()])
            .assign(&coup.get_tdpij().slice(nd::s![.., 0, .., bup_nac.clone(), bup_nac.clone()]));
        pij_t.slice_mut(nd::s![.., .., bdn_hamil.clone(), bdn_hamil.clone()])
            .assign(&coup.get_tdpij().slice(nd::s![.., 1, .., bdn_nac.clone(), bdn_nac.clone()]));

        rij_t.slice_mut(nd::s![.., .., bup_hamil.clone(), bup_hamil.clone()])
            .assign(&coup.get_tdrij().slice(nd::s![.., 0, .., bup_nac.clone(), bup_nac.clone()]));
        rij_t.slice_mut(nd::s![.., .., bdn_hamil.clone(), bdn_hamil.clone()])
            .assign(&coup.get_tdrij().slice(nd::s![.., 1, .., bdn_nac.clone(), bdn_nac.clone()]));

        proj_t.slice_mut(nd::s![.., bup_hamil.clone(), .., ..])
            .assign(&coup.get_tdproj().slice(nd::s![.., 0, bup_hamil.clone(), .., ..]));
        proj_t.slice_mut(nd::s![.., bdn_hamil.clone(), .., ..])
            .assign(&coup.get_tdproj().slice(nd::s![.., 1, bdn_hamil.clone(), .., ..]));

        let efield: Option<Efield> = cfg.get_efield_fname()
            .map(|fname| Efield::from_file(fname).unwrap());

        Ok(Self {
            ikpoint,
            basis_up,
            basis_dn,
            nbasis,
            potim,
            nsw,
            temperature,
            propmethod,

            eig_t,
            nac_t,
            pij_t,
            rij_t,
            proj_t,
            efield
        })
    }


}

fn idx_convert(brange: &[usize; 2], idx: usize, shift: usize) -> Result<usize> {
    ensure!(idx >= brange[0] && idx <= brange[1], "Index out of bounds.");
    Ok(idx - brange[0] + shift)
}
