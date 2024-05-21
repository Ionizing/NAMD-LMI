use std::path::Path;
//use std::fmt;

use pathfinding::prelude::{
    kuhn_munkres,
    Matrix as pfMatrix,
};
use ordered_float::OrderedFloat;
use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::{
    info,
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
use crate::core::constants::*;


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
    reorder:     bool,
    scissor:     Option<f64>,

    eig_t:     nd::Array2<f64>,     // [nsw-1, nbasis]
    nac_t:     nd::Array3<c64>,     // [nsw-1, nbasis, nbasis]
    pij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    rij_t:     nd::Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    proj_t:    nd::Array4<f64>,     // [nsw-1, nbasis, nions, nproj]
    efield:    Option<Efield<'a>>,

    // pre-calculated hamiltonian = eig_t in diag + nac_t in off-diag
    hamil0:    nd::Array3<c64>,     // [nsw-1, nbasis, nbasis]
}


impl<'a> Hamiltonian for SPHamiltonian<'a> {
    type ConfigType   = HamilConfig;
    type CouplingType = Nac;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nsw(&self) -> usize { self.nsw }
    fn get_temperature(&self) -> f64 { self.temperature }

    /// Get the hamiltonian without electron-photon interaction.
    fn get_hamil(&self, iion: usize) -> nd::ArrayView2<c64> {
        self.hamil0.slice(nd::s![iion, .., ..])
    }

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
        let reorder  = f.dataset("reorder")?.read_scalar::<bool>()?;
        let scissor  = f.dataset("scissor")?.read_scalar::<f64>().ok();

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

        let hamil0 = Self::calculate_hamil0(&eig_t, &nac_t);

        Ok(SPHamiltonian {
            ikpoint,
            basis_up,
            basis_dn,
            nbasis,
            potim,
            nsw,
            temperature,
            propmethod,
            reorder,
            scissor,

            eig_t,
            nac_t,
            pij_t,
            rij_t,
            proj_t,
            efield,

            hamil0,
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
        f.new_dataset::<bool>().create("reorder")?.write_scalar(&self.reorder)?;
        if let Some(scissor) = self.scissor.as_ref() {
            f.new_dataset::<f64>().create("scissor")?.write_scalar(scissor)?;
        }

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
        let reorder              = cfg.get_reorder();

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

        if let Some(scissor) = cfg.get_scissor() {
            apply_scissor(&mut eig_t, scissor);
        }

        let hamil0 = Self::calculate_hamil0(&eig_t, &nac_t);

        Ok(Self {
            ikpoint,
            basis_up,
            basis_dn,
            nbasis,
            potim,
            nsw,
            temperature,
            propmethod,
            reorder,
            scissor: cfg.get_scissor(),

            eig_t,
            nac_t,
            pij_t,
            rij_t,
            proj_t,
            efield,

            hamil0,
        })
    }


    pub fn get_rtime_xtime(iion: usize, nsw: usize, namdinit: usize) -> [usize; 2] {
        let rtime = (iion + namdinit) % (nsw - 2);
        let xtime = rtime + 1;
        [rtime, xtime]
    }


    pub fn get_hamil0_rtime(&self, iion: usize, namdinit: usize) -> nd::ArrayView2<c64> {
        let [rtime, _] = Self::get_rtime_xtime(iion, self.nsw, namdinit);
        self.get_hamil(rtime)
    }


    pub fn get_pij_t(&self) -> nd::ArrayView4<c64> {
        self.pij_t.view()
    }


    pub fn get_pij_rtime(&self, iion: usize, namdinit: usize) -> nd::ArrayView3<c64> {
        let [rtime, _] = Self::get_rtime_xtime(iion, self.nsw, namdinit);
        self.pij_t.slice(nd::s![rtime, .., .., ..])
    }


    pub fn get_rij_t(&self) -> nd::ArrayView4<c64> {
        self.rij_t.view()
    }


    pub fn get_rij_rtime(&self, iion: usize, namdinit: usize) -> nd::ArrayView3<c64> {
        let [rtime, _] = Self::get_rtime_xtime(iion, self.nsw, namdinit);
        self.rij_t.slice(nd::s![rtime, .., .., ..])
    }


    pub fn get_efield(&self) -> Option<&Efield> {
        self.efield.as_ref()
    }


    pub fn get_propmethod(&self) -> PropagateMethod {
        self.propmethod
    }


    // H_diag = eig_t
    // H_offdiag = nac_t * -i hbar
    fn calculate_hamil0(eig_t: &nd::Array2<f64>, nac_t: &nd::Array3<c64>) -> nd::Array3<c64> {
        // off-diag = -i * \hbar * NAC
        let mut ret = -IMGUNIT * HBAR * nac_t;
        let nsw = nac_t.shape()[0];

        // diag = eig
        for i in 0 .. nsw {
            ret.slice_mut(nd::s![i, .., ..])
                .diag_mut()
                .assign(&eig_t.slice(nd::s![i, ..]).mapv(|v| c64::new(v, 0.0)));
        }

        ret
    }


    fn find_order(cij: nd::ArrayView2<c64>) -> nd::Array1<usize> {
        let uij = cij.mapv(|v| OrderedFloat(v.norm_sqr()));
        let weights = pfMatrix::square_from_vec(uij.into_raw_vec()).unwrap();
        let (_maxcoup, order) = kuhn_munkres(&weights);
        return nd::Array1::<usize>::from(order);
    }


    fn find_all_orders(cij: nd::ArrayView3<c64>) -> nd::Array2<usize> {
        let (nsw, nbasis, _) = cij.dim();
        let mut orders = nd::Array2::<usize>::zeros((nsw, nbasis));
        for isw in 0 .. nsw {
            orders.slice_mut(nd::s![isw, ..])
                .assign(&Self::find_order(cij.slice(nd::s![isw, .., ..])));
        }

        for isw in 1 .. nsw {
            let perm_ij = orders.slice(nd::s![isw, ..]).to_owned();
            for iband in 0 .. nbasis {
                orders[(isw, iband)] = orders[(isw-1, perm_ij[iband])];
            }
        }

        return orders;
    }


    fn apply_reorder(
        eig_t: &mut nd::ArrayViewMut2<f64>,
        nac_t: &mut nd::ArrayViewMut3<c64>,
        pij_t: &mut nd::ArrayViewMut4<c64>,
        rij_t: &mut nd::ArrayViewMut4<c64>,
        proj_t: &mut nd::ArrayViewMut4<f64>,
        ) {
        let (nsw, nbasis) = eig_t.dim();
        assert!((nsw, nbasis, nbasis) == nac_t.dim());
        assert!((nsw, 3, nbasis, nbasis) == pij_t.dim());
        assert!((nsw, 3, nbasis, nbasis) == rij_t.dim());
        assert!(nsw == proj_t.shape()[0] && nbasis == proj_t.shape()[1]);

        let orders = Self::find_all_orders(nac_t.view());
        let (nsw, nbasis) = orders.dim();

        let orders_sorted = nd::Array1::<usize>::from_iter(0 .. nbasis);

        for isw in 0 .. nsw {
            let iorder = orders.slice(nd::s![isw, ..]);

            // sort the eigen values and projectors
            let proj_t_org = proj_t.slice(nd::s![isw, .., .., ..]).to_owned();
            for iband in 0 .. nbasis {
                eig_t[(isw, iorder[iband])] = eig_t[(isw, iband)];
                proj_t.slice_mut(nd::s![isw, iorder[iband], .., ..])
                    .assign(&proj_t_org.slice(nd::s![iband, .., ..]));
            }

            // sort NAC = <i| d/dt |j> and momentum matrix <i| p |j> and dipole matrix <i| r |j>
            //  prepare the indices
            let iorder = if 0 == isw { orders_sorted.view() } else { orders.slice(nd::s![isw-1, ..]) };
            let jorder = orders.slice(nd::s![isw, ..]);

            let nac_t_org = nac_t.slice(nd::s![isw, .., ..]).to_owned();
            let pij_t_org = pij_t.slice(nd::s![isw, .., .., ..]).to_owned();
            let rij_t_org = rij_t.slice(nd::s![isw, .., .., ..]).to_owned();
            for iband in 0 .. nbasis {
                for jband in 0 .. nbasis {
                    nac_t[(isw, iorder[iband], jorder[jband])] = nac_t_org[(iband, jband)];
                    pij_t[(isw, 0, iorder[iband], jorder[jband])] = pij_t_org[(0, iband, jband)];
                    pij_t[(isw, 1, iorder[iband], jorder[jband])] = pij_t_org[(1, iband, jband)];
                    pij_t[(isw, 2, iorder[iband], jorder[jband])] = pij_t_org[(2, iband, jband)];
                    rij_t[(isw, 0, iorder[iband], jorder[jband])] = rij_t_org[(0, iband, jband)];
                    rij_t[(isw, 1, iorder[iband], jorder[jband])] = rij_t_org[(1, iband, jband)];
                    rij_t[(isw, 2, iorder[iband], jorder[jband])] = rij_t_org[(2, iband, jband)];

                    // alternative form, may be much slower
                    //pij_t.slice_mut(nd::s![isw, .., iorder[iband], jorder[jband]])
                        //.assign(&pij_t_org.slice(nd::s![.., iband, jband]));
                    //rij_t.slice_mut(nd::s![isw, .., iorder[iband], jorder[jband]])
                        //.assign(&rij_t_org.slice(nd::s![.., iband, jband]));
                }
            }
        }
        return;
    }
}

fn idx_convert(brange: &[usize; 2], idx: usize, shift: usize) -> Result<usize> {
    ensure!(idx >= brange[0] && idx <= brange[1], "Index out of bounds.");
    Ok(idx - brange[0] + shift)
}


fn apply_scissor(eigs: &mut nd::Array2<f64>, scissor: f64) {
    let eigs_avg = eigs.mean_axis(nd::Axis(0)).unwrap().to_vec();
    let cbm = eigs_avg.partition_point(|&x| x < 0.0);
    let vbm = cbm - 1;
    let gap = eigs_avg[cbm] - eigs_avg[vbm];

    let gaps = eigs.slice(nd::s![.., cbm]).to_owned() - eigs.slice(nd::s![.., vbm]);
    let gap_min = gaps.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let gap_max = gaps.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();

    let shift = scissor - gap;
    let newgap_min = gap_min + shift;
    let newgap_max = gap_max + shift;

    // Without the asterisk, `cbm_min` would be fucking &f64 instead of f64
    let cbm_min: f64 = *eigs.slice(nd::s![.., cbm]).iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let vbm_max: f64 = *eigs.slice(nd::s![.., vbm]).iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let gap_t_min    = cbm_min - vbm_max + shift;

    eigs.slice_mut(nd::s![.., cbm ..]).mapv_inplace(|x| x + shift);
    info!("Found scissor opeartor of {scissor:.4} eV. current system has gap of \
           {gap_min:.4} .. {gap:.4} .. {gap_max:.4} (min .. avg .. max) (eV).
                Now the gap is set to {newgap_min:.4} .. {scissor:.4} .. {newgap_max:.4} \
                (min .. avg .. max) (eV). min(CBM_t) - max(VBM_t) = {:8.4}", gap_t_min);
}
