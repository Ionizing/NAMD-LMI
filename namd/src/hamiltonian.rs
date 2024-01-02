use std::path::Path;
use std::fmt;

use shared::{
    Result,
    ndarray::s,
    c64,
    Array1,
    Array2,
    Array3,
    ndarray::Array4,
    ndarray_linalg::{
        EighInplace,
        UPLO,
        expm,
    },
    tracing::{self, instrument, info},
};
use hdf5::File as H5File;

use crate::{
    constants::*,
    efield::Efield,
    nac::Nac,
    input::Input,
};


#[derive(PartialEq, Eq, Clone, Copy)]
pub enum PropagateMethod {
    FiniteDifference,
    Exact,
    Expm,
    LiouvilleTrotter,
}


impl fmt::Display for PropagateMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use PropagateMethod::*;
        write!(f, "{}", match self {
            FiniteDifference => "FiniteDifference",
            Exact            => "Exact",
            Expm             => "Expm",
            LiouvilleTrotter => "LiouvilleTrotter",
        })
    }
}


pub struct Hamiltonian {
    pub basis_up:      [usize; 2],
    pub basis_dn:      [usize; 2],
    pub nbasis:        usize,
    pub dt:            f64,
    pub edt:           f64,
    pub basisini:      usize,
    pub namdinit:      usize,
    pub namdtime:      usize,
    pub nsw:           usize,
    pub nelm:          usize,
    pub lreal:         bool,
    pub temperature:   f64,
    //pub efield_lcycle: bool,

    pub psi_p:         Array1<c64>,     // [nbasis]
    pub psi_c:         Array1<c64>,     // [nbasis]
    pub psi_n:         Array1<c64>,     // [nbasis]
    pub psi_t:         Array2<c64>,     // [namdtime, nbasis]
    pub pop_t:         Array2<f64>,     // [nbasis]
    pub psi_h:         Array1<c64>,     // [nbasis]
    
    pub hamil:         Array2<c64>,     // [nbasis, nbasis]
    pub eig_t:         Array2<f64>,     // [nsw-1, nbasis]
    pub prop_eigs:     Array1<f64>,     // [namdtime]
    pub nac_t:         Array3<c64>,     // [nsw-1, nbasis, nbasis]
    pub pij_t:         Array4<c64>,     // [nsw-1, 3, nbasis, nbasis]
    pub efield:        Option<Efield>,
    pub proj:          Array4<f64>,     // [nsw-1, nbasis, nions, nproj]

    // Auxiliary variables used by `make_hamil method`
    delta_eig:         Array1<f64>,     // [nbasis]
    delta_nac:         Array2<c64>,     // [nbasis, nbasis]
    delta_pij:         Array3<c64>,     // [3, nbasis, nbasis]
    pub lmi_t:         Array3<c64>,     // [namdtime, nbasis, nbasis]
    pub ham_t:         Array3<c64>,     // [namdtime, nbasis, nbasis]
}


impl Hamiltonian {
    pub fn from_input(nac: &Nac, inp: &Input, inicon_idx: usize) -> Self {
        Self::init_with_nac(
            nac,
            inp.efield.clone().map(|x| x.2),
            inp.basis_up.clone(),
            inp.basis_dn.clone(),
            inp.dt,
            inp.inisteps[inicon_idx],
            inp.inibands[inicon_idx],
            inp.inispins[inicon_idx] - 1,
            inp.namdtime,
            inp.nelm,
            inp.temperature,
            inp.scissor,
            //inp.lcycle,
        )
    }

    fn init_with_nac(
        nac:           &Nac,
        efield:        Option<Efield>,
        basis_up:      [usize; 2],
        basis_dn:      [usize; 2],
        dt:            f64,
        namdinit:      usize,
        iniband:       usize,
        inispin:       usize,
        namdtime:      usize,
        nelm:          usize,
        temperature:   f64,
        scissor:       Option<f64>,
        //efield_lcycle: bool,
        ) -> Self {
        
        assert!(namdtime > 1,        "!!!!!\nnamdtime should be greater than 1.\n!!!!!");
        assert!(inispin < nac.nspin, "!!!!!\ninispin should be less or equal than ISPIN of WAVECAR.\n!!!!!");
        assert!(dt >= 0.01,          "!!!!!\ndt should be equal or greater than 0.01fs.\n!!!!!");
        assert!(nelm >= 1,           "!!!!!\nnelm should be equal or greater than 1.\n!!!!!");


        // check basis
        let mut nb = [0usize; 2];   // number of selected bands from each spin
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
        if nb[0] != 0 {
            assert!(
                (basis_up[0] >= nac.brange[0]) && (basis_up[1] <= nac.brange[1]),
                "Selected spin up band range overflows: ({}, {}) not in ({}, {})",
                basis_up[0],   basis_up[1],
                nac.brange[0], nac.brange[1],
                );
        }
        if nb[1] != 0 {
            assert!(
                (basis_dn[0] >= nac.brange[0]) && (basis_dn[1] <= nac.brange[1]),
                "Selected spin down band range overflows: ({}, {}) not in ({}, {})",
                basis_dn[0],   basis_dn[1],
                nac.brange[0], nac.brange[1],
                );
        }
        let nbasis = nb[0] + nb[1];
        assert!(nbasis > 1, "nbasis too small: {}, please use a larger basis set.", nbasis);

        assert!(temperature > 0.001, "Temperature too low: {}, please use higher temperature.", temperature);

        let nsw   = nac.nsw;
        let lreal = nac.lreal;
        let edt   = dt / nelm as f64;
        let nions = nac.proj.shape()[3];
        let nproj = nac.proj.shape()[4];

        let     psi_p = Array1::<c64>::zeros(nbasis);
        let mut psi_c = psi_p.clone();
        let     psi_n = psi_p.clone();
        let mut psi_t = Array2::<c64>::zeros((namdtime, nbasis));
        let     pop_t = Array2::<f64>::zeros((namdtime, nbasis));
        let     psi_h = psi_p.clone();

        let     hamil     = Array2::<c64>::zeros((nbasis, nbasis));
        let mut eig_t     = Array2::<f64>::zeros((nsw-1, nbasis));
        let     prop_eigs = Array1::<f64>::zeros(namdtime);
        let mut nac_t     = Array3::<c64>::zeros((nsw-1, nbasis, nbasis));
        let mut pij_t     = Array4::<c64>::zeros((nsw-1, 3, nbasis, nbasis));
        let mut proj      = Array4::<f64>::zeros((nsw-1, nbasis, nions, nproj));

        let basisini = Self::iniband_index_convert(&basis_up, &basis_dn, inispin, iniband);
        psi_c[basisini] = c64::new(1.0, 0.0);
        psi_t.slice_mut(s![0, ..]).assign(&psi_c);

        let bup = [basis_up[0] - nac.brange[0], basis_up[1] - nac.brange[0]];
        let bdn = [basis_dn[0] - nac.brange[0], basis_dn[1] - nac.brange[0]];

        if nb[1] == 0 {         // spin up only
            eig_t.assign(&nac.eigs .slice(s![.., 0,     bup[0] ..= bup[1]]));
            nac_t.assign(&nac.olaps.slice(s![.., 0,     bup[0] ..= bup[1], bup[0] ..= bup[1]]));
            pij_t.assign(&nac.pij  .slice(s![.., 0, .., bup[0] ..= bup[1], bup[0] ..= bup[1]]));
            proj .assign(&nac.proj .slice(s![.., 0,     bup[0] ..= bup[1], .., ..]));
        } else if nb[0] == 0 {  // spin down only
            eig_t.assign(&nac.eigs .slice(s![.., 1,     bdn[0] ..= bdn[1]]));
            nac_t.assign(&nac.olaps.slice(s![.., 1,     bdn[0] ..= bdn[1], bdn[0] ..= bdn[1]]));
            pij_t.assign(&nac.pij  .slice(s![.., 1, .., bdn[0] ..= bdn[1], bdn[0] ..= bdn[1]]));
            proj .assign(&nac.proj .slice(s![.., 1,     bdn[0] ..= bdn[1], .., ..]));
        } else {                // spin up and spin down
            eig_t.slice_mut(s![.., .. nb[0]]).assign(&nac.eigs.slice(s![.., 0, bup[0] ..= bup[1]]));
            eig_t.slice_mut(s![.., nb[0] ..]).assign(&nac.eigs.slice(s![.., 1, bdn[0] ..= bdn[1]]));

            nac_t.slice_mut(s![.., .. nb[0], .. nb[0]]).assign(&nac.olaps.slice(s![.., 0, bup[0] ..= bup[1], bup[0] ..= bup[1]]));
            nac_t.slice_mut(s![.., nb[0] .., nb[0] ..]).assign(&nac.olaps.slice(s![.., 1, bdn[0] ..= bdn[1], bdn[0] ..= bdn[1]]));

            pij_t.slice_mut(s![.., .., .. nb[0], .. nb[0]]).assign(&nac.pij.slice(s![.., 0, .., bup[0] ..= bup[1], bup[0] ..= bup[1]]));
            pij_t.slice_mut(s![.., .., nb[0] .., nb[0] ..]).assign(&nac.pij.slice(s![.., 1, .., bdn[0] ..= bdn[1], bdn[0] ..= bdn[1]]));

            proj.slice_mut(s![.., .. nb[0], .., ..]).assign(&nac.proj.slice(s![.., 0, bup[0] ..= bup[1], .., ..]));
            proj.slice_mut(s![.., nb[0] .., .., ..]).assign(&nac.proj.slice(s![.., 1, bdn[0] ..= bdn[1], .., ..]));
        }

        // NAC = -ihbar * <i|d/dt|j> / 2, in eV
        nac_t *= c64::new(0.0, -1.0) * HBAR / (2.0 * dt);
        eig_t -= nac.efermi;

        // apply the scissor operator
        if let Some(scissor) = scissor {
            info!("Applying scissor operator of {} eV ...", scissor);
            eig_t.mapv_inplace(|e| if e > 0.0 { e + scissor } else { e });
        }

        // initial auxiliary vars
        let delta_eig    = Array1::<f64>::zeros(nbasis);
        let delta_nac    = Array2::<c64>::zeros((nbasis, nbasis));
        let delta_pij    = Array3::<c64>::zeros((3, nbasis, nbasis));
        let lmi_t        = Array3::<c64>::zeros((namdtime, nbasis, nbasis));
        let ham_t        = Array3::<c64>::zeros((namdtime, nbasis, nbasis));

        Self {
            basis_up,
            basis_dn,
            nbasis,
            dt,
            edt,
            basisini,
            namdinit,
            namdtime,
            nsw,
            nelm,
            lreal,
            temperature,
            //efield_lcycle,

            psi_p,
            psi_c,
            psi_n,
            psi_t,
            pop_t,
            psi_h,

            hamil,
            eig_t,
            prop_eigs,
            nac_t,
            pij_t,
            efield,
            proj,

            delta_eig,
            delta_nac,
            delta_pij,
            lmi_t,
            ham_t,
        }
    }


    #[instrument(skip(self, method))]
    pub fn propagate(&mut self, iion: usize, method: PropagateMethod) {
        let (rtime, _xtime) = self.get_rtime_xtime(iion);

        self.pop_t.slice_mut(s![iion, ..]).assign(&self.psi_c.mapv(|v| v.norm_sqr()));
        self.psi_t.slice_mut(s![iion, ..]).assign(&self.psi_c);
        self.prop_eigs[iion] = (self.pop_t.slice(s![iion, ..]).to_owned() * self.eig_t.slice(s![rtime, ..])).sum();

        match method {
            PropagateMethod::FiniteDifference => self.propagate_finite_difference(iion),
            PropagateMethod::Exact            => self.propagate_exact(iion),
            PropagateMethod::Expm             => self.propagate_expm(iion),
            PropagateMethod::LiouvilleTrotter => self.propagate_liouvilletrotter(iion),
        }

        let norm = self.psi_c.mapv(|x| x.norm_sqr()).sum();
        assert!( (norm - 1.0).abs() < 1E-3 , "Propagation failed, norm not conserved: norm = {}", norm);
    }


    fn propagate_finite_difference(&mut self, iion: usize) {
        for iele in 0 .. self.nelm {
            self.make_hamil(iion, iele);
            self.psi_h = self.hamil.dot(&self.psi_c);

            if 0 == iion && 0 == iele {
                self.psi_n = self.psi_c.clone() -       IMGUNIT * self.edt / HBAR * self.psi_h.clone();
            } else {
                self.psi_n = self.psi_p.clone() - 2.0 * IMGUNIT * self.edt / HBAR * self.psi_h.clone();
            }

            self.psi_p.assign(&self.psi_c);
            self.psi_c.assign(&self.psi_n);
        }
    }


    // Perform |psi'> = exp(-iHt/hbar) |psi>
    fn propagate_exact(&mut self, iion: usize) {
        for iele in 0 .. self.nelm {
            self.make_hamil(iion, iele);
            self.hamil.mapv_inplace(|v| v * (-self.edt / HBAR)); // -edt*H/hbar is still hermitian

            // P, Lambda = eigh(-edt*H/hbar)
            let (eigvals, eigvecs) = self.hamil.eigh_inplace(UPLO::Upper).unwrap();
            let expie = eigvals.mapv(|v| c64::new(v.cos(), v.sin()));
            self.psi_c.assign(&eigvecs.dot(&Array2::from_diag(&expie))
                                      .dot(&eigvecs.t().mapv(|v| v.conj()))
                                      .dot(&self.psi_c));
        }
    }


    fn propagate_expm(&mut self, iion: usize) {
        for iele in 0 .. self.nelm {
            self.make_hamil(iion, iele);
            self.hamil.mapv_inplace(|x| x * (-self.edt / HBAR) * IMGUNIT);
            self.psi_c.assign(&expm(&self.hamil).expect("Matrix exponentiation failed during propagate_expm.")
                                                .dot(&self.psi_c));
        }
    }


    // WARN: this method is not tested yet
    fn propagate_liouvilletrotter(&mut self, iion: usize) {
        assert!(self.lreal, "LiouvilleTrotter method can be used for REAL NAC only.");
        assert!(self.efield.is_none(), "LiouvilleTrotter method cannot be used with EFIELD present.");

        // Don't initialize them, in order to make rustc happy
        let mut cjj:     c64;
        let mut ckk:     c64;
        let mut phi:     f64;
        let mut cos_phi: f64;
        let mut sin_phi: f64;

        for iele in 0 .. self.nelm {
            self.make_hamil(iion, iele);
            for jj in 0 .. self.nbasis {        // the traversal order is not tested
                for kk in jj+1 .. self.nbasis {
                    phi     = -(IMGUNIT * self.hamil[(kk, jj)]).re * 0.5 * self.edt / HBAR;
                    cos_phi = phi.cos();
                    sin_phi = phi.sin();
                    cjj     = self.psi_c[jj];
                    ckk     = self.psi_c[kk];
                    self.psi_c[jj] =  cos_phi * cjj + sin_phi * ckk;
                    self.psi_c[kk] = -sin_phi * cjj + cos_phi * ckk;
                }
            }

            for jj in 0 .. self.nbasis {
                phi = -(IMGUNIT * self.hamil[(jj, jj)]).re * 0.5 * self.edt / HBAR;
                self.psi_c[jj] = self.psi_c[jj] * phi.exp();
            }


            for jj in (0 .. self.nbasis).rev() {
                for kk in (jj+1 .. self.nbasis).rev() {
                    phi     = -(IMGUNIT * self.hamil[(kk, jj)]).re * 0.5 * self.edt / HBAR;
                    cos_phi = phi.cos();
                    sin_phi = phi.sin();
                    cjj     = self.psi_c[jj];
                    ckk     = self.psi_c[kk];
                    self.psi_c[jj] =  cos_phi * cjj + sin_phi * ckk;
                    self.psi_c[kk] = -sin_phi * cjj + cos_phi * ckk;
                }
            }
        } // iele
    }


    fn _propagate_runge_kutta_4(&mut self, _iion: usize) {
        todo!("Runge-Kutta method not implemented");
    }


    fn _propagate_crank_nicolson(&mut self, _iion: usize) {
        todo!("Crank-Nicolson method not implemented");
    }


    #[instrument(skip(self), level="info")]
    fn make_hamil(&mut self, iion: usize, iele: usize) {
        let (rtime, xtime) = self.get_rtime_xtime(iion);

        // first electronic step inside ionic step
        if 0 == iele {
            self.delta_eig = (self.eig_t.slice(s![xtime, ..]).to_owned() -
                              self.eig_t.slice(s![rtime, ..])) / self.nelm as f64 ;
            self.delta_nac = (self.nac_t.slice(s![xtime, .., ..]).to_owned() -
                              self.nac_t.slice(s![rtime, .., ..])) / self.nelm as f64;
            self.delta_pij = (self.pij_t.slice(s![xtime, .., .., ..]).to_owned() -
                              self.pij_t.slice(s![rtime, .., .., ..]) ) / self.nelm as f64;
        }

        // non-diagonal part: NAC
        // nac_t is in eV alrady
        self.hamil = self.nac_t.slice(s![rtime, .., ..]).to_owned() +
                     self.delta_nac.to_owned() * (iele as f64);

        //info!("LINE = {}, RTIME = {}, XTIME = {}", line!(), rtime, xtime);

        // if electric field exists
        if self.efield.is_some() {
            // light matter interaction
            // afield .dot. <i|p|j> / m_e, in eV
            let afield  = self.efield.as_ref().map(|e| e.get_afield(iion, iele)).unwrap_or([0.0; 3]);
            let mut lmi = Array2::<c64>::zeros((self.nbasis, self.nbasis));
            for ii in 0 .. 3 {
                lmi += &((self.pij_t.slice(s![rtime, ii, .., ..]).to_owned() +
                          self.delta_pij.slice(s![ii, .., ..]).to_owned() * (iele as f64)) * afield[ii]);
            }

            lmi.map_inplace(|x| { *x *= FACT_PA; });

            self.hamil += &lmi;

            if 0 == iele {
                self.lmi_t.slice_mut(s![iion, .., ..]).assign(&lmi);
            }
        }

        if 0 == iele {
            self.ham_t.slice_mut(s![iion, .., ..]).assign(&self.hamil);
        }

        // dagonal part: eigenvalue of ks orbits, in eV
        // rustc refuses to compile `struct.method() = somethingelse;`
        self.hamil.diag_mut().assign(&(
            self.eig_t.slice(s![rtime, ..]).mapv(|v| c64::new(v, 0.0)) +
            self.delta_eig.to_owned() * (iele as f64)
            ));
    }


    pub fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let f = H5File::create(fname)?;

        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.re)).create("nac_t_r")?;
        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.im)).create("nac_t_i")?;

        f.new_dataset_builder().with_data(&self.eig_t).create("eig_t")?;

        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_r")?;
        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_i")?;

        f.new_dataset_builder().with_data(&self.proj).create("proj")?;

        Ok(())
    }


    #[instrument(ret, level="debug")]
    fn iniband_index_convert(
        basis_up: &[usize; 2], basis_dn: &[usize; 2], inispin: usize, iniband: usize
        ) -> usize {

        if 0 == inispin {
            iniband - basis_up[0]
        } else {
            let nbup = basis_up[1] - basis_up[0] + 1;
            iniband - basis_dn[0] + nbup
        }
    }


    pub fn get_rtime_xtime(&self, iion: usize) -> (usize, usize) {
        let rtime = (iion + self.namdinit) % (self.nsw - 2);
        let xtime = rtime + 1;
        (rtime, xtime)
    }
}
