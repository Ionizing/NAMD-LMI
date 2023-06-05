use std::path::Path;

use shared::{
    Result,
    ndarray::s,
    c64,
    Array1,
    Array2,
    Array3,
    ndarray::Array4,
};
use hdf5::{
    File as H5File,
    H5Type,
    Selection,
};

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


pub struct Hamiltonian {
    pub basis_up:      [usize; 2],
    pub basis_dn:      [usize; 2],
    pub nbasis:        usize,
    pub dt:            f64,
    pub basisini:      usize,
    pub namdinit:      usize,
    pub namdtime:      usize,
    pub nsw:           usize,
    pub nelm:          usize,
    pub lreal:         bool,
    pub temperature:   f64,
    pub efield_lcycle: bool,

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

    // Auxiliary variables used by `make_hamil method`
    delta_eig:         Array1<f64>,     // [nbasis]
    delta_nac:         Array2<c64>,     // [nbasis, nbasis]
    delta_pij:         Array3<c64>,     // [3, nbasis, nbasis]
}


impl Hamiltonian {
    pub fn from_input(nac: &Nac, inp: &Input, inicon_idx: usize) -> Self {
        Self::init_with_nac(
            nac,
            inp.efield.clone().map(|x| x.1),
            inp.basis_up.clone(),
            inp.basis_dn.clone(),
            inp.dt,
            inp.inisteps[inicon_idx],
            inp.inibands[inicon_idx],
            inp.inispins[inicon_idx],
            inp.namdtime,
            inp.nelm,
            inp.temperature,
            inp.scissor,
            inp.lcycle,
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
        scissor:       f64,
        efield_lcycle: bool,
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
            basis_up[1] - basis_up[0]
        };
        nb[1] = if basis_dn.contains(&0) {
            0
        } else {
            basis_dn[1] - basis_dn[0]
        };
        if nb[0] != 0 {
            assert!(
                (basis_up[0] >= nac.brange[0]) && (basis_up[1] <= nac.brange[1]),
                "Selected spin up band range overflows: ({}, {}) not in ({}, {})",
                basis_up[0]+1, basis_up[1],
                nac.brange[0]+1, nac.brange[1],
                );
        }
        if nb[1] != 0 {
            assert!(
                (basis_dn[0] >= nac.brange[0]) && (basis_dn[1] <= nac.brange[1]),
                "Selected spin down band range overflows: ({}, {}) not in ({}, {})",
                basis_dn[0]+1, basis_dn[1],
                nac.brange[0]+1, nac.brange[1],
                );
        }
        let nbasis = nb[0] + nb[1];
        assert!(nbasis <= 1, "nbasis too small: {}, please use a larger basis set.", nbasis);

        assert!(temperature > 0.001, "Temperature too low: {}, please use higher temperature.", temperature);

        let nsw   = nac.nsw;
        let lreal = nac.lreal;

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

        let basisini = Self::iniband_index_convert(&basis_up, &basis_dn, inispin, iniband);
        psi_c[basisini] = c64::new(1.0, 0.0);
        psi_t.slice_mut(s![0, ..]).assign(&psi_c);

        let bup = [basis_up[0] - nac.brange[0], basis_up[1] - nac.brange[1]];
        let bdn = [basis_dn[0] - nac.brange[0], basis_dn[1] - nac.brange[1]];

        if nb[1] == 0 {         // spin up only
            eig_t.assign(&nac.eigs .slice(s![.., 0,     bup[0] .. bup[1]]));
            nac_t.assign(&nac.olaps.slice(s![.., 0,     bup[0] .. bup[1], bup[0] .. bup[1]]));
            pij_t.assign(&nac.pij  .slice(s![.., 0, .., bup[0] .. bup[1], bup[0] .. bup[1]]));
        } else if nb[0] == 0 {  // spin down only
            eig_t.assign(&nac.eigs .slice(s![.., 1,     bdn[0] .. bdn[1]]));
            nac_t.assign(&nac.olaps.slice(s![.., 1,     bdn[0] .. bdn[1], bdn[0] .. bdn[1]]));
            pij_t.assign(&nac.pij  .slice(s![.., 1, .., bdn[0] .. bdn[1], bdn[0] .. bdn[1]]));
        } else {                // spin up and spin down
            eig_t.slice_mut(s![.., .. nb[0]]).assign(&nac.eigs.slice(s![.., 0, bup[0] .. bup[1]]));
            eig_t.slice_mut(s![.., nb[0] ..]).assign(&nac.eigs.slice(s![.., 1, bdn[0] .. bdn[1]]));

            nac_t.slice_mut(s![.., .. nb[0], .. nb[0]]).assign(&nac.olaps.slice(s![.., 0, bup[0] .. bup[1], bup[0] .. bup[1]]));
            nac_t.slice_mut(s![.., nb[0] .., nb[0] ..]).assign(&nac.olaps.slice(s![.., 1, bdn[0] .. bdn[1], bdn[0] .. bdn[1]]));

            pij_t.slice_mut(s![.., .., .. nb[0], .. nb[0]]).assign(&nac.pij.slice(s![.., 0, .., bup[0] .. bup[1], bup[0] .. bup[1]]));
            pij_t.slice_mut(s![.., .., nb[0] .., nb[0] ..]).assign(&nac.pij.slice(s![.., 1, .., bdn[0] .. bdn[1], bdn[0] .. bdn[1]]));
        }

        nac_t *= c64::new(0.0, -1.0).scale(HBAR / (2.0 * dt));
        eig_t -= nac.efermi;

        // apply the scissor operator
        eig_t.mapv_inplace(|e| if e > 0.0 { e + scissor } else { e });

        // initial auxiliary vars
        let delta_eig    = Array1::<f64>::zeros(nbasis);
        let delta_nac    = Array2::<c64>::zeros((nbasis, nbasis));
        let delta_pij    = Array3::<c64>::zeros((3, nbasis, nbasis));

        Self {
            basis_up,
            basis_dn,
            nbasis,
            dt,
            basisini,
            namdinit,
            namdtime,
            nsw,
            nelm,
            lreal,
            temperature,
            efield_lcycle,

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

            delta_eig,
            delta_nac,
            delta_pij,
        }
    }


    pub fn propagate(&mut self, iion: usize, method: PropagateMethod) {
        todo!()
    }


    fn propagate_finite_difference(&mut self, iion: usize) {
        todo!()
    }


    fn propagate_exact(&mut self, iion: usize) {
        todo!()
    }


    fn propagate_expm(&mut self, iion: usize) {
        todo!()
    }


    fn propagate_liouvilletrotter(&mut self, iion: usize) {
        todo!()
    }


    fn make_hamil(&mut self, iion: usize, iele: usize) {
        let rtime: usize = (iion + self.namdinit) % (self.nsw - 1);
        let xtime: usize = rtime + 1;

        // first electronic step inside ionic step
        if 0 == iele {
            self.delta_eig = (self.eig_t.slice(s![xtime, ..]).to_owned() -
                              self.eig_t.slice(s![rtime,   ..]) ) / self.nelm as f64 ;
            self.delta_nac = (self.nac_t.slice(s![xtime, .., ..]).to_owned() -
                              self.nac_t.slice(s![rtime,   .., ..]) ) / self.nelm as f64;
            self.delta_pij = (self.pij_t.slice(s![xtime, .., .., ..]).to_owned() -
                              self.pij_t.slice(s![rtime,   .., .., ..]) ) / self.nelm as f64;
        }

        // non-diagonal part: NAC
        self.hamil = self.nac_t.slice(s![rtime, .., ..]).to_owned() +
                     self.delta_nac.mapv(|v| v.scale(iele as f64));

        // light matter interaction
        // vecpot dot <i|p|j>
        let vecpot = self.get_vecpot(iion, iele);
        for i in 0 .. 3 {
            self.hamil += &(
                self.pij_t.slice(s![rtime, i, .., ..]).to_owned() +
                self.delta_pij.mapv(|v| v.scale(iele as f64 * vecpot[i]))
                );
        }

        // dagonal part: eigenvalue of ks orbits
        // rustc refuses to compile `struct.method() = somethingelse;`
        self.hamil.diag_mut().assign(&(
            self.eig_t.slice(s![rtime, ..]).mapv(|v| c64::new(v, 0.0)) +
            self.delta_eig.mapv(|v| c64::new(v * iele as f64, 0.0))
            ));
    }


    pub fn save_to_h5<P>(&self, fname: &P) -> Result<()>
    where
        P: AsRef<Path> + ?Sized,
    {
        let f = H5File::create(fname)?;

        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.re)).create("nac_t_r")?;
        f.new_dataset_builder().with_data(&self.nac_t.mapv(|v| v.im)).create("nac_t_i")?;

        f.new_dataset_builder().with_data(&self.eig_t).create("eig_t")?;

        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_r")?;
        f.new_dataset_builder().with_data(&self.pij_t.mapv(|v| v.im)).create("pij_t_i")?;

        Ok(())
    }


    fn iniband_index_convert(
        basis_up: &[usize; 2], basis_dn: &[usize; 2], inispin: usize, iniband: usize
        ) -> usize {

        if 0 == inispin {
            iniband - basis_up[0]
        } else {
            let nbup = basis_up[1] - basis_up[0];
            iniband - basis_dn[0] + nbup
        }
    }


    fn get_vecpot(&self, iion: usize, iele: usize) -> [f64; 3] {
        match self.efield.as_ref().unwrap() {
            Efield::Vector3(v) => {
                let mut ret = [0.0f64; 3];
                for i in 0 .. 3 {
                    ret[i] = (v[iion+1][i] - v[iion][i]) * (iele as f64) / (self.nelm as f64);
                }
                ret
            },
            Efield::Function3(f) => {
                // convert indices to real time value
                let t = iion as f64 * self.dt + (iele as f64 / self.nelm as f64) * self.dt;
                [
                    f[0].eval(t),
                    f[1].eval(t),
                    f[2].eval(t),
                ]
            }
        }
    }

}
