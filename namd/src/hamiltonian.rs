use shared::{
    ndarray::s,
    c64,
    Array1,
    Array2,
    Array3,
    ndarray::Array4,
};

use crate::{
    constants::*,
    efield::Efield,
    nac::Nac,
    input::Input,
};


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

    pub psi_p:         Array1<c64>,
    pub psi_c:         Array1<c64>,
    pub psi_n:         Array1<c64>,
    pub psi_t:         Array2<c64>,
    pub pop_t:         Array2<f64>,
    pub psi_h:         Array1<c64>,
    
    pub hamil:         Array2<c64>,
    pub eig_t:         Array2<f64>,
    pub prop_eigs:     Array1<f64>,
    pub nac_t:         Array3<c64>,
    pub pij_t:         Array4<c64>,
    pub efield:        Option<Efield>,
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

        let basisini = Self::_iniband_index_convert(&basis_up, &basis_dn, inispin, iniband);
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
        }
    }


    fn _iniband_index_convert(
        basis_up: &[usize; 2], basis_dn: &[usize; 2], inispin: usize, iniband: usize
        ) -> usize {

        if 0 == inispin {
            iniband - basis_up[0]
        } else {
            let nbup = basis_up[1] - basis_up[0];
            iniband - basis_dn[0] + nbup
        }
    }

}
