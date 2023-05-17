use shared::{
    c64,
    Array1,
    Array2,
    Array3,
    ndarray::Array4,
};


use crate::nac::Nac;


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
    pub prop_eigs:     Array2<f64>,
    pub nac_t:         Array3<c64>,
    pub pij_t:         Array4<c64>,
    pub efield:        Array2<f64>,
}


impl Hamiltonian {
    fn from_params(
        nac:           &Nac,
        efield:        &Array2<f64>,
        basis_up:      &[usize; 2],
        basis_dn:      &[usize; 2],
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
        todo!()
    }

}
