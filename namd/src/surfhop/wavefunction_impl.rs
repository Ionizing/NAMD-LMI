use shared::Result;
use shared::c64;
use shared::ndarray as nd;
use shared::MatX3;
use shared::ndarray_linalg::{
    EighInplace,
    UPLO,
    expm,
};

use crate::core::constants::*;
use crate::core::Wavefunction;
use crate::core::Hamiltonian;
use crate::hamil::{
    PropagateMethod,
    Efield,
    SPHamiltonian,
};

pub struct SPWavefunction {
    nbasis: usize,
    basisini: usize,
    namdinit: usize,
    namdtime: usize,
    potim: f64,
    nelm: usize,

    psi0: nd::Array1<c64>,  // [nbasis]
    psi_t: nd::Array2<c64>, // [namdtime, nbasis]
    pop_t: nd::Array2<f64>, // [namdtime, nbasis]
    eig_t: nd::Array1<f64>, // [namdtime]

    efield_array: Option<[MatX3<f64>;2]>,
}


impl Wavefunction for SPWavefunction {
    type DataType = c64;
    type ArrayType<'a> = nd::ArrayView1<'a, c64>;
    type TdArrayType<'a> = nd::ArrayView2<'a, c64>;
    type PopArraryType<'a> = nd::ArrayView1<'a, f64>;
    type TdPopArrayType<'a> = nd::ArrayView2<'a, f64>;
    type TdEigArrayType<'a> = nd::ArrayView1<'a, f64>;
    type HamiltonianType = SPHamiltonian;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_basisini(&self) -> usize { self.basisini }
    fn get_namdinit(&self) -> usize { self.namdinit }
    fn get_namdtime(&self) -> usize { self.namdtime }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nelm(&self) -> usize { self.nelm }

    fn propagate_full(&mut self, hamil: &SPHamiltonian) {
        let edt = self.potim / self.nelm as f64;
        self.psi_t.slice_mut(nd::s![0, ..]).assign(&self.psi0);
        self.pop_t.slice_mut(nd::s![0, ..]).assign(&self.psi0.mapv(|v| v.norm_sqr()));
        self.eig_t[0] = {
            let hamil_diag = hamil.get_hamil0_rtime(0, self.namdinit).diag().mapv(|v| v.re);
            (hamil_diag * self.pop_t.slice(nd::s![0, ..])).sum()
        };


        // initial psi
        let mut psi = self.psi0.clone();

        // construct hamiltonian with electron-photon interaction
        for iion in 0 .. self.namdtime - 1 {
            let hamil_i = hamil.get_hamil0_rtime(iion, self.namdinit).to_owned();
            let hamil_j = hamil.get_hamil0_rtime(iion + 1, self.namdinit).to_owned();
            let delta_hamil = (hamil_j - hamil_i.view()) / (self.nelm as f64);

            for ielm in 0 .. self.nelm {
                let lmi   = self.get_lmi(hamil, iion, ielm);

                // This matrix is re-constructed in every loop, thus propagate_s can take its
                // ownersip.
                let hamil_t = hamil_i.to_owned() + delta_hamil.to_owned() * ielm as f64 + lmi;
                Self::propagate_dispatch(hamil_t, &mut psi, edt, hamil.get_propmethod());
            }

            self.psi_t.slice_mut(nd::s![iion+1, ..]).assign(&psi);
            let pop_c = psi.mapv(|v| v.norm_sqr());
            self.pop_t.slice_mut(nd::s![iion+1, ..]).assign(&pop_c);
        }
    }

    fn get_psi(&self, iion: usize) -> Self::ArrayType<'_> {
        self.psi_t.slice(nd::s![iion, ..])
    }
    fn get_psi_t(&self) -> Self::TdArrayType<'_> {
        self.psi_t.view()
    }

    fn get_pop(&self, iion: usize) -> Self::PopArraryType<'_> {
        self.pop_t.slice(nd::s![iion, ..])
    }
    fn get_pop_t(&self) -> Self::TdPopArrayType<'_> {
        self.pop_t.view()
    }

    fn get_prop_eigs(&self) -> Self::TdEigArrayType<'_> {
        self.eig_t.view()
    }
}


impl SPWavefunction {
    pub fn calculate_efield_array(&mut self, efield: &mut Efield) {
        if self.efield_array.is_some() {
            return;
        }
        let (_tt, efield_array) = efield.get_eafield_array(self.namdtime, self.potim, self.nelm);
        self.efield_array = Some(efield_array);
    }


    pub fn get_efield_array(&self) -> Option<&[MatX3<f64>; 2]> {
        self.efield_array.as_ref()
    }


    fn propagate_dispatch(hamil: nd::Array2<c64>, psi: &mut nd::Array1<c64>, edt: f64, method: PropagateMethod) {
        use PropagateMethod::*;

        match method {
            Expm => Self::propagate_expm(hamil, psi, edt),
            Exact => Self::propagate_exact(hamil, psi, edt),
            FiniteDifference => Self::propagate_fd(hamil, psi, edt),
            LiouvilleTrotter => Self::propagate_lt(hamil, psi, edt),
        }
    }


    fn propagate_expm(mut hamil: nd::Array2<c64>, psi: &mut nd::Array1<c64>, edt: f64) {
        hamil *= edt / HBAR * IMGUNIT;
        psi.assign(&expm(&hamil)
            .expect("Matrix exponentiation failed during propagate_expm.")
            .dot(psi));
    }


    fn propagate_fd(mut hamil: nd::Array2<c64>, psi: &mut nd::Array1<c64>, edt: f64) {
        todo!("To be implemented.")
    }


    fn propagate_exact(mut hamil: nd::Array2<c64>, psi: &mut nd::Array1<c64>, edt: f64) {
        hamil *= edt / HBAR * IMGUNIT;
        let (eigvals, eigvecs) = hamil.eigh_inplace(UPLO::Upper).unwrap();
        let expie = eigvals.mapv(|v| c64::new(v.cos(), v.sin()));
        psi.assign(&eigvecs.dot(&nd::Array2::from_diag(&expie))
                           .dot(&eigvecs.t().mapv(|v| v.conj()))
                           .dot(psi));
    }


    fn propagate_lt(mut hamil: nd::Array2<c64>, psi: &mut nd::Array1<c64>, edt: f64) {
        todo!("To be implemented.")
    }


    // get ep*A/m where p is <i|p|j>, A is vector potential; e and m are absorbed in FACT_PA
    pub fn get_lmi(&self, hamil: &SPHamiltonian, iion: usize, ielm: usize) -> nd::Array2<c64> {
        let itime = iion * self.nelm + ielm;

        let pij_i  = hamil.get_pij_rtime(iion, self.namdinit).to_owned();
        let pij_j  = hamil.get_pij_rtime(iion, self.namdinit).to_owned();

        let pij = pij_i.clone() + (pij_j - pij_i) / (self.nelm as f64) * (ielm  as f64);
        let afield = self.get_efield_array()
            .map(|[_e, a]| a[itime])
            .unwrap_or([0.0; 3]);

        let mut lmi = nd::Array2::<c64>::zeros((self.nbasis, self.nbasis));
        for i in 0 .. self.nbasis {
            lmi[(i, i)] = c64::new(0.0, 0.0);
            for j in i+1 .. self.nbasis {
                lmi[(i, j)] = pij[(0, i, j)] * afield[0]
                            + pij[(1, i, j)] * afield[1]
                            + pij[(2, i, j)] * afield[2];
                lmi[(j, i)] = lmi[(i, j)].conj();
            }
        }

        lmi * FACT_PA       // unit: eV
    }
}


impl SPWavefunction {
    pub fn from_hamil_and_params(
        hamil: &SPHamiltonian,
        iniband: usize, inispin: usize,
        namdtime: usize, nelm: usize, namdinit: usize
    ) -> Result<Self> {
        let nbasis = hamil.get_nbasis();
        let basisini = hamil.get_converted_index(iniband, inispin)?;
        // namdinit
        // namdtime
        let potim = hamil.get_potim();
        // nelm

        let mut psi0 = nd::Array1::<c64>::zeros(nbasis);
        psi0[basisini] = c64::new(1.0, 0.0);

        let psi_t = nd::Array2::<c64>::zeros((namdtime, nbasis));
        let pop_t = nd::Array2::<f64>::zeros((namdtime, nbasis));
        let eig_t = nd::Array1::<f64>::zeros(namdtime);

        let efield_array = if let Some(efield) = hamil.get_efield() {
            Some(efield.lock().unwrap().get_eafield_array(namdtime, hamil.get_potim(), nelm).1)
        } else {
            None
        };

        Ok(Self {
            nbasis,
            basisini,
            namdinit,
            namdtime,
            potim,
            nelm,

            psi0,
            psi_t,
            pop_t,
            eig_t,
            efield_array,
        })
    }
}
