use std::fmt::Write;

use shared::Result;
use shared::anyhow;
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
    SPHamiltonian,
};
use itertools::Itertools;

pub struct SPWavefunction {
    nbasis: usize,
    basisini: nd::Array1<usize>,
    nbasisini: usize,
    namdinit: usize,
    namdtime: usize,
    potim: f64,
    nelm: usize,
    nspin: usize,
    lncl: bool,
    maxocc: usize,

    psi0: nd::Array1<c64>,  // [nbasis]
    psi_t: nd::Array2<c64>, // [namdtime, nbasis]
    pop_t: nd::Array2<f64>, // [namdtime, nbasis]
    eig_t: nd::Array1<f64>, // [namdtime]

    eafield_array: Option<[MatX3<f64>;2]>,
}


impl Wavefunction for SPWavefunction {
    type DataType = c64;
    type ArrayType<'a> = nd::ArrayView1<'a, c64>;
    type TdArrayType<'a> = nd::ArrayView2<'a, c64>;
    type PopArraryType<'a> = nd::ArrayView1<'a, f64>;
    type TdPopArrayType<'a> = nd::ArrayView2<'a, f64>;
    type TdEigArrayType<'a> = nd::ArrayView1<'a, f64>;
    type BasisIni<'a> = nd::ArrayView1<'a, usize>;
    type HamiltonianType = SPHamiltonian;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_basisini(&self) -> Self::BasisIni<'_> { self.basisini.view() }
    fn get_namdinit(&self) -> usize { self.namdinit }
    fn get_namdtime(&self) -> usize { self.namdtime }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nelm(&self) -> usize { self.nelm }
    fn get_nspin(&self) -> usize { self.nspin }
    fn get_lncl(&self) -> bool { self.lncl }

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
            let eig = (hamil_i.diag().mapv(|v| v.re) * pop_c).sum();
            self.eig_t[iion+1] = eig;
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
    pub fn get_eafield_array(&self) -> Option<&[MatX3<f64>; 2]> {
        self.eafield_array.as_ref()
    }


    pub fn get_nbasisini(&self) -> usize {
        self.basisini.len()
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
        let afield = self.get_eafield_array()
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
        iniband: &[i32],
        namdtime: usize, nelm: usize, namdinit: usize,
        eafield_array: Option<[MatX3<f64>; 2]>,
    ) -> Result<Self> {
        let nspin = hamil.get_nspin();
        let lncl  = hamil.get_lncl();
        let max_occupation = if lncl || nspin == 2 { 1usize } else { 2 };

        // check occupation of initial state
        let mut err = format!("Invalid inibands, some bands have more than {} electrons: ", max_occupation);
        let mut iniband_dedup = iniband.to_owned();
        iniband_dedup.sort();
        let mut fail = false;
        for cnt in iniband_dedup.into_iter().dedup_with_count() {
            if cnt.0 > max_occupation {
                fail = true;
                write!(&mut err, " occ({})={}", cnt.1, cnt.0)?;
            }
        }
        if fail {
            write!(&mut err, ", please check.")?;
            anyhow::bail!(err);
        }
        

        let nbasis = hamil.get_nbasis();
        let basisini = iniband.iter()
            .map(|&ib| -> Result<usize> {
                hamil.get_converted_index(ib)
            })
            .collect::<Result<nd::Array1<usize>>>()?;

        let nbasisini = basisini.len();
        anyhow::ensure!(nbasisini > 0);

        let potim = hamil.get_potim();

        let mut psi0 = nd::Array1::<c64>::zeros(nbasis);
        for &ib in basisini.iter() {
            psi0[ib] += c64::new(1.0, 0.0);
        }
        psi0.mapv(|x| x / (nbasisini as f64).sqrt());       // normalize

        let psi_t = nd::Array2::<c64>::zeros((namdtime, nbasis));
        let pop_t = nd::Array2::<f64>::zeros((namdtime, nbasis));
        let eig_t = nd::Array1::<f64>::zeros(namdtime);

        Ok(Self {
            nbasis,
            basisini,
            nbasisini,
            namdinit,
            namdtime,
            potim,
            nelm,
            nspin,
            lncl,
            maxocc: max_occupation,

            psi0,
            psi_t,
            pop_t,
            eig_t,
            eafield_array,
        })
    }

    pub fn get_max_occupation(&self) -> usize {
        self.maxocc
    }
}
