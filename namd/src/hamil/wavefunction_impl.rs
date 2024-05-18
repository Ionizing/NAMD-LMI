use shared::c64;
use shared::ndarray as nd;

use crate::core::Hamiltonian;
use crate::core::Wavefunction;

pub struct SPWavefunction {
    nbasis: usize,
    basisini: usize,
    namdinit: usize,
    namdtime: usize,
    potim: f64,
    nelm: usize,

    psi: nd::Array1<c64>,   // [nbasis]
    psi_t: nd::Array2<c64>, // [namdtime, nbasis]
    pop_t: nd::Array2<f64>, // [namdtime, nbasis]
    eig_t: nd::Array1<f64>, // [namdtime]
}


impl Wavefunction for SPWavefunction {
    type DataType = c64;
    type ArrayType<'a> = nd::ArrayView1<'a, c64>;
    type TdArrayType<'a> = nd::ArrayView2<'a, c64>;
    type PopArraryType<'a> = nd::ArrayView1<'a, f64>;
    type TdPopArrayType<'a> = nd::ArrayView2<'a, f64>;
    type TdEigArrayType<'a> = nd::ArrayView1<'a, f64>;

    fn get_nbasis(&self) -> usize { self.nbasis }
    fn get_basisini(&self) -> usize { self.basisini }
    fn get_namdinit(&self) -> usize { self.namdinit }
    fn get_namdtime(&self) -> usize { self.namdtime }
    fn get_potim(&self) -> f64 { self.potim }
    fn get_nelm(&self) -> usize { self.nelm }

    fn propagate_one_step(&mut self, hamil: &impl Hamiltonian, iion: usize) {
        todo!()
    }
    fn propagate_full(&mut self, hamil: &impl Hamiltonian){
        todo!()
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
