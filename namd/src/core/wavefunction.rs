use shared::ndarray as nd;
use crate::core::Hamiltonian;

pub trait Wavefunction {
    type DataType;
    type ArrayType<'a>: nd::AsArray<'a, Self::DataType, nd::Ix1> where Self: 'a;     // [nbasis]
    type TdArrayType<'a>: nd::AsArray<'a, Self::DataType, nd::Ix2> where Self: 'a;     // [time, nbasis]
    type PopArraryType<'a>: nd::AsArray<'a, f64, nd::Ix1> where Self: 'a;
    type TdPopArrayType<'a>: nd::AsArray<'a, f64, nd::Ix2> where Self: 'a;
    type TdEigArrayType<'a>: nd::AsArray<'a, f64, nd::Ix1> where Self: 'a;

    fn get_nbasis(&self) -> usize;
    fn get_basisini(&self) -> usize;
    fn get_namdinit(&self) -> usize;
    fn get_namdtime(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_nelm(&self) -> usize;

    fn propagate_one_step(&mut self, hamil: &impl Hamiltonian, iion: usize);
    fn propagate_full(&mut self, hamil: &impl Hamiltonian);
    fn get_psi<'a>(&'a self, iion: usize) -> Self::ArrayType<'a>;       // [time, nbasis]
    fn get_psi_t<'a>(&'a self) -> Self::TdArrayType<'a>;                // [time, nbasis]

    fn get_pop<'a>(&'a self, iion: usize) -> Self::PopArraryType<'a>;   // [nbasis]
    fn get_pop_t<'a>(&'a self) -> Self::TdPopArrayType<'a>;             // [time, nbasis]

    fn get_prop_eigs<'a>(&'a self) -> Self::TdEigArrayType<'a>;         // [time, nbasis]
}
