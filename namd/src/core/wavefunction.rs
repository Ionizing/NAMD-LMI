use shared::ndarray as nd;
use crate::core::Hamiltonian;

pub trait Wavefunction {
    type DataType;
    type ArrayType<'a>: nd::AsArray<'a, Self::DataType, nd::Ix1> where Self: 'a;     // [nbasis]
    type TdArrayType<'a>: nd::AsArray<'a, Self::DataType, nd::Ix2> where Self: 'a;     // [time, nbasis]
    type PopArraryType<'a>: nd::AsArray<'a, f64, nd::Ix1> where Self: 'a;
    type TdPopArrayType<'a>: nd::AsArray<'a, f64, nd::Ix2> where Self: 'a;
    type TdEigArrayType<'a>: nd::AsArray<'a, f64, nd::Ix1> where Self: 'a;
    type BasisIni<'a>: nd::AsArray<'a, usize, nd::Ix1> where Self: 'a;
    type HamiltonianType: Hamiltonian;

    fn get_nbasis(&self) -> usize;
    fn get_basisini(&self) -> Self::BasisIni<'_>;
    fn get_namdinit(&self) -> usize;
    fn get_namdtime(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_nelm(&self) -> usize;
    fn get_nspin(&self) -> usize;
    fn get_lncl(&self) -> bool;

    fn propagate_full(&mut self, hamil: &Self::HamiltonianType);
    fn get_psi(&self, iion: usize) -> Self::ArrayType<'_>;       // [time, nbasis]
    fn get_psi_t(&self) -> Self::TdArrayType<'_>;                // [time, nbasis]

    fn get_pop(&self, iion: usize) -> Self::PopArraryType<'_>;   // [nbasis]
    fn get_pop_t(&self) -> Self::TdPopArrayType<'_>;             // [time, nbasis]

    fn get_prop_eigs(&self) -> Self::TdEigArrayType<'_>;         // [time, nbasis]
}
