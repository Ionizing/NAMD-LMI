use shared::ndarray as nd;
use crate::core::Hamiltonian;

pub trait Wavefunction {
    type DataType;
    type ArrayType: for<'a>  nd::AsArray<'a, Self::DataType, nd::Ix1>;     // [nbasis]
    type TdArrayType: for<'a> nd::AsArray<'a, Self::DataType, nd::Ix2>;     // [time, nbasis]

    fn get_nbasis(&self) -> usize;
    fn get_basisini(&self) -> usize;
    fn get_namdinit(&self) -> usize;

    fn propagate(&mut self, hamil: &impl Hamiltonian);
    fn get_psi(&self, iion: usize) -> &Self::ArrayType;
    fn get_psi_t(&self) -> &Self::TdArrayType;

    fn get_pop<A>(&self, iion: usize) -> &A
        where A: for<'a> nd::AsArray<'a, f64, nd::Ix1>;     // [nbasis]
    fn get_pop_t<A>(&self, iion: usize) -> &A
        where A: for<'a> nd::AsArray<'a, f64, nd::Ix2>;     // [time, nbasis]

    fn get_namdtime(&self) -> usize;
    fn get_prop_eigs<A>(&self) -> &A
        where A: for<'a> nd::AsArray<'a, f64, nd::Ix1>;     // [time, nbasis]
}
