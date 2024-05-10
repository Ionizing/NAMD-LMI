use std::path::Path;
use std::fmt;

use shared::ndarray as nd;
use shared::{c64, Result};

use crate::core::couplings::Couplings;


pub trait Hamiltonian {
    fn get_nbasis(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_basisini(&self) -> usize;
    fn get_namdinit(&self) -> usize;
    fn get_namdtime(&self) -> usize;
    fn get_nsw(&self) -> usize;
    fn get_nelm(&self) -> usize;
    fn get_temperature(&self) -> f64;

    fn get_tdhamil(&self) -> nd::ArrayView3<c64>;  // [time, nbasis, nbasis]
    fn get_tdproj(&self) -> nd::ArrayView4<f64>;   // [time, nbasis. nion, norbit]
    fn get_hamil(&self, iion: usize) -> nd::ArrayView2<c64>;

    fn propagate<M>(&mut self, iion: usize, method: M)
        where M: PartialEq + Eq + Clone + Copy + fmt::Display;

    fn from_h5<P>(fname: P) -> Result<Self>
        where P: AsRef<Path>,
              Self: Sized;
    fn from_coupling<C>(coup: &C) -> Result<Self>
        where C: Couplings,
              Self: Sized;
    fn save_to_h5<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;

    fn get_edt(&self) -> f64 {
        self.get_potim() / self.get_nelm() as f64
    }

}
