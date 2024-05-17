use std::path::Path;

use shared::ndarray as nd;
use shared::{c64, Result};

use crate::core::Couplings;
use crate::core::NamdConfig;


pub trait Hamiltonian {
    type ConfigType: NamdConfig;
    type CouplingType: Couplings;

    fn get_nbasis(&self) -> usize;
    fn get_nsw(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_temperature(&self) -> f64;

    fn get_hamil(&self, iion: usize) -> nd::ArrayView2<c64>;

    fn from_config(fname: &Self::ConfigType) -> Result<Self>
        where Self: Sized;
    fn from_h5<P>(fname: P) -> Result<Self>
        where P: AsRef<Path>,
              Self: Sized;
    fn save_to_h5<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;
}
