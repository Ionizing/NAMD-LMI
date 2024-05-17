use std::path::Path;

use shared::Result;

use crate::core::NamdConfig;

pub trait Couplings {
    type TdCoupType<'x> where Self: 'x;
    type TdPijType<'x> where Self: 'x;
    type TdProjType<'x> where Self: 'x;
    type TdEigsType<'x> where Self: 'x;
    type ConfigType: NamdConfig;

    fn get_nspin(&self) -> usize;
    fn get_nbands(&self) -> usize;
    fn get_ikpoint(&self) -> usize;
    fn get_brange(&self) -> [usize; 2];
    fn get_nsw(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_temperature(&self) -> f64;
    fn get_efermi(&self) -> f64;

    fn get_tdcoup<'a>(&'a self) -> Self::TdCoupType<'a>;
    fn get_tdpij<'a>(&'a self) -> Self::TdPijType<'a>;
    fn get_tdrij<'a>(&'a self) -> Self::TdPijType<'a>;
    fn get_tdproj<'a>(&'a self) -> Self::TdProjType<'a>;
    fn get_tdeigs<'a>(&'a self) -> Self::TdEigsType<'a>;

    fn from_config(cfg: &Self::ConfigType) -> Result<Self>
        where Self: Sized;
    fn from_h5<P>(fname: P) -> Result<Self>
        where P: AsRef<Path>,
              Self: Sized;
    fn save_to_h5<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;

    fn get_nbrange(&self) -> usize {
        let range = self.get_brange();
        range[1] - range[0] + 1
    }
}
