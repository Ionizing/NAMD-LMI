use std::path::Path;

use shared::Result;

use crate::core::NamdConfig;

pub trait Couplings {
    type TdCoupType<'a>;
    type TdPijType<'a>;
    type TdProjType<'a>;
    type TdEigsType<'a>;
    type ConfigType<'a>: NamdConfig<'a>;

    fn get_nspin(&self) -> usize;
    fn get_nbands(&self) -> usize;
    fn get_ikpoint(&self) -> usize;
    fn get_brange(&self) -> [usize; 2];
    fn get_nsw(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_efermi(&self) -> f64;

    fn get_tdcoup<'a>(&self) -> Self::TdCoupType<'a>;
    fn get_tdpij<'a>(&self) -> Self::TdPijType<'a>;
    fn get_tdrij<'a>(&self) -> Self::TdPijType<'a>;
    fn get_tdproj<'a>(&self) -> Self::TdProjType<'a>;
    fn get_tdeigs<'a>(&self) -> Self::TdEigsType<'a>;

    fn from_config<'a>(cfg: &Self::ConfigType<'a>) -> Result<Self>
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
