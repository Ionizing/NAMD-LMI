use std::path::Path;

use shared::Result;

use crate::core::NamdConfig;

pub trait Couplings<'a> {
    type TdCoupType;
    type TdPijType;
    type TdProjType;
    type TdEigsType;
    type ConfigType: NamdConfig<'a>;

    fn get_nspin(&self) -> usize;
    fn get_nbands(&self) -> usize;
    fn get_ikpoint(&self) -> usize;
    fn get_brange(&self) -> [usize; 2];
    fn get_nsw(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_efermi(&self) -> f64;

    fn get_tdcoup(&self) -> Self::TdCoupType;
    fn get_tdpij(&self) -> Self::TdPijType;
    fn get_tdrij(&self) -> Self::TdPijType;
    fn get_tdproj(&self) -> Self::TdProjType;
    fn get_tdeigs(&self) -> Self::TdEigsType;

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
