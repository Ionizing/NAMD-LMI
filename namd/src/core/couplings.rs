use std::path::Path;

use shared::Result;
use shared::ndarray as nd;
use shared::c64;

pub trait Couplings {
    fn get_nspin(&self) -> usize;
    fn get_nbands(&self) -> usize;
    fn get_nkpoint(&self) -> usize;
    fn get_brange(&self) -> [usize; 2];
    fn get_nsw(&self) -> usize;
    fn get_potim(&self) -> f64;
    fn get_efermi(&self) -> f64;

    fn get_tdcoup<D>(&self) -> nd::ArrayView<'_, c64, D>;
    fn get_tdpij<D>(&self) -> nd::ArrayView<'_, c64, D>;
    fn get_tdrij<D>(&self) -> nd::ArrayView<'_, c64, D>;
    fn get_tdproj<D>(&self) -> nd::ArrayView<'_, f64, D>;
    fn get_tdeigs<D>(&self) -> nd::ArrayView<'_, f64, D>;

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
