use std::path::Path;

use shared::ndarray as nd;
use shared::Result;

use crate::core::hamiltonian::Hamiltonian;

pub trait SurfaceHopping {
    fn get_ntraj(&self) -> usize;
    fn get_lexcitation(&self) -> usize;
    fn get_tdpops(&self) -> nd::ArrayView2<f64>;
    fn get_tdrecomb(&self) -> nd::ArrayView2<f64>;
    fn get_tdenergy(&self) -> nd::ArrayView1<f64>;

    fn run(&mut self) -> Result<()>;

    fn from_h5<P>(fname: P) -> Result<Self>
        where P: AsRef<Path>,
              Self: Sized;
    fn from_hamil<H>(hamil: &H) -> Result<Self>
        where H: Hamiltonian,
              Self: Sized;
    fn save_to_h5<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;
}
