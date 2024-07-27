use std::path::Path;

use shared::ndarray as nd;
use shared::Result;

use crate::core::NamdConfig;
use crate::core::hamiltonian::Hamiltonian;

pub trait SurfaceHopping {
    type ConfigType: NamdConfig;
    type HamiltonianType: Hamiltonian;

    fn get_ntraj(&self) -> usize;
    fn get_tdpops(&self) -> nd::ArrayView2<f64>;
    fn get_tdenergy(&self) -> nd::ArrayView1<f64>;

    fn run(&mut self) -> Result<()>;

    fn from_config(cfg: &Self::ConfigType) -> Result<Vec<Self>>
        where Self: Sized;
    fn save_to_h5<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;
}
