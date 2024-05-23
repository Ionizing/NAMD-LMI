use std::path::{
    Path,
    PathBuf,
};

use shared::ndarray as nd;
use shared::Result;

use crate::core::{
    NamdConfig,
    Hamiltonian,
    SurfaceHopping,
};
use crate::hamil::{
    config::HamilConfig,
    config::PropagateMethod,
    hamil_impl::SPHamiltonian,
};
use crate::surfhop::wavefunction_impl::SPWavefunction;
use crate::surfhop::config::{
    SHMethod,
    SurfhopConfig,
};


pub struct Surfhop {
    shmethod: SHMethod,
    hamil: SPHamiltonian,
    wfn: SPWavefunction,
    outdir: PathBuf,
    lexcitation: bool,

    ntraj: usize,
    inistep: usize,
    namdtime: usize,
    tdpops: nd::Array2<f64>,            // [namdtime, nbasis]
    tdenergy: nd::Array1<f64>,          // [namdtime]
    tdphotons: nd::Array3<usize>,       // [namdtime, nbasis, nbasis]
    tdphonons: nd::Array3<usize>,       // [namdtime, nbasis, nbasis]
}


impl<'a> SurfaceHopping for Surfhop {
    type ConfigType = SurfhopConfig;
    type HamiltonianType = SPHamiltonian;

    fn get_ntraj(&self) -> usize { self.ntraj }
    fn get_lexcitation(&self) -> bool { self.lexcitation }
    fn get_tdpops(&self) -> nd::ArrayView2<f64> { self.tdpops.view() }
    fn get_tdenergy(&self) -> nd::ArrayView1<f64> { self.tdenergy.view() }
    
    fn run(&mut self) -> Result<()> {
        todo!()
    }
    fn from_config(cfg: &Self::ConfigType) -> Result<Vec<Self>> {
        //let hamil = SPHamiltonian::from_h5(cfg.get_hamil_fname())?;
        //let iband = hamil.get_converted_index(cfg.get_iniband(), cfg.get_inispin())?;

        todo!()
    }
    fn save_to_h5<P>(&self, _fname: P) -> Result<()>
    where P: AsRef<Path> {
        todo!()
    }
}
