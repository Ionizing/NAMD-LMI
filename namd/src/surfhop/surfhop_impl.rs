use std::path::{
    Path,
    PathBuf,
};

use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::Result;

use crate::core::{
    NamdConfig,
    Hamiltonian,
    SurfaceHopping,
    Wavefunction,
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
        let shmethod = cfg.get_shmethod();
        let hamil = SPHamiltonian::from_h5(cfg.get_hamil_fname())?;
        // wfn constructed inside the closure
        let outdir = cfg.get_outdir().clone();
        let lexcitation = cfg.get_lexcitation();
        let ntraj = cfg.get_ntraj();
        // inistep got in the closure
        let namdtime = cfg.get_namdtime();
        let nbasis = hamil.get_nbasis();
        let tdpops = nd::Array2::<f64>::zeros((namdtime, nbasis));
        let tdenergy = nd::Array1::<f64>::zeros(namdtime);
        let tdphotons = nd::Array3::<usize>::zeros((namdtime, nbasis, nbasis));
        let tdphonons = nd::Array3::<usize>::zeros((namdtime, nbasis, nbasis));

        cfg.get_inisteps().iter()
            .map(|&istep| -> Result<Self> {
                let hamil = hamil.clone();
                let inistep = istep;
                let wfn = SPWavefunction::from_hamil_and_params(
                    &hamil, cfg.get_iniband(), cfg.get_inispin(),
                    cfg.get_namdtime(), cfg.get_nelm(), istep)?;

                Ok(Self {
                    shmethod,
                    hamil,
                    wfn,
                    outdir: outdir.clone(),
                    lexcitation,

                    ntraj,
                    inistep,
                    namdtime,
                    tdpops: tdpops.clone(),
                    tdenergy: tdenergy.clone(),
                    tdphotons: tdphotons.clone(),
                    tdphonons: tdphonons.clone(),
                })
            })
            .collect()
    }

    fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        let f = H5File::create(fname)?;

        f.new_dataset::<usize>().create("inistep")?.write_scalar(&self.inistep)?;
        f.new_dataset::<usize>().create("namdtime")?.write_scalar(&self.namdtime)?;
        f.new_dataset::<usize>().create("ntraj")?.write_scalar(&self.ntraj)?;
        f.new_dataset::<bool>().create("lexcitation")?.write_scalar(&self.lexcitation)?;

        f.new_dataset_builder().with_data(&self.shmethod.to_string()).create("shmethod")?;
        f.new_dataset_builder().with_data(&self.wfn.get_psi_t().mapv(|v| v.re)).create("psi_t_r")?;
        f.new_dataset_builder().with_data(&self.wfn.get_psi_t().mapv(|v| v.im)).create("psi_t_i")?;
        f.new_dataset_builder().with_data(&self.wfn.get_prop_eigs()).create("prop_eigs_t")?;

        f.new_dataset_builder().with_data(&self.tdpops).create("sh_pops_t")?;
        f.new_dataset_builder().with_data(&self.tdenergy).create("sh_eigs_t")?;
        f.new_dataset_builder().with_data(&self.tdphotons).create("sh_photons_t")?;
        f.new_dataset_builder().with_data(&self.tdphonons).create("sh_phonons_t")?;

        Ok(())
    }
}
