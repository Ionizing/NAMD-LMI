use std::path::Path;

use shared::{
    ndarray::s,
    Result,
    Array1,
    Array2,
    Array3,
};
use hdf5::File as H5File;

use crate::{
    input::Input,
    hamiltonian::{
        Hamiltonian,
        PropagateMethod,
    },
    constants::{
        IMGUNIT,
        HBAR,
        BOLKEV,
    },
};


#[derive(Clone, Copy, PartialEq, Eq)]
pub enum SHMethod {
    FSSH,
    DISH,
    DCSH,
    GFSH,
}

pub struct SurfaceHopping {
    pub shmethod:    SHMethod,
    pub propmethod:  PropagateMethod,
    pub ntraj:       usize,
    pub lexcitation: bool,
    pub hamil:       Hamiltonian,

    //pub prob:        Array3<f64>,       // [namdtime, nbasis, nbasis]
    pub pops:        Array2<f64>,       // [namdtime, nbasis]
    pub recomb:      Array2<f64>,       // [namdtime, nbasis]
    pub energy:      Array1<f64>,       // [namdtime]
}


impl SurfaceHopping {
    pub fn from_input(
        hamil: Hamiltonian,
        inp:   &Input,
    ) -> Self {
        Self::init_with_hamiltonian(
            hamil,
            inp.propmethod,
            inp.shmethod,
            inp.ntraj,
            inp.lexcitation,
        )
    }


    fn init_with_hamiltonian(
        hamil:       Hamiltonian,
        propmethod:  PropagateMethod,
        shmethod:    SHMethod,
        ntraj:       usize,
        lexcitation: bool,
    ) -> Self {
        let namdtime = hamil.namdtime;
        let nbasis   = hamil.nbasis;

        //let prob   = Array3::<f64>::zeros((namdtime, nbasis, nbasis));
        let pops   = Array2::<f64>::zeros((namdtime, nbasis));
        let recomb = Array2::<f64>::zeros((namdtime, nbasis));
        let energy = Array1::<f64>::zeros(namdtime);

        Self {
            shmethod,
            propmethod,
            ntraj,
            lexcitation,
            hamil,

            //prob,
            pops,
            recomb,
            energy
        }
    }


    pub fn run(&mut self) {
        todo!()
    }


    pub fn save_to_h5<P>(&self, fname: &P) -> Result<()>
    where
        P: AsRef<Path> + ?Sized,
    {
        let f = H5File::create(fname)?;

        let time_idx = (0 .. self.hamil.namdtime)
            .map(|v| v as f64 * self.hamil.dt)
            .collect::<Vec<f64>>();

        f.new_dataset_builder().with_data(&time_idx).create("time")?;
        f.new_dataset_builder().with_data(&self.hamil.prop_eigs).create("prop_energy")?;
        f.new_dataset_builder().with_data(&self.hamil.psi_t.mapv(|v| v.re)).create("psi_t_r")?;
        f.new_dataset_builder().with_data(&self.hamil.psi_t.mapv(|v| v.im)).create("psi_t_i")?;

        f.new_dataset_builder().with_data(&self.energy).create("sh_energy")?;
        f.new_dataset_builder().with_data(&self.pops).create("shpops")?;

        if self.shmethod == SHMethod::DISH {
            f.new_dataset_builder().with_data(&self.recomb).create("dish_recomb")?;
        }

        Ok(())
    }


    fn fssh(&mut self) {
        self.pops.fill(0.00);
        let mut prob = Array3::<f64>::zeros((self.hamil.namdtime, self.hamil.nbasis, self.hamil.nbasis));

        for iion in 0 .. self.hamil.namdtime {
            self.hamil.propagate(iion, self.propmethod);
        }

        for iion in 0 .. self.hamil.namdtime {

        }

        todo!()
    }

    fn fssh_hop_prob(&self, iion: usize, istate: usize) -> Array1<f64> {
        let rtime: usize = (iion + self.hamil.namdinit) % (self.hamil.nsw - 1);

        // |phi(j)|^2
        let rho_jj = self.hamil.psi_t[(iion, istate)].norm_sqr();
        // phi(j).conj() * phi(k)
        let rho_jk = self.hamil.psi_t[(iion, istate)].conj() * self.hamil.psi_t.slice(s![iion, ..]).to_owned();
        // \max[\frac{2*\int_t^{t+\Delta t} Re(\rho_{jk}*H_{jk} / -ihbar) dt}{\rho_{jj}}, 0]
        let mut prob = (IMGUNIT * rho_jk * self.hamil.ham_t.slice(s![iion, istate, ..]) * (2.0 * self.hamil.dt / (HBAR * rho_jj)))
            .mapv(|v| if v.re > 0.0 { v.re } else { 0.0 });

        if !self.lexcitation {
            let thermal_factor = (self.hamil.eig_t[(rtime, istate)] - self.hamil.eig_t.slice(s![rtime, ..]).to_owned())
                .mapv(|v| f64::exp(
                        f64::min(v, 0.0) / (BOLKEV * self.hamil.temperature)
                ));
            prob *= &thermal_factor;
        }

        prob
    }


    fn dish(&mut self) {
        todo!()
    }


    fn dcsh(&mut self) {
        todo!()
    }


    fn gfsh(&mut self) {
        todo!()
    }
}
