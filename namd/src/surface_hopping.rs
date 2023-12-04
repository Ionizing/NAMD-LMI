use std::path::Path;
use std::fmt;

use shared::{
    ndarray::s,
    Result,
    Array1,
    Array2,
    Array3,
    info,
};
use hdf5::File as H5File;
use rand::{thread_rng, Rng};

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
        EPS,
    },
};


#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SHMethod {
    FSSH,
    DISH,
    DCSH,
    GFSH,
}

impl fmt::Display for SHMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use SHMethod::*;
        write!(f, "{}", match self {
            FSSH => "FSSH",
            DISH => "DISH",
            DCSH => "DCSH",
            GFSH => "GFSH",
        })
    }
}


pub struct SurfaceHopping {
    pub shmethod:    SHMethod,
    pub propmethod:  PropagateMethod,
    pub ntraj:       usize,
    pub lexcitation: bool,
    pub hamil:       Hamiltonian,

    pub ndigit:      usize,
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
            inp.ndigit,
            inp.ntraj,
            inp.lexcitation,
        )
    }


    fn init_with_hamiltonian(
        hamil:       Hamiltonian,
        propmethod:  PropagateMethod,
        shmethod:    SHMethod,
        ndigit:      usize,
        ntraj:       usize,
        lexcitation: bool,
    ) -> Self {
        let namdtime = hamil.namdtime;
        let nbasis   = hamil.nbasis;

        let pops   = Array2::<f64>::zeros((namdtime, nbasis));
        let recomb = Array2::<f64>::zeros((namdtime, nbasis));
        let energy = Array1::<f64>::zeros(namdtime);

        Self {
            shmethod,
            propmethod,
            ntraj,
            lexcitation,
            hamil,

            ndigit,
            pops,
            recomb,
            energy
        }
    }


    pub fn run(&mut self) {
        info!("Running surface hopping with namdini = {} ...", self.hamil.namdinit);
        match self.shmethod {
            SHMethod::FSSH => self.fssh(),
            SHMethod::DISH => self.dish(),
            SHMethod::DCSH => self.dcsh(),
            SHMethod::GFSH => self.gfsh(),
        }

        let ndigit = self.ndigit;
        let fname = format!("result_{:0ndigit$}.h5", self.hamil.namdinit);
        self.save_to_h5(&fname).unwrap();
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
            for istate in 0 .. self.hamil.nbasis {
                prob.slice_mut(s![iion, istate, ..]).assign(&self.fssh_hop_prob(iion, istate));
            }
        }

        let mut rng = thread_rng();

        for _ in 0 .. self.ntraj {
            let mut curstate = self.hamil.basisini;
            for iion in 0 .. self.hamil.namdtime {
                let randnum: f64 = rng.gen();
                let hop_dest = prob.slice(s![iion, curstate, ..])
                    .as_slice().unwrap()
                    .partition_point(|v| v < &randnum);
                //info!("curstate = {}, hop_dest = {}, randnum = {}", curstate, hop_dest, randnum);
                curstate = if hop_dest < self.hamil.nbasis { hop_dest } else { curstate };
                self.pops[(iion, curstate)] += 1.0;
            }
        }

        self.pops /= self.ntraj as f64;
    }

    fn fssh_hop_prob(&self, iion: usize, istate: usize) -> Array1<f64> {
        let rtime: usize = (iion + self.hamil.namdinit) % (self.hamil.nsw - 1);

        // |phi(j)|^2
        let rho_jj   = self.hamil.psi_t[(iion, istate)].norm_sqr();
        // phi(j).conj() * phi(k)
        let rho_jk   = self.hamil.psi_t[(iion, istate)].conj() * self.hamil.psi_t.slice(s![iion, ..]).to_owned();
        // \max[\frac{2*\int_t^{t+\Delta t} Re(\rho_{jk}*H_{jk} / -ihbar) dt}{\rho_{jj}}, 0]
        let mut prob = (IMGUNIT * rho_jk * self.hamil.ham_t.slice(s![iion, istate, ..]) * (2.0 * self.hamil.dt / (HBAR * rho_jj)))
            .mapv(|v| if v.re > 0.0 { v.re } else { 0.0 });

        // determine if the electric field is still present
        let has_efield = self.hamil.get_efield(iion).map(|v| v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).unwrap_or(0.0) > EPS;

        // balance factor to restrict upward hops
        // upward hops is restricted only when
        //     - not excitation process
        //     - |EFIELD| != 0
        if !self.lexcitation && !has_efield {
            let thermal_factor = (self.hamil.eig_t[(rtime, istate)] - self.hamil.eig_t.slice(s![rtime, ..]).to_owned())
                .mapv(|v| f64::exp(
                        f64::min(v, 0.0) / (BOLKEV * self.hamil.temperature)
                ));
            prob *= &thermal_factor;
        }

        // cumulative sum
        prob.into_iter()
            .scan(0.0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect()
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
