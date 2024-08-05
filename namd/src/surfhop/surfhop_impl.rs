use std::path::{
    Path,
    PathBuf,
};

use rand::{Rng,thread_rng};
use hdf5::File as H5File;
use shared::ndarray as nd;
use shared::log;
use shared::Result;
use shared::link_file_to;
use rayon::prelude::*;

use crate::core::{
    constants::*,
    Hamiltonian,
    SurfaceHopping,
    Wavefunction,
};
use crate::hamil::{
    Efield,
    SPHamiltonian,
};
use crate::surfhop::SPWavefunction;
use crate::surfhop::config::{
    SHMethod,
    DetailedBalance as DB,
    SurfhopConfig,
};


pub struct Surfhop {
    shmethod: SHMethod,
    hamil: SPHamiltonian,
    wfn: SPWavefunction,
    outdir: PathBuf,
    detailed_balance: DB,

    basis_list: Vec<i32>,
    basis_labels: Option<Vec<String>>,

    ntraj: usize,
    namdinit: usize,
    namdtime: usize,
    time: nd::Array1<f64>,
    tdpops: nd::Array2<f64>,            // [namdtime, nbasis]
    tdenergy: nd::Array1<f64>,          // [namdtime]
    tdphotons: nd::Array3<f64>,         // [namdtime, nbasis, nbasis], emit => plus, absorb => minus
    tdphonons: nd::Array3<f64>,         // [namdtime, nbasis, nbasis]
}


impl<'a> SurfaceHopping for Surfhop {
    type ConfigType = SurfhopConfig;
    type HamiltonianType = SPHamiltonian;

    fn get_ntraj(&self) -> usize { self.ntraj }
    fn get_tdpops(&self) -> nd::ArrayView2<f64> { self.tdpops.view() }
    fn get_tdenergy(&self) -> nd::ArrayView1<f64> { self.tdenergy.view() }
    
    fn run(&mut self) -> Result<()> {
        use SHMethod::*;
        let namdinit = self.namdinit;
        let ndigit = self.hamil.get_ndigit();

        log::info!("Running surface hopping with namdinit = {} ...", namdinit);

        match self.shmethod {
            FSSH => self.fssh(),
            DISH => self.dish(),
            DCSH => self.dcsh(),
        }
        let fname = self.outdir.join(format!("result_{:0ndigit$}.h5", namdinit));
        self.save_to_h5(&fname)
    }

    fn from_config(cfg: &Self::ConfigType) -> Result<Vec<Self>> {
        let shmethod = cfg.get_shmethod();
        let hamil = SPHamiltonian::from_h5(cfg.get_hamil_fname())?;
        // wfn constructed inside the closure
        let outdir = cfg.get_outdir().clone();
        let detailed_balance = cfg.get_detailed_balance();

        let basis_list = hamil.get_basis_list().to_owned();
        let basis_labels = hamil.get_basis_labels().map(|v| v.clone());
        let ntraj = cfg.get_ntraj();
        // namdinit got in the closure
        let namdtime = cfg.get_namdtime();
        let nbasis = hamil.get_nbasis();
        let tdpops = nd::Array2::<f64>::zeros((namdtime, nbasis));
        let tdenergy = nd::Array1::<f64>::zeros(namdtime);
        let tdphotons = nd::Array3::<f64>::zeros((namdtime, nbasis, nbasis));
        let tdphonons = nd::Array3::<f64>::zeros((namdtime, nbasis, nbasis));

        log::info!("Linking Hamiltonian file {:?} to {:?}", cfg.get_hamil_fname(), &outdir);
        link_file_to(cfg.get_hamil_fname(), &outdir)?;

        let potim = hamil.get_potim();
        let nelm = cfg.get_nelm();
        let time: nd::Array1<f64> = (0 .. namdtime).map(|i| i as f64 * potim).collect();
        let efield = hamil.get_efield().map(|x| Efield::singleton_from_str(x).unwrap());
        let eafield = efield.map(|x| {
            let mut e = x.write().unwrap();
            e.print_eafield_tofile(&outdir, namdtime, potim, nelm).unwrap();
            let (_tt, arr) = e.get_eafield_array(namdtime, potim, nelm);
            arr
        });

        cfg.get_inisteps().par_iter()
            .map(|&istep| -> Result<Self> {
                let hamil = hamil.clone();
                let namdinit = istep;
                let wfn = SPWavefunction::from_hamil_and_params(
                    &hamil, cfg.get_iniband(),
                    cfg.get_namdtime(), cfg.get_nelm(), istep, eafield.clone())?;

                Ok(Self {
                    shmethod,
                    hamil,
                    wfn,
                    outdir: outdir.clone(),
                    detailed_balance,

                    basis_list: basis_list.clone(),
                    basis_labels: basis_labels.clone(),

                    ntraj,
                    namdinit,
                    namdtime,
                    time: time.clone(),
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

        f.new_dataset::<usize>().create("namdinit")?.write_scalar(&self.namdinit)?;
        f.new_dataset::<usize>().create("namdtime")?.write_scalar(&self.namdtime)?;
        f.new_dataset::<usize>().create("ntraj")?.write_scalar(&self.ntraj)?;

        f.new_dataset_builder().with_data(&self.basis_list).create("basis_list")?;
        if let Some(labels) = self.basis_labels.as_ref() {
            let data = labels.join("\n");
            f.new_dataset_builder().with_data(&data).create("basis_labels")?;
        }

        f.new_dataset_builder().with_data(&self.detailed_balance.to_string()).create("detailed_balance")?;
        f.new_dataset_builder().with_data(&self.time).create("time")?;
        f.new_dataset_builder().with_data(&self.shmethod.to_string()).create("shmethod")?;
        f.new_dataset_builder().with_data(&self.wfn.get_psi_t().mapv(|v| v.re)).create("psi_t_r")?;
        f.new_dataset_builder().with_data(&self.wfn.get_psi_t().mapv(|v| v.im)).create("psi_t_i")?;
        f.new_dataset_builder().with_data(&self.wfn.get_prop_eigs()).create("prop_energy")?;

        f.new_dataset_builder().with_data(&self.tdpops).create("sh_pops")?;
        f.new_dataset_builder().with_data(&self.tdenergy).create("sh_energy")?;
        f.new_dataset_builder().with_data(&self.tdphotons).create("sh_photons_t")?;
        f.new_dataset_builder().with_data(&self.tdphonons).create("sh_phonons_t")?;

        Ok(())
    }
}


impl Surfhop {
    pub fn get_detailed_balance(&self) -> DB { self.detailed_balance }

    fn fssh(&mut self) {
        let namdtime = self.namdtime;
        let namdinit = self.namdinit;
        let nbasis = self.hamil.get_nbasis();

        self.wfn.propagate_full(&self.hamil);

        let mut cumprob = nd::Array3::<f64>::zeros((namdtime, nbasis, nbasis));
        let mut epc_normsqr = nd::Array3::<f64>::zeros((namdtime, nbasis, nbasis));
        let mut lmi_normsqr = nd::Array3::<f64>::zeros((namdtime, nbasis, nbasis));
        let mut td_eigs = nd::Array2::<f64>::zeros((namdtime, nbasis));
        for iion in 0 .. namdtime {
            for istate in 0 .. nbasis {
                cumprob.slice_mut(nd::s![iion, istate, ..])
                    .assign(&self.fssh_hop_prob(iion, istate));
            }
            epc_normsqr.slice_mut(nd::s![iion, .., ..])
                .assign(&self.hamil.get_hamil0_rtime(iion, namdinit).mapv(|v| v.norm_sqr()));
            lmi_normsqr.slice_mut(nd::s![iion, .., ..])
                .assign(&self.wfn.get_lmi(&self.hamil, iion, 0).mapv(|v| v.norm_sqr()));
            td_eigs.slice_mut(nd::s![iion, ..])
                .assign(&self.hamil.get_eigs_rtime(iion, namdinit));
        }

        // cancel mutability
        let cumprob = cumprob;
        let epc_normsqr = epc_normsqr;
        let lmi_normsqr = lmi_normsqr;
        let td_eigs = td_eigs;


        let mut rng = thread_rng();
        self.tdpops.fill(0.0);
        for _ in 0 .. self.ntraj {
            let mut curstate: usize = self.wfn.get_basisini();
            let mut nxtstate: usize;
            for iion in 0 .. namdtime {
                let randnum: f64 = rng.gen();
                let hop_dest = cumprob.slice(nd::s![iion, curstate, ..])
                    .as_slice().unwrap()
                    .partition_point(|v| v < &randnum);
                nxtstate = if hop_dest < nbasis { hop_dest } else { curstate };

                // when hopping happens
                if curstate != nxtstate {
                    let mut epc: f64 = epc_normsqr[(iion, curstate, nxtstate)];
                    let     lmi: f64 = lmi_normsqr[(iion, curstate, nxtstate)];

                    // normalize
                    let scale = 1.0 / (epc + lmi);
                    epc *= scale;
                    //lmi *= scale;

                    let randnum2: f64 = rng.gen();

                    // downward hop: cur > nxt => tdxxx[cur, nxt] +1
                    if td_eigs[(iion, curstate)] > td_eigs[(iion, nxtstate)] {
                        if randnum2 < epc {     // phonon emitted
                            self.tdphonons[(iion, curstate, nxtstate)] += 1.0;
                        } else {                // photon emitted
                            self.tdphotons[(iion, curstate, nxtstate)] += 1.0;
                        }
                    // upward hop: cur < nxt => tdxxx[cur, nxt] -1
                    } else {
                        if randnum2 < epc {     // phonon absorbed
                            self.tdphonons[(iion, curstate, nxtstate)] -= 1.0;
                        } else {                // photon absorbed
                            self.tdphotons[(iion, curstate, nxtstate)] -= 1.0;
                        }
                    }
                }

                curstate = nxtstate;
                self.tdpops[(iion, curstate)] += 1.0;
            }
        }

        self.tdphotons /= self.ntraj as f64;
        self.tdphonons /= self.ntraj as f64;
        self.tdpops /= self.ntraj as f64;

        for iion in 0 .. namdtime {
            self.tdenergy[iion] = (
                self.tdpops.slice(nd::s![iion, ..]).to_owned() *
                self.hamil.get_eigs_rtime(iion, namdinit)
            ).sum();
        }
    }

    fn fssh_hop_prob(&self, iion: usize, istate: usize) -> nd::Array1<f64> {
        let namdinit = self.wfn.get_namdinit();

        // |phi(j)|^2
        let rho_jj = self.wfn.get_pop_t()[(iion, istate)];

        // phi(j).conj() * phi(k)
        let rho_jk = self.wfn.get_psi_t()[(iion, istate)].conj() * self.wfn.get_psi(iion).to_owned();

        // -i hbar <phi(j) | d/dt | phi(k)>
        let epc = self.hamil.get_hamil0_rtime(iion, namdinit);       // view

        // e A <phi(j) | p | phi(k)>/m
        let lmi = self.wfn.get_lmi(&self.hamil, iion, 0);            // owned

        // total hamiltonian
        let ham = lmi + epc;

        // hopping probability
        let mut prob = (
            rho_jk * ham.slice(nd::s![istate, ..]) * (-2.0 * self.hamil.get_potim() / (HBAR * rho_jj))
        ).mapv(|v| v.im.max(0.0));

        // determine if the electric field is still present
        let has_efield = self.wfn.get_eafield_array()
            .map(|[e, _a]| {
                let idx = iion * self.wfn.get_nelm();
                let v = e[idx];
                v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
            })
            .unwrap_or(0.0) > EPS;

        // balance factor to restrict upward hops
        // upward hops is restricted only when
        //     - not excitation process
        //     - |EFIELD| != 0
        if DB::Always == self.detailed_balance ||
            (DB::DependsOnEField == self.detailed_balance && !has_efield) {
            let eig = self.hamil.get_eigs_rtime(iion, namdinit);
            let thermal_factor = (eig[istate] - eig.to_owned())
                .mapv(|v| f64::exp(
                        f64::min(v, 0.0) / (BOLKEV * self.hamil.get_temperature())
                ));
            prob *= &thermal_factor;
        }

        // cumulative sum
        prob.into_iter()
            .scan(0.0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect::<nd::Array1<f64>>()
    }


    fn dish(&mut self) {
        todo!()
    }

    fn dcsh(&mut self) {
        todo!()
    }
}
