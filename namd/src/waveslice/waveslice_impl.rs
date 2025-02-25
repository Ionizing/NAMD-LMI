use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};

use hdf5::File as H5File;
use rayon::prelude::*;
use shared::{
    c64,
    ndarray as nd,
    Context,
    Result,
    anyhow,
    log,
};

use vasp_parsers::{
    procar::Procar,
    Wavecar,
    WavecarType,
    Xdatcar,
    Poscar,
};
use crate::waveslice::WavesliceConfig;


struct SliceIRet {
    eigs_i:     nd::Array3<f64>, // [nspin, nkpoints, nbrange]
    fweights_i: nd::Array3<f64>, // [nspin, nkpoints, nbrange]
    coeffs_i:  nd::Array4<c64>, // [nspin, nkpoints, nbrange, nplwmax]
    projs_i:   nd::Array5<f64>, // [nsw, nkpoints, nions, nspinor, nbrange]
    poscar:    Poscar,          // parsed from CONTCAR
    efermi:    f64,
}


struct SliceTotRet {
    eigs:     nd::Array4<f64>,  // [nsw, nspin, nkpoints, nbrange]
    fweights: nd::Array4<f64>,  // [nsw, nspin, nkpoints, nbrange]
    coeffs:   nd::Array5<c64>,  // [nsw, nspin, nkpoints, nbrange, nplwmax]
    projs:    nd::Array6<f64>,  // [nsw, nkpoints, nspinor, nbrange, nions, nspd]
    efermis:  nd::Array1<f64>,  // [nsw,]
    xdatcar:  Xdatcar,          // [nsw,]
}


pub struct Waveslice {
    /// Should count from 1
    ikpoints: Vec<usize>,

    /// number of spin channels
    nspin:  usize,

    /// check if use spin diabatics representation
    spin_diabatics: bool,

    /// check if this waveslice contains normalcar
    lnormalcar: bool,

    /// check if this waveslice contains soccar
    lsoccar: bool,

    /// Total bands of WAVECAR
    nbands: usize,

    /// stores `(brange[0] ..= brange[1])` where `brange[1]` is included, counts from 1
    brange:  [usize; 2],

    /// brange[1] - brange[0] + 1
    nbrange: usize,

    /// Number of digits for each step index.
    ///
    /// ndigit("0001") == 4, ndigit("00001") == 5
    ndigit: usize,

    /// Total number of WAVECARs
    nsw:     usize,

    /// ENCUT from INCAR
    encut:   f64,
    
    phasecorrection:   bool,
    unitary_transform: bool,
    rearrangement:     bool,

    /// Type of WAVECAR, should be "std", "gamx", "gamz" or "ncl".
    wavetype: String,

    /// Lattice vectors in real space, [3, 3]
    real_cell: nd::Array2<f64>,

    /// Lattice vectors in reciprocal space, [3, 3]
    reci_cell: nd::Array2<f64>,

    /// Number of grids in reciprocal space, [3,]
    ngrid:     [usize; 3],

    /// Fermi level of each WAVECAR
    efermis: nd::Array1<f64>,

    /// K-vector of selected K points, [nkpoints, 3]
    kvecs:    nd::Array2<f64>,

    /// Number of plane wave coefficients for selected K points, [nkpoints,]
    num_plws: Vec<usize>,

    /// G vectors of each plane wave for selected K points, [nkpoints, [nplw, 3]]
    gvecs:    Vec<nd::Array2<i64>>,

    /// Eigenvalue of each band for each Waveslice, [nkpoints, [nspin, nsw, nbrange]]
    eigs:     Vec<nd::Array3<f64>>,

    /// Fermi occupation of each band for each Waveslice, [nkpoints, [nspin, nsw, nbrange]]
    fweights: Vec<nd::Array3<f64>>,

    /// Coefficients for each k-point, [nkpoints, [nsw, nspin, nbrange, nplwmax]]
    coeffs:   Vec<nd::Array4<c64>>,

    /// Slice of PROCAR for each WAVECAR, [nkpoints, [nsw, nspinor, nbrange, nions, nspd]]
    projs:    Vec<nd::Array5<f64>>,

    /// Slice of NormalCAR, projection coefficient of PAW projectors, [nkpoints, [nsw, 2, nbands, nproj]]
    cprojs:   Option<Vec<nd::Array4<c64>>>,

    /// SocCars, [nsw, 4, nproj, nproj]
    soccars:  Option<nd::Array4<c64>>,

    /// Spin orbit matrix, [nsw, 4, nbrange, nbrange]
    hmms: Option<nd::Array4<c64>>,

    /// Atomic trajectory of each step, [nsw]
    ///
    /// In .h5 , it will be converted to string with XDATCAR format.
    xdatcar:  Xdatcar,
}


impl Waveslice {
    pub fn get_ikpoints(&self) -> &[usize] { &self.ikpoints }
    pub fn get_nspin(&self) -> usize { self.nspin }
    pub fn get_spin_diabatics(&self) -> bool { self.spin_diabatics }
    pub fn get_lnormalcar(&self) -> bool { self.lnormalcar }
    pub fn get_lsoccar(&self) -> bool { self.lsoccar }
    pub fn get_nbands(&self) -> usize { self.nbands }
    pub fn get_brange(&self) -> [usize;2] { self.brange }
    pub fn get_nbrange(&self) -> usize { self.nbrange }
    pub fn get_ndigit(&self) -> usize { self.ndigit }
    pub fn get_nsw(&self) -> usize { self.nsw }
    pub fn get_encut(&self) -> f64 { self.encut }
    pub fn get_phasecorrection(&self) -> bool { self.phasecorrection }
    pub fn get_unitary_transform(&self) -> bool { self.unitary_transform }
    pub fn get_rearrangement(&self) -> bool { self.rearrangement }
    pub fn get_wavetype(&self) -> &str { &self.wavetype }
    pub fn get_real_cell(&self) -> nd::ArrayView2<f64> { self.real_cell.view() }
    pub fn get_reci_cell(&self) -> nd::ArrayView2<f64> { self.reci_cell.view() }
    pub fn get_ngrid(&self) -> [usize; 3] { self.ngrid }
    pub fn get_efermis(&self) -> nd::ArrayView1<f64> { self.efermis.view() }
    pub fn get_kvecs(&self) -> nd::ArrayView2<f64> { self.kvecs.view() }
    pub fn get_num_plws(&self) -> &[usize] { &self.num_plws }
    pub fn get_gvecs(&self) -> &[nd::Array2<i64>] { &self.gvecs }
    pub fn get_eigs(&self) -> &[nd::Array3<f64>] { self.eigs.as_ref() }
    pub fn get_fweights(&self) -> &[nd::Array3<f64>] { self.fweights.as_ref() }
    pub fn get_coeffs(&self) -> &[nd::Array4<c64>] { self.coeffs.as_ref() }
    pub fn get_projs(&self) -> &[nd::Array5<f64>] { self.projs.as_ref() }
    pub fn get_cprojs(&self) -> Option<&[nd::Array4<c64>]> { self.cprojs.as_ref().map(|x| x.as_ref()) }
    pub fn get_soccars(&self) -> Option<nd::ArrayView4<c64>> { self.soccars.as_ref().map(|x| x.view()) }
    pub fn get_hmms(&self) -> Option<nd::ArrayView4<c64>> { self.hmms.as_ref().map(|x| x.view()) }
    pub fn get_xdatcar(&self) -> &Xdatcar { &self.xdatcar }


    pub fn from_config(cfg: &WavesliceConfig) -> Result<Self> {
        if cfg.get_waveslicefname().is_file() {
            log::info!("Found pre-calculated NAC available in {:?}, reading WaveSlice from it ...",
                cfg.get_waveslicefname());
            let waveslice = Self::from_h5(cfg.get_waveslicefname())?;

            let ikpoints = cfg.get_ikpoints();
            if ikpoints != &[0] {
                anyhow::ensure!(ikpoints == waveslice.ikpoints);
            } else {
                anyhow::ensure!(ikpoints.iter().last().unwrap() + 1 == ikpoints.len());
            }

            anyhow::ensure!(cfg.get_nsw() == waveslice.nsw, "Incompatible nsw.");
            anyhow::ensure!(cfg.get_ndigit() == waveslice.ndigit, "Incompatible ndigit.");
            anyhow::ensure!(cfg.get_brange() == &waveslice.brange, "Incompatible brange.");
            anyhow::ensure!(cfg.get_phasecorrection() == waveslice.phasecorrection, "Incompatible phasecorrection.");
            anyhow::ensure!(cfg.get_unitary_transform() == waveslice.unitary_transform, "Incompatible unitary_transform.");
            anyhow::ensure!(cfg.get_rearrangement() == waveslice.rearrangement, "Incompatible rearrangement.");

            return Ok(waveslice);
        }

        Self::calculate_from_scratch(&cfg)
    }


    pub fn from_h5<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        let f = H5File::open(fname)?;

        let ikpoints: Vec<usize> = f.dataset("ikpoints")?.read_raw()?;

        let nspin  = f.dataset("nspin")?.read_scalar::<usize>()?;
        let spin_diabatics = f.dataset("spin_diabatics")?.read_scalar::<bool>()?;
        let lnormalcar = f.dataset("lnormalcar")?.read_scalar::<bool>()?;
        let lsoccar = f.dataset("lsoccar")?.read_scalar::<bool>()?;
        let nbands = f.dataset("nbands")?.read_scalar::<usize>()?;
        let brange = f.dataset("brange")?.read_scalar::<[usize;2]>()?;
        let nbrange = f.dataset("nbrange")?.read_scalar::<usize>()?;
        let ndigit = f.dataset("ndigit")?.read_scalar::<usize>()?;
        let nsw    = f.dataset("nsw")?.read_scalar::<usize>()?;
        let encut  = f.dataset("encut")?.read_scalar::<f64>()?;

        let phasecorrection = f.dataset("phasecorrectio")?.read_scalar::<bool>()?;
        let unitary_transform = f.dataset("unitary_transform")?.read_scalar::<bool>()?;
        let rearrangement = f.dataset("rearrangement")?.read_scalar::<bool>()?;

        let wavetype = {
            let raw = f.dataset("wavetype")?.read_raw::<u8>()?;
            String::from_utf8(raw)?
        };
        let real_cell: nd::Array2<f64> = f.dataset("real_cell")?.read()?;
        let reci_cell: nd::Array2<f64> = f.dataset("reci_cell")?.read()?;

        let ngrid = f.dataset("ngrid")?.read_scalar::<[usize;3]>()?;

        let efermis: nd::Array1<f64> = f.dataset("efermis")?.read()?;
        let kvecs: nd::Array2<f64>   = f.dataset("kvecs")?.read()?;
        let num_plws: Vec<usize>     = f.dataset("num_plws")?.read_raw()?;

        let nkpoints = ikpoints.len();
        let mut gvecs = Vec::<nd::Array2<i64>>::new();
        let mut eigs = Vec::<nd::Array3<f64>>::new();
        let mut fweights = Vec::<nd::Array3<f64>>::new();
        let mut coeffs = Vec::<nd::Array4<c64>>::new();
        let mut projs = Vec::<nd::Array5<f64>>::new();

        let mut cprojs = if lnormalcar {
            Some(Vec::<nd::Array4<c64>>::new())
        } else { None };

        for (i, &ik) in ikpoints.iter().enumerate() {
            let grp = f.group(format!("k{ik}").as_ref())?;

            let gvec: nd::Array2<i64> = grp.dataset("gvecs")?.read()?;
            gvecs.push(gvec);

            let eig: nd::Array3<f64> = grp.dataset("eigs")?.read()?;
            eigs.push(eig);

            let fweight: nd::Array3<f64> = grp.dataset("fweights")?.read()?;
            fweights.push(fweight);

            let coeff_r: nd::Array4<f64> = grp.dataset("coeffs_r")?.read()?;
            let coeff_i: nd::Array4<f64> = grp.dataset("coeffs_i")?.read()?;
            coeffs.push(coeff_r.mapv(|v| c64::new(v, 0.0)) + coeff_i.mapv(|v| c64::new(0.0, v)));

            let proj: nd::Array5<f64> = grp.dataset("projs")?.read()?;
            projs.push(proj);

            if let Some(cproj) = cprojs.as_mut() {
                let cproj_r: nd::Array4<f64> = grp.dataset("cprojs_r")?.read()?;
                let cproj_i: nd::Array4<f64> = grp.dataset("cprojs_i")?.read()?;
                cproj[i] = cproj_r.mapv(|x| c64::new(x, 0.0)) + cproj_i.mapv(|x| c64::new(0.0, x));
            }
        }

        let soccars = if lsoccar {
            let soccar_r: nd::Array4<f64> = f.dataset("soccars_r")?.read()?;
            let soccar_i: nd::Array4<f64> = f.dataset("soccars_i")?.read()?;
            Some(soccar_r.mapv(|x| c64::new(x, 0.0)) + soccar_i.mapv(|x| c64::new(0.0, x)))
        } else {
            None
        };

        let hmms = if spin_diabatics {
            let hmm_r: nd::Array4<f64> = f.dataset("hmms_r")?.read()?;
            let hmm_i: nd::Array4<f64> = f.dataset("hmms_i")?.read()?;
            Some(hmm_r.mapv(|x| c64::new(x, 0.0)) + hmm_i.mapv(|x| c64::new(0.0, x)))
        } else {
            None
        };

        let xdatcar = {
            let bytes = f.dataset("xdatcar")?.read_raw::<u8>()?;
            let xdatcar_str = String::from_utf8(bytes)?;
            Xdatcar::from_txt(&xdatcar_str)?
        };

        Ok(Self {
            ikpoints,
            nspin,
            spin_diabatics,
            lnormalcar,
            lsoccar,
            nbands,
            brange,
            nbrange,
            ndigit,
            nsw,
            encut,
            phasecorrection,
            unitary_transform,
            rearrangement,
            wavetype,
            real_cell,
            reci_cell,
            ngrid,
            efermis,
            kvecs,
            num_plws,
            gvecs,
            eigs,
            fweights,
            coeffs,
            projs,
            cprojs,
            soccars,
            hmms,
            xdatcar,
        })
    }


    pub fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        let f = H5File::create(fname)?;

        f.new_dataset_builder().with_data(&self.ikpoints).create("ikpoints")?;

        f.new_dataset::<usize>().create("nspin")?.write_scalar(&self.nspin)?;
        f.new_dataset::<bool>().create("spin_diabatics")?.write_scalar(&self.spin_diabatics)?;
        f.new_dataset::<bool>().create("lnormalcar")?.write_scalar(&self.lnormalcar)?;
        f.new_dataset::<bool>().create("lsoccar")?.write_scalar(&self.lsoccar)?;
        f.new_dataset::<usize>().create("nbands")?.write_scalar(&self.nbands)?;
        f.new_dataset::<[usize;2]>().create("brange")?.write_scalar(&self.brange)?;
        f.new_dataset::<usize>().create("nbrange")?.write_scalar(&self.nbrange)?;
        f.new_dataset::<usize>().create("ndigit")?.write_scalar(&self.ndigit)?;
        f.new_dataset::<usize>().create("nsw")?.write_scalar(&self.nsw)?;
        f.new_dataset::<f64>().create("encut")?.write_scalar(&self.encut)?;

        f.new_dataset::<bool>().create("phasecorrection")?.write_scalar(&self.phasecorrection)?;
        f.new_dataset::<bool>().create("unitary_transform")?.write_scalar(&self.unitary_transform)?;
        f.new_dataset::<bool>().create("rearrangement")?.write_scalar(&self.rearrangement)?;

        f.new_dataset_builder().with_data(&self.wavetype.as_bytes()).create("wavetype")?;
        f.new_dataset_builder().with_data(&self.real_cell).create("real_cell")?;
        f.new_dataset_builder().with_data(&self.reci_cell).create("reci_cell")?;

        f.new_dataset::<[usize;3]>().create("ngrid")?.write_scalar(&self.ngrid)?;

        f.new_dataset_builder().with_data(&self.efermis).create("efermis")?;
        f.new_dataset_builder().with_data(&self.kvecs).create("kvecs")?;
        f.new_dataset_builder().with_data(&self.num_plws).create("num_plws")?;

        for (i, &ik) in self.ikpoints.iter().enumerate() {
            let grp = f.create_group(format!("k{ik}").as_ref())?;

            grp.new_dataset_builder().with_data(&self.kvecs.slice(nd::s![i, ..])).create("kvec")?;

            grp.new_dataset::<usize>().create("nplw")?.write_scalar(&self.num_plws[i])?;
            grp.new_dataset_builder().with_data(&self.gvecs[i]).create("gvecs")?;

            grp.new_dataset_builder().with_data(&self.eigs[i]).create("eigs")?;
            grp.new_dataset_builder().with_data(&self.fweights[i]).create("fweights")?;

            grp.new_dataset_builder().with_data(&self.coeffs[i].mapv(|x| x.re)).create("coeffs_r")?;
            grp.new_dataset_builder().with_data(&self.coeffs[i].mapv(|x| x.im)).create("coeffs_i")?;

            grp.new_dataset_builder().with_data(&self.projs[i]).create("projs")?;

            if let Some(cprojs) = self.cprojs.as_ref() {
                grp.new_dataset_builder().with_data(&cprojs[i].mapv(|x| x.re)).create("cprojs_r")?;
                grp.new_dataset_builder().with_data(&cprojs[i].mapv(|x| x.im)).create("cprojs_i")?;
            }
        }

        if let Some(soccars) = self.soccars.as_ref() {
            f.new_dataset_builder().with_data(&soccars.mapv(|x| x.re)).create("soccars_r")?;
            f.new_dataset_builder().with_data(&soccars.mapv(|x| x.im)).create("soccars_i")?;
        }

        if let Some(hmm) = self.hmms.as_ref() {
            f.new_dataset_builder().with_data(&hmm.mapv(|x| x.re)).create("hmms_r")?;
            f.new_dataset_builder().with_data(&hmm.mapv(|x| x.im)).create("hmms_i")?;
        }

        let xdatcar_str = format!("{}", self.xdatcar);
        f.new_dataset_builder().with_data(&xdatcar_str.as_bytes()).create("xdatcar")?;

        Ok(())
    }


    pub fn calculate_from_scratch(cfg: &WavesliceConfig) -> Result<Self> {
        log::info!("Slicing WAVECAR from scratch in {:?}/.../WAVECARs ...", cfg.get_rundir());

        let rundir = Path::new(cfg.get_rundir());
        let nsw    = cfg.get_nsw();
        let brange = Range { start: cfg.get_brange()[0] - 1, end: cfg.get_brange()[1] };
        let nbrange = brange.len();
        let ndigit = cfg.get_ndigit();
        

        let path_1 = rundir.join(format!("{:0ndigit$}", 1)).join("WAVECAR");
        let w1 = Wavecar::from_file(&path_1)
            .with_context(|| format!("Failed to parse {:?} as WAVECAR.", &path_1))
            .unwrap();

        let nspin   = w1.nspin as usize;
        let nbands  = w1.nbands as usize;

        let encut = w1.encut;
        let wavetype = match w1.wavecar_type {
            WavecarType::Standard => "std",
            WavecarType::GammaHalf(shared::Axis::X) => "gamx",
            WavecarType::GammaHalf(shared::Axis::Z) => "gamz",
            WavecarType::NonCollinear => "ncl",
            WavecarType::GammaHalf(shared::Axis::Y) => panic!("Impossible branch"),
        }.to_string();

        let nkpoints_wavecar = w1.nkpoints;

        let ikpoints = cfg.get_ikpoints();
        let ikpoints = if ikpoints.len() == 1 && ikpoints[0] == 0 {
            (0 .. nkpoints_wavecar as usize).collect::<Vec<usize>>()
        } else {
            ikpoints.into_iter().map(|x| x - 1).collect::<Vec<usize>>()
        };  // ikpoints now counts from 0
        let nkpoints = ikpoints.len();
        anyhow::ensure!(nkpoints != 0, "At least one kpoint should be included");

        let real_cell = w1.acell.iter().flatten().cloned().collect::<nd::Array1<f64>>().into_shape((3, 3)).unwrap();
        let reci_cell = w1.bcell.iter().flatten().cloned().collect::<nd::Array1<f64>>().into_shape((3, 3)).unwrap();
        let ngrid = [w1.ngrid[0] as usize, w1.ngrid[1] as usize, w1.ngrid[2] as usize];

        let kvecs    = {
            let mut kvecs = nd::Array2::<f64>::zeros((nkpoints, 3));
            for (i, &ik) in ikpoints.iter().enumerate() {
                kvecs.row_mut(i).assign(&w1.kvecs.row(ik));
            }
            kvecs
        };
        let num_plws = ikpoints.iter().map(|ik| w1.nplws[*ik] as usize).collect::<Vec<usize>>();
        let lncl = w1.wavecar_type == WavecarType::NonCollinear;

        let gvecs = ikpoints.iter().enumerate().map(|(i, &ik)| {
            let nspinor = if lncl { 2usize } else { 1 };
            w1.generate_fft_grid((ik) as u64)
                .into_iter().flatten()
                .collect::<nd::Array1::<i64>>()
                .into_shape((num_plws[i]/nspinor, 3))
                .unwrap()
        }).collect::<Vec<nd::Array2<i64>>>();

        let procar_1 = rundir.join(format!("{:0ndigit$}", 1)).join("PROCAR");
        let p1 = Procar::from_file(&procar_1)
            .with_context(|| format!("Failed to parse {:?} as PROCAR.", &procar_1))
            .unwrap();
        let nions = p1.pdos.nions as usize;
        let nspd  = p1.pdos.projected.shape()[4];


        let phasecorrection = cfg.get_phasecorrection();
        let unitary_transform = cfg.get_unitary_transform();
        let rearrangement = cfg.get_rearrangement();


        let SliceTotRet {
            eigs, fweights, coeffs, projs, efermis, xdatcar
        } = Self::from_wavecars(&rundir, nsw, &ikpoints, brange.clone(), ndigit,
            nspin, &num_plws, lncl, nions, nspd)?;

        Ok(Self {
            ikpoints,
            nspin,
            nbands,
            brange: cfg.get_brange().clone(),
            nbrange,
            ndigit,
            nsw,
            encut,

            phasecorrection,
            unitary_transform,
            rearrangement,

            wavetype,
            real_cell,
            reci_cell,
            ngrid,

            efermis,
            kvecs,
            num_plws,
            gvecs,

            eigs,
            fweights,
            coeffs,
            projs,
            xdatcar,
        })
    }


    fn from_wavecars(rundir: &Path, nsw: usize, ikpoints: &[usize], brange: Range<usize>,
        ndigit: usize, nspin: usize, num_plws: &[usize], lncl: bool, nions: usize, nspd: usize,
        ) -> Result<SliceTotRet> {
        
        let nbrange   = brange.clone().count();
        let nkpoints  = ikpoints.len();
        let nplws_max = num_plws.iter().cloned().max().unwrap();
        let nspinors  = if lncl { 4 } else { nspin };

        let ret_eigs = Arc::new(Mutex::new(
                nd::Array4::<f64>::zeros((nsw, nspin, nkpoints, nbrange))
                ));
        let ret_fweights = Arc::new(Mutex::new(
                nd::Array4::<f64>::zeros((nsw, nspin, nkpoints, nbrange))
                ));
        let ret_coeffs = Arc::new(Mutex::new(
                nd::Array5::<c64>::zeros((nsw, nspin, nkpoints, nbrange, nplws_max))
                ));
        let ret_projs = Arc::new(Mutex::new(
                nd::Array6::<f64>::zeros((nsw, nkpoints, nspinors, nbrange, nions, nspd))
                ));
        let ret_efermis = Arc::new(Mutex::new(
                nd::Array1::<f64>::zeros((nsw,))
                ));
        let ret_xdatcar = Arc::new(Mutex::new(vec![Poscar::default(); nsw]));

        
        let remain_count = Arc::new(Mutex::new(nsw - 1));
        (0 .. nsw).into_par_iter().for_each(|isw| {
            let path_i = rundir.join(format!("{:0ndigit$}", isw+1));
            {
                let mut remain_now = remain_count.lock().unwrap();
                let remains: usize = *remain_now;
                log::info!(" Slicing WAVECAR and PROCAR from {:?} ..., remains: {:6}", &path_i, &remains);
                *remain_now -= 1;
            }

            let SliceIRet { eigs_i, fweights_i, coeffs_i, projs_i, poscar, efermi } =
                Self::slice_i(
                    &path_i, ikpoints, brange.clone(), nspin, lncl, nplws_max, nions, nspd
                ).with_context(|| format!("Failed to slicing WAVECAR or PROCAR from {:?}.", &path_i))
                .unwrap();

            ret_eigs.lock().unwrap().slice_mut(nd::s![isw, .., .., ..]).assign(&eigs_i);
            ret_fweights.lock().unwrap().slice_mut(nd::s![isw, .., .., ..]).assign(&fweights_i);
            ret_coeffs.lock().unwrap().slice_mut(nd::s![isw, .., .., .., ..]).assign(&coeffs_i);
            ret_projs.lock().unwrap().slice_mut(nd::s![isw, .., .., .., .., ..]).assign(&projs_i);
            ret_efermis.lock().unwrap()[isw] = efermi;
            ret_xdatcar.lock().unwrap()[isw] = poscar;
        });


        Ok( SliceTotRet {
            eigs: Arc::try_unwrap(ret_eigs).unwrap().into_inner()?,
            fweights: Arc::try_unwrap(ret_fweights).unwrap().into_inner()?,
            coeffs: Arc::try_unwrap(ret_coeffs).unwrap().into_inner()?,
            projs: Arc::try_unwrap(ret_projs).unwrap().into_inner()?,
            xdatcar: Xdatcar::from(Arc::try_unwrap(ret_xdatcar).unwrap().into_inner()?),
            efermis: Arc::try_unwrap(ret_efermis).unwrap().into_inner()?,
        })
    }


    fn slice_i(path_i: &Path, ikpoints: &[usize], brange: Range<usize>, nspin: usize, lncl: bool,
        nplws_max: usize, nions: usize, nspd: usize) -> Result<SliceIRet> {
        let wav = Wavecar::from_file(&path_i.join("WAVECAR"))?;
        let proj = Procar::from_file(&path_i.join("PROCAR"))?;
        let poscar = Poscar::from_file(&path_i.join("CONTCAR"))?;

        let nbrange = brange.clone().count();
        let nkpoints = ikpoints.len();
        let nspinors = if lncl { 4 } else { nspin };

        anyhow::ensure!(nkpoints <= wav.nkpoints as usize);
        anyhow::ensure!(nbrange <= wav.nbands as usize);

        let mut eigs_i = nd::Array3::<f64>::zeros((nspin, nkpoints, nbrange));
        let mut fweights_i = eigs_i.clone();
        let mut coeffs_i = nd::Array4::<c64>::zeros((nspin, nkpoints, nbrange, nplws_max));
        let mut projs_i = nd::Array5::<f64>::zeros((nkpoints, nspinors, nbrange, nions, nspd));

        for (ii, &ik) in ikpoints.iter().enumerate() {
            eigs_i.slice_mut(nd::s![.., ik, ..])
                .assign(&wav.band_eigs.slice(nd::s![.., ik, brange.clone()]));
            fweights_i.slice_mut(nd::s![.., ik, ..])
                .assign(&wav.band_fweights.slice(nd::s![.., ik, brange.clone()]));
            projs_i.slice_mut(nd::s![ik, .., .., .., ..])
                .assign(&proj.pdos.projected.slice(nd::s![.., ik, brange.clone(), .., ..]));

            let nplw = wav.nplws[ik] as usize;
            let nspinor = if lncl { 2 } else { 1usize };
            for ispin in 0 .. nspin {
                for (jj, iband) in brange.clone().into_iter().enumerate() {
                    let coeff = wav._wav_kspace(ispin as u64, ik as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .context("Wavefunction reshape failed.")?;
                    coeffs_i.slice_mut(nd::s![ispin, ii, jj, 0..nplw])
                        .assign(&coeff);
                }
            }
        }

        let efermi = wav.efermi;

        Ok( SliceIRet {
            eigs_i,
            fweights_i,
            coeffs_i,
            projs_i,
            poscar,
            efermi,
        })
    }
}
