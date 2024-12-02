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
};
use crate::waveslice::WavesliceConfig;


struct SliceIRet {
    eigs_i:     nd::Array3<f64>, // [nspin, nkpoints, nbrange]
    fweights_i: nd::Array3<f64>, // [nspin, nkpoints, nbrange]
    coeffs_i:  nd::Array4<c64>, // [nspin, nkpoints, nbrange, nplwmax]
    projs_i:   nd::Array5<f64>, // [nsw, nkpoints, nions, nspinor, nbrange]
    efermi:    f64,
}


struct SliceTotRet {
    eigs:     nd::Array4<f64>,  // [nsw, nspin, nkpoints, nbrange]
    fweights: nd::Array4<f64>,  // [nsw, nspin, nkpoints, nbrange]
    coeffs:   nd::Array5<c64>,  // [nsw, nspin, nkpoints, nbrange, nplwmax]
    projs:    nd::Array6<f64>,  // [nsw, nkpoints, nspinor, nbrange, nions, nspd]
    efermis:  nd::Array1<f64>,  // [nsw,]
}


pub struct Waveslice {
    /// Should count from 1
    ikpoints: Vec<usize>,

    /// number of spin channels
    nspin:  usize,

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

    /// Eigenvalue of each band for each Waveslice, [nsw, nkpoints, nbrange]
    eigs:     nd::Array4<f64>,

    /// Fermi occupation of each band for each Waveslice, [nsw, nkpoints, nbrange]
    fweights: nd::Array4<f64>,

    /// Coefficients for each k-point, [nsw, nspin, nkpoints, nbrange, nplwmax]
    coeffs:   nd::Array5<c64>,

    /// Slice of PROCAR for each WAVECAR, [nsw, nkpoints, nspinor, nbrange, nions, nspd]
    projs:    nd::Array6<f64>,
}


impl Waveslice {
    pub fn from_config(cfg: &WavesliceConfig) -> Result<Self> {
        if cfg.get_waveslicefname().is_file() {
            log::info!("Found pre-calculated NAC available in {:?}, reading WaveSlice from it ...",
                cfg.get_waveslicefname());
            let waveslice = Self::from_h5(cfg.get_waveslicefname())?;

            // TODO: add checkings

            return Ok(waveslice);
        }

        Self::calculate_from_scratch(&cfg)
    }


    pub fn from_h5<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        let f = H5File::open(fname)?;

        let ikpoints: Vec<usize> = f.dataset("ikpoints")?.read_raw()?;

        let nspin  = f.dataset("nspin")?.read_scalar::<usize>()?;
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

        let gvecs = {
            let gvecs_group = f.group("gvecs_group")?;

            let mut gvecs: Vec<nd::Array2<i64>> = Vec::new();
            for &ik in &ikpoints {
                let ikpath = format!("k{}", ik+1);
                let gvec: nd::Array2<i64> = gvecs_group.dataset(&ikpath)?.read()?;
                gvecs.push(gvec);
            }
            gvecs
        };

        let eigs: nd::Array4<f64> = f.dataset("eigs")?.read()?;
        let fweights: nd::Array4<f64> = f.dataset("fweights")?.read()?;
        let coeffs: nd::Array5<c64> = {
            let coeffs_r: nd::Array5<f64> = f.dataset("coeffs_r")?.read()?;
            let coeffs_i: nd::Array5<f64> = f.dataset("coeffs_i")?.read()?;
            coeffs_r.mapv(|x| c64::new(x, 0.0)) + coeffs_i.mapv(|x| c64::new(0.0, x))
        };

        let projs: nd::Array6<f64> = f.dataset("projs")?.read()?;

        Ok(Self {
            ikpoints,
            nspin,
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
        })
    }


    pub fn save_to_h5<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        let f = H5File::create(fname)?;

        f.new_dataset_builder().with_data(&self.ikpoints).create("ikpoints")?;

        f.new_dataset::<usize>().create("nspin")?.write_scalar(&self.nspin)?;
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

        let gvecs_group = f.create_group("gvecs_group")?;
        for (i, &ik) in self.ikpoints.iter().enumerate() {
            let ikpath = format!("k{}", ik+1);
            gvecs_group.new_dataset_builder().with_data(&self.gvecs[i]).create(ikpath.as_ref())?;
        }

        f.new_dataset_builder().with_data(&self.eigs).create("eigs")?;
        f.new_dataset_builder().with_data(&self.fweights).create("fweights")?;

        f.new_dataset_builder().with_data(&self.coeffs.mapv(|x| x.re)).create("coeffs_r")?;
        f.new_dataset_builder().with_data(&self.coeffs.mapv(|x| x.im)).create("coeffs_i")?;

        f.new_dataset_builder().with_data(&self.projs).create("projs")?;

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

        let gvecs = ikpoints.iter().map(|ik| {
            w1.generate_fft_grid((*ik) as u64)
                .into_iter().flatten()
                .collect::<nd::Array1::<i64>>()
                .into_shape((3, 3))
                .unwrap()
        }).collect::<Vec<nd::Array2<i64>>>();

        let lncl = w1.wavecar_type == WavecarType::NonCollinear;

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
            eigs, fweights, coeffs, projs, efermis
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

        
        let remain_count = Arc::new(Mutex::new(nsw - 1));
        (0 .. nsw).into_par_iter().for_each(|isw| {
            let path_i = rundir.join(format!("{:0ndigit$}", isw+1));
            {
                let mut remain_now = remain_count.lock().unwrap();
                let remains: usize = *remain_now;
                log::info!(" Slicing WAVECAR and PROCAR from {:?} ..., remains: {:6}", &path_i, &remains);
                *remain_now -= 1;
            }

            let SliceIRet { eigs_i, fweights_i, coeffs_i, projs_i, efermi } =
                Self::slice_i(
                    &path_i, ikpoints, brange.clone(), nspin, lncl, nplws_max, nions, nspd
                ).with_context(|| format!("Failed to slicing WAVECAR or PROCAR from {:?}.", &path_i))
                .unwrap();

            ret_eigs.lock().unwrap().slice_mut(nd::s![isw, .., .., ..]).assign(&eigs_i);
            ret_fweights.lock().unwrap().slice_mut(nd::s![isw, .., .., ..]).assign(&fweights_i);
            ret_coeffs.lock().unwrap().slice_mut(nd::s![isw, .., .., .., ..]).assign(&coeffs_i);
            ret_projs.lock().unwrap().slice_mut(nd::s![isw, .., .., .., .., ..]).assign(&projs_i);
            ret_efermis.lock().unwrap()[isw] = efermi;
        });
        //let SliceIRet { eigs_i, fweights_i, coeffs_i, projs_i } = Self::slice_i()


        Ok( SliceTotRet {
            eigs: Arc::try_unwrap(ret_eigs).unwrap().into_inner()?,
            fweights: Arc::try_unwrap(ret_fweights).unwrap().into_inner()?,
            coeffs: Arc::try_unwrap(ret_coeffs).unwrap().into_inner()?,
            projs: Arc::try_unwrap(ret_projs).unwrap().into_inner()?,
            efermis: Arc::try_unwrap(ret_efermis).unwrap().into_inner()?,
        })
    }


    fn slice_i(path_i: &Path, ikpoints: &[usize], brange: Range<usize>, nspin: usize, lncl: bool,
        nplws_max: usize, nions: usize, nspd: usize) -> Result<SliceIRet> {
        let wav = Wavecar::from_file(&path_i.join("WAVECAR"))?;
        let proj = Procar::from_file(&path_i.join("PROCAR"))?;

        let nbrange = brange.clone().count();
        let nkpoints = ikpoints.len();
        let nspinors = if lncl { 4 } else { nspin };

        anyhow::ensure!(nkpoints <= wav.nkpoints as usize);
        anyhow::ensure!(nbrange <= wav.nbands as usize);

        let mut eigs_i = nd::Array3::<f64>::zeros((nspin, nkpoints, nbrange));
        let mut fweights_i = eigs_i.clone();
        let mut coeffs_i = nd::Array4::<c64>::zeros((nspin, nkpoints, nbrange, nplws_max));
        let mut projs_i = nd::Array5::<f64>::zeros((nkpoints, nspinors, nbrange, nions, nspd));

        for &ik in ikpoints {
            eigs_i.slice_mut(nd::s![.., ik, ..])
                .assign(&wav.band_eigs.slice(nd::s![.., ik, brange.clone()]));
            fweights_i.slice_mut(nd::s![.., ik, ..])
                .assign(&wav.band_fweights.slice(nd::s![.., ik, brange.clone()]));
            projs_i.slice_mut(nd::s![ik, .., .., .., ..])
                .assign(&proj.pdos.projected.slice(nd::s![.., ik, brange.clone(), .., ..]));

            let nplw = wav.nplws[ik] as usize;
            let nspinor = if lncl { 2 } else { 1usize };
            for ispin in 0 .. nspin {
                for iband in brange.clone().into_iter() {
                    let coeff = wav._wav_kspace(ispin as u64, ik as u64, iband as u64, nplw / nspinor)
                        .into_shape((nplw,))
                        .context("Wavefunction reshape failed.")?;
                    coeffs_i.slice_mut(nd::s![ispin, ik, iband, 0..nplw])
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
            efermi,
        })
    }
}
