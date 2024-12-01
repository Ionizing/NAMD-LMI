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
    eig_i:     nd::Array3<f64>, // [nspin, nkpoints, nbrange]
    fweight_i: nd::Array3<f64>, // [nspin, nkpoints, nbrange]
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
    pub fn from_h5<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        todo!()
    }

    pub fn save_to_h5<P>(&self, fname: P) -> Result<Self>
    where P: AsRef<Path> {
        todo!()
    }


    fn calculate_from_scratch(cfg: &WavesliceConfig) -> Result<Self> {
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
        };
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


        todo!();


        Ok( SliceTotRet {
            eigs: Arc::try_unwrap(ret_eigs).unwrap().into_inner()?,
            fweights: Arc::try_unwrap(ret_fweights).unwrap().into_inner()?,
            coeffs: Arc::try_unwrap(ret_coeffs).unwrap().into_inner()?,
            projs: Arc::try_unwrap(ret_projs).unwrap().into_inner()?,
            efermis: Arc::try_unwrap(ret_efermis).unwrap().into_inner()?,
        })
    }
}
