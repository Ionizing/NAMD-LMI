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
};

use vasp_parsers::{
    procar::Procar,
    Wavecar,
    WavecarType,
};
use crate::waveslice::WavesliceConfig;


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
    ngrid:     Vec<usize>,

    /// Fermi level of each WAVECAR
    efermis: Vec<f64>,

    /// K-vector of selected K points, [3, 3]
    kvecs:    nd::Array2<f64>,

    /// Number of plane wave coefficients for selected K points, [nkpoints,]
    num_plws: Vec<usize>,

    /// G vectors of each plane wave for selected K points, [nkpoints, [nplw, 3]]
    gvecs:    Vec<nd::Array2<i64>>,

    /// Slice of PROCAR for each WAVECAR, [nkpoints, [nsw, nspinor, nions, nbrange, nspd]]
    projs:    Vec<nd::Array5<f64>>,

    /// Eigenvalue of each band for each Waveslice, [nkpoints, [nsw, nbrange]]
    eigs:     Vec<nd::Array3<f64>>,

    /// Fermi occupation of each band for each Waveslice, [nkpoints, [nsw, nbrange]]
    fweights: Vec<nd::Array3<f64>>,
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
        todo!()
    }
}
