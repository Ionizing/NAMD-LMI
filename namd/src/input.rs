use shared::{
    anyhow,
    Result,
};

use toml;
use serde::{
    Deserialize,
    Serialize,
};
use crate::efield::Efield;

pub struct Input {
    pub rundir:       String,
    pub wavetype:     String,
    pub ikpoint:      usize,
    pub brange:       [usize; 2],
    pub basis_up:     [usize; 2],
    pub basis_dn:     [usize; 2],
    pub nsw:          usize,
    pub ndigit:       usize,
    pub namdtime:     usize,
    pub dt:           f64,
    pub nsample:      usize,
    pub ntraj:        usize,
    pub propmethod:   String,
    pub shmethod:     String,
    pub nelm:         usize,
    pub lreal:        bool,
    pub lprint_input: bool,
    pub lexcitation:  bool,
    pub lreorder:     bool,
    pub fname:        String,
    pub temperature:  f64,
    pub scissor:      f64,

    pub efield:       Option<Efield>,
    pub lcycle:       bool,
}
