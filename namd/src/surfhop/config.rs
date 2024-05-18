use std::fs;
use std::fmt;
use std::path::{Path, PathBuf};

use serde::{de::Error, Deserialize, Deserializer};
use toml;
use shared::{
    log,
    Result,
};

use crate::core::NamdConfig;


#[derive(Clone, Copy)]
pub enum SHMethod {
    FSSH,
    DISH,
    DCSH,
}


#[derive(Clone, Deserialize)]
pub struct Config {
    hamil_fname: PathBuf,
    nelm: usize,
    shmethod: SHMethod,
    outdir: PathBuf,

    iniband: usize,
    inispin: usize,
    inisteps: Vec<usize>,
}


impl Default for Config {
    fn default() -> Self {
        Self {
            hamil_fname: "HAMIL.h5".into(),
            nelm: 10,
            shmethod: SHMethod::FSSH,
            outdir: ".".into(),
            iniband: 0,
            inispin: 1,
            inisteps: vec![1, 2, 3],
        }
    }
}


impl fmt::Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
    }
}


impl NamdConfig for Config {
    fn from_file<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        todo!()
    }

    fn to_file<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        Ok(())
    }
}
