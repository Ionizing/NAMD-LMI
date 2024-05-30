use std::fs;
use std::fmt;
use std::path::{Path, PathBuf};

use serde::Deserialize;
use toml;
use shared::{
    log,
    Result,
    anyhow::ensure,
};

use crate::core::NamdConfig;


#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
pub enum SHMethod {
    #[serde(alias="fssh", alias="FSSH")]
    FSSH,
    #[serde(alias="dish", alias="DISH")]
    DISH,
    #[serde(alias="dcsh", alias="DCSH")]
    DCSH,
}

const INISTEP_PY: &str = include_str!("./inisteps.py");


impl fmt::Display for SHMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use SHMethod::*;
        let s = match self {
            FSSH => "FSSH",
            DISH => "DISH",
            DCSH => "DCSH",
        };

        if f.alternate() {
            writeln!(f, "\"{}\"", s)
        } else {
            writeln!(f, "{}", s)
        }
    }
}


#[derive(Clone, Debug, Deserialize, PartialEq, Eq)]
pub struct SurfhopConfig {
    hamil_fname: PathBuf,
    namdtime: usize,
    nelm: usize,
    ntraj: usize,
    shmethod: SHMethod,
    outdir: PathBuf,
    lexcitation: bool,

    iniband: usize,
    inispin: usize,
    inisteps: Vec<usize>,
}


impl SurfhopConfig {
    pub fn get_hamil_fname(&self) -> &PathBuf { &self.hamil_fname }
    pub fn get_namdtime(&self) -> usize { self.namdtime }
    pub fn get_nelm(&self) -> usize { self.nelm }
    pub fn get_ntraj(&self) -> usize { self.ntraj }
    pub fn get_shmethod(&self) -> SHMethod { self.shmethod }
    pub fn get_outdir(&self) -> &PathBuf { &self.outdir }
    pub fn get_outdir_mut(&mut self) -> &mut PathBuf { &mut self.outdir }
    pub fn get_lexcitation(&self) -> bool { self.lexcitation }

    pub fn get_iniband(&self) -> usize { self.iniband }
    pub fn get_inispin(&self) -> usize { self.inispin }
    pub fn get_inisteps(&self) -> &[usize] { &self.inisteps }

    pub fn print_to_log(&self) {
        let config_print = format!("{}", self);
        let hashtag_line = "#".repeat(120);
        log::info!("SurfhopConfig file loaded. The formatted config is:\n\n{hashtag_line}\n{}{hashtag_line}\n\n", config_print);
    }

    pub fn write_inistep_py<P>(fname: P) -> Result<()>
    where P: AsRef<Path> {
        fs::write(fname, INISTEP_PY)?;
        Ok(())
    }
}


impl Default for SurfhopConfig {
    fn default() -> Self {
        Self {
            hamil_fname: "HAMIL.h5".into(),
            namdtime: 1000,
            nelm: 10,
            ntraj: 10000,
            shmethod: SHMethod::FSSH,
            outdir: ".".into(),
            lexcitation: true,

            iniband: 0,
            inispin: 1,
            inisteps: vec![1, 2, 3],
        }
    }
}


impl fmt::Display for SurfhopConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "####            NAMD-lumi config for surface-hopping calculation            ####")?;
        writeln!(f, "#### YOUT NEED TO CHANGE THE PARAMETERS IN THE FOLLOWING TO FIT YOUR SYSTEM ####")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "hamil_fname", self.hamil_fname)?;
        writeln!(f, " {:>20} = {:?}", "namdtime", self.namdtime)?;
        writeln!(f, " {:>20} = {:?}", "nelm", self.nelm)?;
        writeln!(f, " {:>20} = {:?}", "ntraj", self.ntraj)?;
        writeln!(f, " {:>20} = {:#}", "shmethod", self.shmethod)?;
        writeln!(f, " {:>20} = {:?}", "outdir", self.outdir)?;
        writeln!(f, " {:>20} = {:?}", "lexcitation", self.lexcitation)?;

        writeln!(f, " {:>20} = {:?}", "iniband", self.iniband)?;
        writeln!(f, " {:>20} = {:?}", "inispin", self.inispin)?;
        writeln!(f, " {:>20} = [", "inisteps")?;
        for step in self.inisteps.iter() {
            writeln!(f, "{:>10} ,", step)?;
        }
        writeln!(f, "]")?;

        Ok(())
    }
}


impl NamdConfig for SurfhopConfig {
    fn from_file<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        ensure!(fname.as_ref().is_file(), "Config file {:?} for SurfhopConfig not available", fname.as_ref());
        let raw = fs::read_to_string(fname)?;
        let cfg = toml::from_str::<Self>(&raw)?;

        Ok(cfg)
    }

    fn to_file<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        if fname.as_ref().is_file() {
            log::warn!("File {:?} exists, overwriting ...", fname.as_ref());
        }
        fs::write(fname, self.to_string())?;
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deserialize() {
        let txt = r#"
        hamil_fname = "HAMIL1.h5"
        namdtime = 2000
        nelm = 20
        ntraj = 20000
        shmethod = "DISH"
        outdir = "shout"
        lexcitation = true

        iniband = 3
        inispin = 2
        inisteps = [
            114,
            514,
        ]
        "#;

        let actual_cfg: SurfhopConfig = toml::from_str(txt).unwrap();
        let expect_cfg: SurfhopConfig = SurfhopConfig {
            hamil_fname: "HAMIL1.h5".into(),
            namdtime: 2000,
            nelm: 20,
            ntraj: 20000,
            shmethod: SHMethod::DISH,
            outdir: "shout".into(),
            lexcitation: true,

            iniband: 3,
            inispin: 2,
            inisteps: vec![114, 514],
        };

        assert_eq!(expect_cfg, actual_cfg);

        println!("{}", &expect_cfg);
    }
}
