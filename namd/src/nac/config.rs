use std::fmt;
use std::path::{Path, PathBuf};
use std::fs;
use serde::Deserialize;
use toml;
use shared::{
    log,
    anyhow::ensure,
    anyhow::Context,
    Result,
};
use crate::core::config::NamdConfig;

#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NacConfig {
    rundir:   PathBuf,

    #[serde(default = "NacConfig::default_ikpoint")]
    ikpoint:  usize,

    brange:   [usize; 2],

    nsw:      usize,

    #[serde(default = "NacConfig::default_ndigit")]
    ndigit:   usize,

    potim:    f64,

    temperature: f64,

    #[serde(default = "NacConfig::default_phasecorrection")]
    phasecorrection: bool,

    #[serde(default = "NacConfig::default_nacfname")]
    nacfname: PathBuf,
}


impl NacConfig {
    fn default_ikpoint() -> usize { 1 }
    fn default_ndigit() -> usize { 4 }
    fn default_phasecorrection() -> bool { true }
    fn default_nacfname() -> PathBuf { PathBuf::from("NAC.h5") }
}


impl NacConfig {
    pub fn get_rundir(&self) -> &PathBuf { &self.rundir }
    pub fn get_ikpoint(&self) -> usize { self.ikpoint }
    pub fn get_brange(&self) -> [usize;2] { self.brange }
    pub fn get_nsw(&self) -> usize { self.nsw }
    pub fn get_ndigit(&self) -> usize { self.ndigit }
    pub fn get_potim(&self) -> f64 { self.potim }
    pub fn get_temperature(&self) -> f64 { self.temperature }
    pub fn get_phasecorrection(&self) -> bool { self.phasecorrection }
    pub fn get_nacfname(&self) -> &PathBuf { &self.nacfname }

    pub fn check_config(&self) -> Result<()> {
        let mut ret = Ok(());
        if !self.rundir.is_file() {
            ret = ret.context("Field 'rundir' is not a valid directory.");
        }

        if self.ikpoint == 0 {
            ret = ret.context("Field 'ikpoint' counts from 1 thus cannot be 0.");
        }

        if self.brange[0] >= self.brange[0] {
            ret = ret.context("Field 'brange' must include the range of two or more bands.");
        }

        if self.brange[0] == 0 || self.brange[1] == 0 {
            ret = ret.context("Field 'brange' counts from 1 thus cannot be 0.");
        }

        if self.nsw < 4 {
            ret = ret.context("Field 'nsw' too short, 2000 or more is recommended.");
        }

        if self.ndigit == 0 {
            ret = ret.context("Field 'ndigit' cannot be 0.");
        }

        if (self.nsw as f64).log10() as i32 > self.ndigit as i32 {
            ret = ret.context("Field 'ndigit' not match with field 'nsw', please check.");
        }

        if self.potim <= 0.0 {
            ret = ret.context("Field 'potim' cannot be less than or equal to 0.");
        }

        if self.temperature < 0.0 {
            ret = ret.context("Temperature cannot be negative.");
        }

        if self.nacfname.as_os_str().is_empty() {
            ret = ret.context("Field 'nacfname' cannot be empy.");
        }

        ret
    }

    pub fn print_to_log(&self) {
        let input_print = format!("{}", self);
        let hashtag_line = "#".repeat(120);
        log::info!("Input file loaded. The formatted input is:\n\n{hashtag_line}\n{}\n{hashtag_line}\n\n", input_print);
    }
}


impl Default for NacConfig {
    fn default() -> Self {
        NacConfig {
            rundir: PathBuf::from("../run"),
            ikpoint: 1,
            brange: [0, 0],
            nsw: 2000,
            ndigit: 4,
            potim: 1.0,
            temperature: 0.0,
            phasecorrection: true,
            nacfname: PathBuf::from("NAC.h5"),
        }
    }
}


impl fmt::Display for NacConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "####        NAMD-lumi config for Non-Adiabatic Coupling (NAC) calculation       ####")?;
        writeln!(f, "####    YOU NEED TO CHANGE THE PARAMETERS IN THE FOLLOWING TO FIT YOU SYSTEM    ####")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "rundir",   self.rundir)?;
        writeln!(f, " {:>20} = {:?}", "ikpoint",  self.ikpoint)?;
        writeln!(f, " {:>20} = {:?}", "brange",   self.brange)?;
        writeln!(f, " {:>20} = {}",   "nsw",      self.nsw)?;
        writeln!(f, " {:>20} = {}",   "ndigit",   self.ndigit)?;
        writeln!(f, " {:>20} = {}",   "potim",    self.potim)?;
        writeln!(f, " {:>20} = {}",   "temperature", self.temperature)?;
        writeln!(f, " {:>20} = {}",   "phasecorrection", self.phasecorrection)?;
        writeln!(f, " {:>20} = {:?}", "nacfname", self.nacfname)?;

        Ok(())
    }
}


impl NamdConfig for NacConfig {
    fn from_file<P>(fname: P) -> Result<Self> 
    where P: AsRef<Path> {
        ensure!(fname.as_ref().is_file(), "Config file {:?} for NacConfig not available.", fname.as_ref());
        let raw = fs::read_to_string(fname)?;
        let cfg = toml::from_str::<Self>(&raw)?;
        cfg.check_config()?;
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
        rundir = "../run"
        ikpoint = 2
        brange = [100, 200]
        nsw = 3000
        ndigit = 5
        potim = 1.5
        temperature = 150
        phasecorrection = false
        nacfname = "NAC2.h5"
        "#;

        let actual_cfg: NacConfig = toml::from_str(txt).unwrap();
        let expect_cfg: NacConfig = NacConfig {
            rundir: PathBuf::from("../run"),
            ikpoint: 2,
            brange: [100, 200],
            nsw: 3000,
            ndigit: 5,
            potim: 1.5,
            temperature: 150.0,
            phasecorrection: false,
            nacfname: PathBuf::from("NAC2.h5"),
        };

        assert_eq!(expect_cfg, actual_cfg);
    }
}
