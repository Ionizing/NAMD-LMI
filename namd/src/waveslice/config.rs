use std::fmt;
use std::fs;
use std::path::{PathBuf, Path};
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
pub struct WavesliceConfig {
    rundir: PathBuf,

    #[serde(default = "WavesliceConfig::default_ikpoints")]
    /// When "0" is the only number in this field, all the k-points are selected.
    ikpoints: Vec<usize>,

    brange: [usize; 2],

    nsw: usize,

    #[serde(default = "WavesliceConfig::default_ndigit")]
    ndigit: usize,

    #[serde(default = "WavesliceConfig::default_phasecorrection")]
    phasecorrection: bool,

    #[serde(default = "WavesliceConfig::default_unitary_transform")]
    unitary_transform: bool,

    #[serde(default = "WavesliceConfig::default_rearrangement")]
    rearrangement: bool,

    #[serde(default = "WavesliceConfig::default_waveslicefname")]
    waveslicefname: PathBuf,
}


impl WavesliceConfig {
    fn default_ikpoints() -> Vec<usize> { vec![1] }
    fn default_ndigit() -> usize { 4 }
    fn default_phasecorrection() -> bool { true }
    fn default_waveslicefname() -> PathBuf { PathBuf::from("waveslice.h5") }
    fn default_unitary_transform() -> bool { false }
    fn default_rearrangement() -> bool { false }
}


impl WavesliceConfig {
    pub fn check_config(&self) -> Result<()> {
        let mut ret = Ok(());

        if !self.rundir.is_file() {
            ret = ret.context("Field 'rundir' is not a valid directory.");
        }

        if self.ikpoints.is_empty() {
            ret = ret.context("Field 'ikpoints' is empty.");
        }

        if self.ikpoints.contains(&0) && self.ikpoints.len() != 1 {
            ret = ret.context("Field 'ikpoints' contains 0 while it should cont from 1.");
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

        if self.waveslicefname.as_os_str().is_empty() {
            ret = ret.context("Field 'waveslicefname' cannot be empy.");
        }

        ret
    }


    pub fn get_brange(&self) -> &[usize; 2] { &self.brange }
    pub fn get_ikpoints(&self) -> &[usize] { &self.ikpoints }
    pub fn get_ndigit(&self) -> usize { self.ndigit }
    pub fn get_nsw(&self) -> usize { self.nsw }
    pub fn get_phasecorrection(&self) -> bool { self.phasecorrection }
    pub fn get_rearrangement(&self) -> bool { self.rearrangement }
    pub fn get_unitary_transform(&self) -> bool { self.unitary_transform }
    pub fn get_rundir(&self) -> &Path { &self.rundir }
    pub fn get_waveslicefname(&self) -> &Path { &self.waveslicefname }
}


impl Default for WavesliceConfig {
    fn default() -> Self {
        WavesliceConfig {
            rundir: PathBuf::from("../run"),
            ikpoints: vec![1],
            brange: [0, 0],
            nsw: 2000,
            ndigit: 4,
            phasecorrection: true,
            unitary_transform: false,
            rearrangement: false,
            waveslicefname: PathBuf::from("waveslice.h5"),
        }
    }
}


impl fmt::Display for WavesliceConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "####          NAMD-LMI config for WAVECARs & PROCARs slicing calculation        ####")?;
        writeln!(f, "####    YOU NEED TO CHANGE THE PARAMETERS IN THE FOLLOWING TO FIT YOU SYSTEM    ####")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "rundir",   self.rundir)?;
        writeln!(f, " {:>20} = {:?}", "ikpoints",  self.ikpoints)?;
        writeln!(f, " {:>20} = {:?}", "brange",   self.brange)?;
        writeln!(f, " {:>20} = {}",   "nsw",      self.nsw)?;
        writeln!(f, " {:>20} = {}",   "ndigit",   self.ndigit)?;
        writeln!(f, " {:>20} = {}",   "phasecorrection", self.phasecorrection)?;
        writeln!(f, " {:>20} = {}",   "unitary_transform", self.unitary_transform)?;
        writeln!(f, " {:>20} = {}",   "rearrangement", self.rearrangement)?;
        writeln!(f, " {:>20} = {:?}", "waveslicefname", self.waveslicefname)?;

        Ok(())
    }
}


impl NamdConfig for WavesliceConfig {
    fn from_file<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        ensure!(fname.as_ref().is_file(), "Config file {:?} for WavesliceConfig not available.", fname.as_ref());
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
        ikpoints = [0]
        brange = [100, 200]
        nsw = 3000
        ndigit = 5
        phasecorrection = false
        unitary_transform = true
        rearrangement = true
        waveslicefname = "Waveslice.h5"
        "#;

        let actual_cfg: WavesliceConfig = toml::from_str(txt).unwrap();
        let expect_cfg: WavesliceConfig = WavesliceConfig {
            rundir: PathBuf::from("../run"),
            ikpoints: vec![0],
            brange: [100, 200],
            nsw: 3000,
            ndigit: 5,
            phasecorrection: false,
            unitary_transform: true,
            rearrangement: true,
            waveslicefname: PathBuf::from("Waveslice.h5"),
        };

        assert_eq!(expect_cfg, actual_cfg);
    }
}
