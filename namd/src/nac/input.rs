use std::fmt;
use std::path::{Path, PathBuf};
use std::fs;
use serde::Deserialize;
use toml;
use shared::{
    log,
    bail,
    Result,
};
use crate::core::input::Input;

#[derive(Clone, Deserialize)]
pub struct NacInput {
    rundir:   PathBuf,
    brange:   [usize; 2],
    nsw:      usize,
    ndigit:   usize,
    potim:    f64,
    nacfname: PathBuf,
}


impl Default for NacInput {
    fn default() -> Self {
        NacInput {
            rundir: PathBuf::from("../run"),
            brange: [0, 0],
            nsw: 2000,
            ndigit: 4,
            potim: 1.0,
            nacfname: PathBuf::from("NAC.h5"),
        }
    }
}


impl fmt::Display for NacInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result<()> {
        writeln!(f, "# NAMD-lumi input for Non-Adiabatic Coupling (NAC) calculation")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "rundir",   self.rundir)?;
        writeln!(f, " {:>20} = {:?}", "brange",   self.brange)?;
        writeln!(f, " {:>20} = {}",   "nsw",      self.nsw)?;
        writeln!(f, " {:>20} = {}",   "ndigit",   self.ndigit)?;
        writeln!(f, " {:>20} = {}",   "potim",    self.potim)?;
        writeln!(f, " {:>20} = {:?}", "nacfname", self.nacfname)?;

        Ok(())
    }
}


impl<'a> Input<'a> for NacInput {
    fn from_file<P>(fname: P) -> Result<Self> 
    where P: AsRef<Path> {
        if !fname.as_ref().is_file() {
            bail!("Input file {:?} for NacInput not available.", fname.as_ref());
        }
        let raw   = fs::read_to_string(fname)?;
        let input = toml::from_str::<Self>(&raw)?;
        Ok(input)
    }

    fn to_file<P>(&self, fname: P) -> Result<()>
    where P: AsRef<Path> {
        if fname.as_ref().is_file() {
            log::warn!("File {:?} exists, overwriting ...", fname.as_ref());
        }
        fs::write(fname.as_ref(), self.to_string())?;
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
        brange = [100, 200]
        nsw = 3000
        ndigit = 5
        potim = 1.5
        nacfname = "NAC2.h5"
        "#;

        let nac_input: NacInput = toml::from_str(txt).unwrap();

        assert_eq!(nac_input.rundir, Path::new("../run"));
        assert_eq!(nac_input.brange, &[100, 200]);
        assert_eq!(nac_input.nsw, 3000);
        assert_eq!(nac_input.ndigit, 5);
        assert_eq!(nac_input.potim, 1.5);
        assert_eq!(nac_input.nacfname, Path::new("NAC2.h5"));
    }
}
