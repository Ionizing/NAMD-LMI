use std::fs;
use std::path::{Path, PathBuf};
use std::fmt;

use serde::{de::Error, Deserialize, Deserializer};
use toml;
use shared::{
    log,
    anyhow::ensure,
    anyhow::Context,
    anyhow::bail,
    Result,
};

use crate::core::NamdConfig;


#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum PropagateMethod {
    FiniteDifference,
    Exact,
    Expm,
    LiouvilleTrotter,
}


impl fmt::Display for PropagateMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use PropagateMethod::*;
        write!(f, "{}", match self {
            FiniteDifference => "FiniteDifference",
            Exact            => "Exact",
            Expm             => "Expm",
            LiouvilleTrotter => "LiouvilleTrotter",
        })
    }
}


impl PropagateMethod {
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
                "finitedifference" | "fd" => Ok(PropagateMethod::FiniteDifference),
                "exact"                   => Ok(PropagateMethod::Exact),
                "expm"                    => Ok(PropagateMethod::Expm),
                "liouvilletrotter" | "lt" => Ok(PropagateMethod::LiouvilleTrotter),
                _ => bail!("Invalid propmethod from input: {}, available methods: \
                                 FiniteDifference(or FD), Exact, Expm, LiouvilleTrotter(or LT)", &s)
        }
    }
}


#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct HamilConfig {
    #[serde(default = "HamilConfig::default_ikpoint")]
    ikpoint: usize,

    #[serde(default = "HamilConfig::default_basis")]
    basis_up: [usize; 2],

    #[serde(default = "HamilConfig::default_basis")]
    basis_dn: [usize; 2],

    #[serde(default = "HamilConfig::default_nelm")]
    nelm: usize,

    nac_fname: PathBuf,

    efield_fname: Option<PathBuf>,

    #[serde(default = "HamilConfig::default_hamil_fname")]
    hamil_fname: PathBuf,

    #[serde(deserialize_with="HamilConfig::parse_propmethod")]
    propmethod: PropagateMethod,

    scissor: Option<f64>,
}


impl HamilConfig {
    fn default_ikpoint() -> usize { 1 }
    fn default_basis() -> [usize; 2] { [0, 0] }
    fn default_nelm() -> usize { 10 }
    fn default_hamil_fname() -> PathBuf { PathBuf::from("HAMIL.h5") }

    fn parse_propmethod<'de, D>(deserializer: D) -> std::result::Result<PropagateMethod, D::Error>
    where D: Deserializer<'de> {
        let s = String::deserialize(deserializer)?;
        match s.to_lowercase().as_str() {
                "finitedifference" | "fd" => Ok(PropagateMethod::FiniteDifference),
                "exact"                   => Ok(PropagateMethod::Exact),
                "expm"                    => Ok(PropagateMethod::Expm),
                "liouvilletrotter" | "lt" => Ok(PropagateMethod::LiouvilleTrotter),
                _ => Err(D::Error::custom(
                        format!("Invalid propmethod from input: {}, available methods: \
                                 FiniteDifference(or FD), Exact, Expm, LiouvilleTrotter(or LT)", &s)
                        )),
        }
    }
}


impl HamilConfig {
    pub fn get_ikpoint(&self) -> usize { self.ikpoint }
    pub fn get_basis_up(&self) -> [usize; 2] { self.basis_up }
    pub fn get_basis_dn(&self) -> [usize; 2] { self.basis_dn }
    pub fn get_nelm(&self) -> usize { self.nelm }
    pub fn get_nac_fname(&self) -> &PathBuf { &self.nac_fname }
    pub fn get_efield_fname(&self) -> Option<&PathBuf> { self.efield_fname.as_ref() }
    pub fn get_hamil_fname(&self) -> &PathBuf { &self.hamil_fname }
    pub fn get_propmethod(&self) -> PropagateMethod { self.propmethod }
    pub fn get_scissor(&self) -> Option<f64> { self.scissor }

    pub fn  check_config(&self) -> Result<()> {
        let mut ret = Ok(());

        if self.ikpoint == 0 {
            ret = ret.context("Field 'ikpoint' counts from 1, thus cannot be 0.");
        }

        if self.basis_up[0] > self.basis_up[1] || self.basis_dn[0] > self.basis_dn[1] {
            ret = ret.context("Field 'basis_up' and 'basis_dn' must include valid and closed interval of bands.")
        }
        
        let nbasis = 
            if self.basis_up.contains(&0) {
                0
            } else {
                self.basis_up[1] - self.basis_up[0] + 1
            }
                +
            if self.basis_dn.contains(&0) {
                0
            } else {
                self.basis_dn[1] - self.basis_dn[0] + 1
            }
        ;

        if nbasis <= 1 {
            ret = ret.context("Field 'basis_up' and 'basis_dn' must include at least two bands to form a valid Hamiltonian.");
        }

        if self.nelm < 1 {
            ret = ret.context("Field 'nelm' cannot be 0.");
        }

        if !self.nac_fname.is_file() {
            ret = ret.context("Field 'nac_fname' does not point to a valid file.");
        }

        if let Some(efield) = self.efield_fname.as_ref() {
            if !efield.is_file() {
                ret = ret.context("Field 'efield_fname' does not point to a valid file.");
            }
        }

        if self.hamil_fname.as_os_str().is_empty() {
            ret = ret.context("Field 'hamil_fname' cannot be empty.");
        }

        if self.hamil_fname.is_file() {
            log::warn!("Field 'hamil_fname' points to an existing file and it will be overwritten.");
        }

        if let Some(s) = self.scissor {
            if s < 0.0 {
                ret = ret.context("Field 'scissor' cannot be a negative value.");
            }
        }

        ret
    }
}


impl fmt::Display for HamilConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "# NAMD-lumi config for Hamiltonian generation")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "ikpoint", self.ikpoint)?;
        writeln!(f, " {:>20} = {:?}", "basis_up", self.basis_up)?;
        writeln!(f, " {:>20} = {:?}", "basis_dn", self.basis_dn)?;
        writeln!(f, " {:>20} = {:?}", "nelm", self.nelm)?;
        writeln!(f, " {:>20} = {:?}", "nac_fname", self.nac_fname)?;
        if let Some(efield) = self.efield_fname.as_ref() {
            writeln!(f, " {:>20} = {:?}", "efield", efield)?;
        }
        writeln!(f, " {:>20} = {:?}", "hamil_fname", self.hamil_fname)?;
        writeln!(f, " {:>20} = {:?}", "propmethod", self.propmethod)?;

        if let Some(s) = self.scissor.as_ref() {
            writeln!(f, " {:>20} = {:?}", "scissor", s)?;
        }

        Ok(())
    }
}


impl Default for HamilConfig {
    fn default() -> Self {
        Self {
            ikpoint: 1,
            basis_up: [0, 0],
            basis_dn: [0, 0],
            nelm: 10,
            nac_fname: PathBuf::from("NAC.h5"),
            efield_fname: None,
            hamil_fname: PathBuf::from("HAMIL.h5"),
            propmethod: PropagateMethod::Expm,
            scissor: None,
        }
    }
}


impl NamdConfig for HamilConfig {
    fn from_file<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        ensure!(fname.as_ref().is_file(), "Config file {:?} for HamilConfig not available.", fname.as_ref());
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
        ikpoint = 2
        basis_up = [114, 514]
        basis_dn = [256, 512]
        nelm = 20
        nac_fname = "NAC_test.h5"
        hamil_fname = "HAMIL_test.h5"
        propmethod = "fd"
        scissor = 1.5
        "#;

        let actual_cfg: HamilConfig = toml::from_str(txt).unwrap();
        let expect_cfg: HamilConfig = HamilConfig {
            ikpoint: 2,
            basis_up: [114, 514],
            basis_dn: [256, 512],
            nelm: 20,
            nac_fname: PathBuf::from("NAC_test.h5"),
            efield_fname: None,
            hamil_fname: PathBuf::from("HAMIL_test.h5"),
            propmethod: PropagateMethod::FiniteDifference,
            scissor: Some(1.5),
        };

        assert_eq!(expect_cfg, actual_cfg);
    }
}
