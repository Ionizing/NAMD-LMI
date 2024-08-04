use std::fs;
use std::path::{Path, PathBuf};
use std::fmt;

use serde::{de::Error, Deserialize, Deserializer};
use toml;
use regex::Regex;
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
#[serde(deny_unknown_fields)]
pub struct HamilConfig {
    #[serde(default = "HamilConfig::default_ikpoint")]
    ikpoint: usize,

    // List of bands to form the basis.
    //
    // - Positive one indicates for the spin-up band;
    // - Negative one indicates for the spin-down band.
    // - Zeros are not allowed
    //
    // You can specify a consecutive bands with multiple `range`s or singletons:
    //
    // A `range` is a pattern of `start..end` where
    //   - start and end are integers with same sign (both + or both -)
    //   - |start| <= |end|
    //   - no other characters (whitespace or any other thinds) around `..` 
    //
    // And multiple ranges and singletons are separated with one or more whitespaces.
    //
    // EXAMPLE:
    //   - `1..4` expands to `[1, 2, 3, 4]`;
    //   - `1..1` expands to `[1]`;
    //   - `4..1` expands to empty list `[]`;
    //
    //   - `-1..-4` expands to `[-1, -2, -3, -4]`.
    //   - `-4..-4` expands to `[-4]`
    //   - `-4..-1` expands to empty list `[]`;
    //
    //   - Tokens like `-4..4` `1 ..3`, `0`, `-3..0`, `-3..3` are not allowed.
    #[serde(deserialize_with="HamilConfig::parse_basis_list")]
    basis_list: Vec<i32>,

    // Corresponding labels to the basis_list, must have same length
    // The labels must not contain newline related characters like CRLF.
    basis_labels: Option<Vec<String>>,

    nac_fname: PathBuf,

    efield_fname: Option<PathBuf>,

    #[serde(default = "HamilConfig::default_hamil_fname")]
    hamil_fname: PathBuf,

    #[serde(deserialize_with="HamilConfig::parse_propmethod")]
    propmethod: PropagateMethod,

    #[serde(default = "HamilConfig::default_reorder")]
    reorder: bool,

    scissor: Option<f64>,
}


impl HamilConfig {
    fn default_ikpoint() -> usize { 1 }
    fn default_hamil_fname() -> PathBuf { PathBuf::from("HAMIL.h5") }
    fn default_reorder() -> bool { false }

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

    fn parse_basis_list<'de, D>(deserializer: D) -> std::result::Result<Vec<i32>, D::Error>
    where D: Deserializer<'de> {
        let input = String::deserialize(deserializer)?;
        let mut ret = vec![];

        let re_range = Regex::new(r"^(-?\d+)\.\.(-?\d+)$").unwrap();
        let re_digit = Regex::new(r"^-?\d+$").unwrap();

        for s in input.split_ascii_whitespace() {
            if re_digit.is_match(s) {
                let num = s.parse::<i32>().unwrap();
                if num == 0 {
                    return Err(D::Error::custom(
                        format!("Invalid integer '{}' from `basis_list`, 0 is not allowed.", s)
                    ));
                }
                ret.push(s.parse().unwrap())
            } else if re_range.is_match(s) {
                let m = re_range.captures(s).unwrap();
                let start = m.get(1).unwrap().as_str().parse::<i32>().unwrap();
                let end   = m.get(2).unwrap().as_str().parse::<i32>().unwrap();

                if start * end <= 0 {
                    return Err(D::Error::custom(
                        format!("Invalid range '{}' from `basis_list`, start and end must have save sign.", s)
                    ));
                }

                let sign = start.signum();
                let to_be_extend = (start.abs() ..= end.abs()).map(|x| x*sign).collect::<Vec<_>>();
                ret.extend(to_be_extend);
            } else {
                return Err(D::Error::custom(
                    format!("Invalid token '{}' from `basis_list`, it should be either range (start..end) or integer.", s)
                ));
            }
        }

        // check duplication
        let mut ret_dedup = ret.clone();
        ret_dedup.sort();
        ret_dedup.dedup();

        if ret.len() != ret_dedup.len() {
            return Err(D::Error::custom(
                format!("Invalid `basis_list`, each band in the basis should be unique. \nThe expanded list is: {:?}", ret_dedup)
            ));
        }
        
        return Ok(ret)
    }


    pub fn get_ikpoint(&self) -> usize { self.ikpoint }
    pub fn get_basis_list(&self) -> &[i32] { &self.basis_list }

    pub fn get_basis_labels(&self) -> Option<&Vec<String>> {
        self.basis_labels.as_ref()
    }

    pub fn get_nac_fname(&self) -> &PathBuf { &self.nac_fname }
    pub fn get_efield_fname(&self) -> Option<&PathBuf> { self.efield_fname.as_ref() }
    pub fn get_hamil_fname(&self) -> &PathBuf { &self.hamil_fname }
    pub fn get_propmethod(&self) -> PropagateMethod { self.propmethod }
    pub fn get_reorder(&self) -> bool { self.reorder }
    pub fn get_scissor(&self) -> Option<f64> { self.scissor }

    pub fn  check_config(&self) -> Result<()> {
        let mut ret = Ok(());

        if self.ikpoint == 0 {
            ret = ret.context("Field 'ikpoint' counts from 1, thus cannot be 0.");
        }

        let nbasis = self.basis_list.len();
        if nbasis <= 1 {
            ret = ret.context("Field `basis_list` must include at least two bands to form a valid Hamiltonian.");
        }

        if let Some(labels) = self.basis_labels.as_ref() {
            if labels.len() != nbasis {
                ret = ret.context(
                    format!("Number of entries in `basis_labels` is not consistent with basis size: {} != {}",
                        labels.len(), nbasis)
                );
            }

            if labels.iter().any(|s| s.is_empty()) {
                ret = ret.context("`basis_labels` contains empty entries, which is not allowed.");
            }

            if labels.iter().any(|s| s.contains('\n') || s.contains('\r')) {
                ret = ret.context("`basis_labels` contains newline character `CRLF`, which is not allowed.");
            }
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


    pub fn print_to_log(&self) {
        let input_print = format!("{}", self);
        let hashtag_line = "#".repeat(120);
        log::info!("Input file loaded. The formatted input is:\n\n{hashtag_line}\n{}\n{hashtag_line}\n\n", input_print);
    }
}


impl fmt::Display for HamilConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "# NAMD-lumi config for Hamiltonian generation")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "ikpoint", self.ikpoint)?;
        writeln!(f, " {:>20} = {:?}", "basis_list", self.basis_list)?;
        if let Some(labels) = self.basis_labels.as_ref() {
            writeln!(f, " {:>20} = {:?}", "basis_labels", labels)?;
        }
        writeln!(f, " {:>20} = {:?}", "nac_fname", self.nac_fname)?;
        if let Some(efield) = self.efield_fname.as_ref() {
            writeln!(f, " {:>20} = {:?}", "efield_fname", efield)?;
        } else {
            writeln!(f, "#{:>20} = # to be filled", "efield_fname")?;
        }
        writeln!(f, " {:>20} = {:?}", "hamil_fname", self.hamil_fname)?;
        writeln!(f, " {:>20} = \"{:?}\"", "propmethod", self.propmethod)?;
        //writeln!(f, " {:>20} = {:?}", "reorder", self.reorder)?;

        if let Some(s) = self.scissor.as_ref() {
            writeln!(f, " {:>20} = {:?}", "scissor", s)?;
        } else {
            writeln!(f, "#{:>20} = 1.5 # unit: eV", "scissor")?;
        }

        Ok(())
    }
}


impl Default for HamilConfig {
    fn default() -> Self {
        Self {
            ikpoint: 1,
            basis_list: vec![0],
            basis_labels: Some(vec!["0".into()]),
            nac_fname: PathBuf::from("NAC.h5"),
            efield_fname: None,
            hamil_fname: PathBuf::from("HAMIL.h5"),
            propmethod: PropagateMethod::Expm,
            reorder: false,
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
        log::info!("Writing config to file {:?}", fname.as_ref());
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
        basis_list = "-1..-4 1..4"
        nac_fname = "NAC_test.h5"
        hamil_fname = "HAMIL_test.h5"
        propmethod = "fd"
        scissor = 1.5
        "#;

        let actual_cfg: HamilConfig = toml::from_str(txt).unwrap();
        let expect_cfg: HamilConfig = HamilConfig {
            ikpoint: 2,
            basis_list: ((-4..=-1).rev().chain(1..=4)).collect(),
            basis_labels: None,
            nac_fname: PathBuf::from("NAC_test.h5"),
            efield_fname: None,
            hamil_fname: PathBuf::from("HAMIL_test.h5"),
            propmethod: PropagateMethod::FiniteDifference,
            reorder: false,
            scissor: Some(1.5),
        };

        assert_eq!(expect_cfg, actual_cfg);
    }
}
