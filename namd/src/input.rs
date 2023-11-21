use std::path::{
    Path,
    PathBuf,
};
use std::fs::read_to_string;

use shared::{
    bail,
    Result,
};

use toml;
use serde::{
    de::Error,
    //ser,
    Deserialize,
    Deserializer,
    //Serialize,
};
use crate::efield::Efield;
use crate::hamiltonian::PropagateMethod;
use crate::surface_hopping::SHMethod;

// All indices starst from 1
#[derive(Deserialize)]
pub struct Input {
    pub rundir:       PathBuf,
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

    #[serde(deserialize_with = "Input::parse_propmethod")]
    pub propmethod:   PropagateMethod,

    #[serde(deserialize_with = "Input::parse_shmethod")]
    pub shmethod:     SHMethod,
    pub nelm:         usize,
    pub lreal:        bool,
    pub lprint_input: bool,
    pub lexcitation:  bool,
    pub lreorder:     bool,
    pub nacfname:     PathBuf,
    pub temperature:  f64,

    #[serde(default)]
    pub scissor:      f64,

    pub inibands:     Vec<usize>,
    pub inispins:     Vec<usize>,
    pub inisteps:     Vec<usize>,

    #[serde(default)]
    #[serde(deserialize_with = "Input::efield_from_file")]
    pub efield:       Option<(PathBuf, Efield)>,
    pub lcycle:       bool,
}


impl Input {
    pub fn from_file<P>(fname: &P) -> Result<Self>
    where
        P: AsRef<Path> + ?Sized
    {
        if !fname.as_ref().is_file() {
            bail!("Input file {:?} not available.", fname.as_ref());
        }
        let raw         = read_to_string(fname)?;
        let input: Self = toml::from_str(&raw)?;
        Ok(input)
    }

    fn efield_from_file<'de, D>(deserializer: D) -> Result<Option<(PathBuf, Efield)>, D::Error>
    where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        let txt = read_to_string(&s).map_err(D::Error::custom)?;
        Ok(Some((
                PathBuf::from(s),
                Efield::from_str(&txt).map_err(D::Error::custom)?
                )))
    }

    fn parse_propmethod<'de, D>(deserializer: D) -> Result<PropagateMethod, D::Error>
    where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        match s.to_lowercase().as_str() {
            "finitedifference" | "fd" => Ok(PropagateMethod::FiniteDifference),
            "exact"                   => Ok(PropagateMethod::Exact),
            "expm"                    => Ok(PropagateMethod::Expm),
            "liouvilletrotter" | "lt" => Ok(PropagateMethod::LiouvilleTrotter),
            _ => Err(D::Error::custom(
                    format!("Invalid propmethod from input: {}, available methods: FiniteDifference(or FD), Exact, Expm, LiouvilleTrotter(or LT)", &s)
                    )),
        }
    }

    fn parse_shmethod<'de, D>(deserializer: D) -> Result<SHMethod, D::Error>
    where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        match s.to_ascii_lowercase().as_str() {
            "fssh" => Ok(SHMethod::FSSH),
            "dish" => Ok(SHMethod::DISH),
            "dcsh" => Ok(SHMethod::DCSH),
            "gfsh" => Ok(SHMethod::GFSH),
            _ => Err(D::Error::custom(
                    format!("Invalid shmethod from input: {}, available methods: FSSH, DISH, DCSH, GFSH", &s)
                    )),
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_input() {
        let raw = r#"
        rundir       = "~/"
        ikpoint      = 1
        brange       = [50, 100]
        basis_up     = [60, 90]
        basis_dn     = [0, 0]
        nsw          = 3000
        ndigit       = 4
        namdtime     = 1000
        dt           = 0.1
        nsample      = 100
        ntraj        = 10000
        propmethod   = "exact"
        shmethod     = "fssh"
        nelm         = 10
        lreal        = false
        lprint_input = true
        lexcitation  = true
        lreorder     = false
        nacfname     = "NAC.h5"
        temperature  = 300
        lcycle       = false
        efield       = "tests/1.5eV.txt"

        inibands     = [114, 114]
        inispins     = [1, 1]
        inisteps     = [514, 810]
            "#;

        let input: Input = toml::from_str(raw).unwrap();
        assert_eq!(input.efield.unwrap().0, PathBuf::from("tests/1.5eV.txt"));
    }


    #[test]
    fn test_parse_input2() {
        let raw = r#"
        rundir       = "~/"
        ikpoint      = 1
        brange       = [50, 100]
        basis_up     = [60, 90]
        basis_dn     = [0, 0]
        nsw          = 3000
        ndigit       = 4
        namdtime     = 1000
        dt           = 0.1
        nsample      = 100
        ntraj        = 10000
        propmethod   = "exact"
        shmethod     = "fssh"
        nelm         = 10
        lreal        = false
        lprint_input = true
        lexcitation  = true
        lreorder     = false
        nacfname     = "NAC.h5"
        temperature  = 300
        lcycle       = false
        efield       = "tests/2.0eV.txt"

        inibands     = [114, 114]
        inispins     = [1, 1]
        inisteps     = [514, 810]
            "#;

        let input: Input = toml::from_str(raw).unwrap();
        assert_eq!(input.efield.unwrap().0, PathBuf::from("tests/2.0eV.txt"));
    }
}
