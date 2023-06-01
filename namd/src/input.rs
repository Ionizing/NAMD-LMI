use std::fs::read_to_string;

use shared::{
    anyhow,
    Result,
};

use toml;
use serde::{
    de::Error,
    ser,
    Deserialize,
    Deserializer,
    Serialize,
};
use crate::efield::Efield;

#[derive(Deserialize)]
pub struct Input {
    pub rundir:       String,
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
    pub nacfname:        String,
    pub temperature:  f64,

    #[serde(default)]
    pub scissor:      f64,

    #[serde(default)]
    #[serde(deserialize_with = "Input::efield_from_file")]
    pub efield:       Option<(String, Efield)>,
    pub lcycle:       bool,
}


impl Input {
    fn efield_from_file<'de, D>(deserializer: D) -> Result<Option<(String, Efield)>, D::Error>
    where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        let txt = read_to_string(&s).map_err(D::Error::custom)?;
        Ok(Some((
                s.to_string(),
                Efield::from_str(&txt).map_err(D::Error::custom)?
                )))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_efield() {
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
            "#;

        let input: Input = toml::from_str(raw).unwrap();
        assert_eq!(input.efield.unwrap().0, "tests/1.5eV.txt");
    }
}
