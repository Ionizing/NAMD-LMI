use std::fmt;
use std::path::{
    Path,
    PathBuf,
};
use std::fs::{
    create_dir_all,
    read_to_string,
};

use shared::{
    bail,
    Result,
    info,
    warn,
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
#[derive(Clone, Deserialize)]
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
    pub scissor:      Option<f64>,

    #[serde(default = "Input::default_outputdir")]
    pub outdir:       PathBuf,

    pub inibands:     Vec<usize>,
    pub inispins:     Vec<usize>,
    pub inisteps:     Vec<usize>,

    #[serde(default)]
    #[serde(deserialize_with = "Input::efield_from_file")]
    pub efield:       Option<(PathBuf, String, Efield)>,
}


impl Input {
    pub fn from_file<P>(fname: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        if !fname.as_ref().is_file() {
            bail!("Input file {:?} not available.", fname.as_ref());
        }
        let raw = read_to_string(fname)?;
        let mut input: Self = toml::from_str(&raw)?;

        if let Some(e) = input.efield.as_mut() {
            e.2.initialize(input.namdtime, input.nelm, input.dt);
        }

        Self::create_outputdir(&mut input.outdir)?;
        Ok(input)
    }


    fn efield_from_file<'de, D>(deserializer: D) -> Result<Option<(PathBuf, String, Efield)>, D::Error>
    where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        let path = PathBuf::from(s);
        Ok(Some((
                path.to_owned(),
                read_to_string(&path).unwrap(),
                Efield::from_file(path),
                )))
    }


    pub fn print_to_log(&self) {
        let input_print = format!("{}", self);
        let hashtag_line = "#".repeat(120);
        info!("Input file loaded. The formatted input is:\n\n{hashtag_line}\n{}\n{hashtag_line}\n\n", input_print);
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


    fn default_outputdir() -> PathBuf {
        PathBuf::from("output")
    }


    fn create_outputdir(dir: &mut PathBuf) -> Result<()> {
        if dir.is_file() {
            bail!("The output dir {:?} exists as a regular file, please change.", dir);
        }

        if dir.file_name().is_none() {
            bail!("The output dir {:?} cannot be current working dir, please change.", dir);
        }

        if dir.is_dir() {
            let parent = dir.parent().unwrap();
            let subdir = dir.file_name().unwrap().to_str().unwrap();
            let mut newdir: Option<PathBuf> = None;
            let mut tmpdir = PathBuf::new();

            for i in 1 ..= 99 {
                let dirstr = format!("{}_{:02}", &subdir, i);
                tmpdir = parent.join(&dirstr);
                if !tmpdir.is_file() && !tmpdir.is_dir() {
                    newdir = Some(tmpdir.clone());
                    break;
                }
            }

            if let Some(newdir) = newdir {
                warn!("The outdir {:?} is already exists and will be switched to {:?} for this run.", dir, newdir);
                *dir = newdir;
            } else {
                bail!("Existed outdir reached maximum homonymy outdirs: {:?}", tmpdir);
            }
        }

        info!("Log and output files will be stored in {:?} .", dir);
        create_dir_all(dir)?;

        Ok(())
    }
}


impl fmt::Display for Input {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "# NAMD-lumi input in toml format.")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}",   "rundir",     self.rundir)?;
        writeln!(f, " {:>20} = {}",     "ikpoint",    self.ikpoint)?;
        writeln!(f, " {:>20} = {:?}",   "brange",     self.brange)?;
        writeln!(f, " {:>20} = {:?}",   "basis_up",   self.basis_up)?;
        writeln!(f, " {:>20} = {:?}",   "basis_dn",   self.basis_dn)?;
        writeln!(f, " {:>20} = {}",     "nsw",        self.nsw)?;
        writeln!(f, " {:>20} = {}",     "ndigit",     self.ndigit)?;
        writeln!(f, " {:>20} = {}",     "namdtime",   self.namdtime)?;
        writeln!(f, " {:>20} = {}",     "dt",         self.dt)?;
        writeln!(f, " {:>20} = {}",     "nsample",    self.nsample)?;
        writeln!(f, " {:>20} = {}",     "ntraj",      self.ntraj)?;
        writeln!(f, " {:>20} = \"{}\"", "propmethod", self.propmethod)?;
        writeln!(f, " {:>20} = \"{}\"", "shmethod",   self.shmethod)?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {}", "nelm",         self.nelm)?;
        writeln!(f, " {:>20} = {}", "lreal",        self.lreal)?;
        writeln!(f, " {:>20} = {}", "lprint_input", self.lprint_input)?;
        writeln!(f, " {:>20} = {}", "lexcitation",  self.lexcitation)?;
        writeln!(f, " {:>20} = {}", "lreorder",     self.lreorder)?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "nacfname",    self.nacfname)?;
        writeln!(f, " {:>20} = {}",   "temperature", self.temperature)?;
        writeln!(f)?;

        if let Some(s) = self.scissor {
            writeln!(f, " {:>20} = {}", "scissor", s)?;
            writeln!(f)?;
        }

        if let Some(e) = self.efield.as_ref() {
            writeln!(f, " {:>20} = {:?}", "efield", e.0)?;
            writeln!(f, "## Content of {:?} :", e.0)?;
            for l in e.1.lines() {
                writeln!(f, "##     {}", l)?;
            }
            writeln!(f)?;
        }

        let nsample = self.inibands.len();
        writeln!(f, "# There are {} samples", nsample)?;
        let mut inibands = String::from("[");
        let mut inispins = String::from("[");
        let mut inisteps = String::from("[");
        for i in 0 .. nsample - 1 {
            inibands += &format!("{:6},", self.inibands[i]);
            inispins += &format!("{:6},", self.inispins[i]);
            inisteps += &format!("{:6},", self.inisteps[i]);
        }
        inibands += &format!("{:6}]", self.inibands[nsample - 1]);
        inispins += &format!("{:6}]", self.inispins[nsample - 1]);
        inisteps += &format!("{:6}]", self.inisteps[nsample - 1]);

        writeln!(f, " {:>10} = {}", "inibands", inibands)?;
        writeln!(f, " {:>10} = {}", "inispins", inispins)?;
        writeln!(f, " {:>10} = {}", "inisteps", inisteps)?;

        Ok(())
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
        efield       = "tests/1.5eV.rhai"
        scissor      = 1.14  # eV

        inibands     = [114, 114]
        inispins     = [1, 1]
        inisteps     = [514, 810]
            "#;

        let input: Input = toml::from_str(raw).unwrap();
        println!("{}", &input);
        assert_eq!(input.efield.unwrap().0, PathBuf::from("tests/1.5eV.rhai"));
        assert_eq!(input.scissor, Some(1.14f64));
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
        efield       = "tests/2.0eV.rhai"

        inibands     = [114, 114]
        inispins     = [1, 1]
        inisteps     = [514, 810]
            "#;

        let input: Input = toml::from_str(raw).unwrap();
        assert_eq!(input.efield.unwrap().0, PathBuf::from("tests/2.0eV.rhai"));
        assert_eq!(input.scissor, None);
    }
}
