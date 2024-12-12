use std::fs;
use std::fmt;
use std::path::{Path, PathBuf};
use std::f64::consts::PI;

use serde::{de::Error, Deserialize, Deserializer};
use toml;
use regex::Regex;
use itertools::Itertools;
use shared::{
    log,
    Result,
    anyhow::ensure,
};

use crate::core::NamdConfig;
use crate::cli::write_script;


#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
pub enum SHMethod {
    #[serde(alias="fssh", alias="FSSH")]
    FSSH,
    #[serde(alias="dish", alias="DISH")]
    DISH,
    #[serde(alias="dcsh", alias="DCSH")]
    DCSH,
}

const INISTEPS_PY: &str = include_str!("./inisteps.py");


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


#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
pub enum SmearingMethod {
    #[serde(alias="Gaussian", alias="gaussian")]
    GaussianSmearing,

    #[serde(alias="Lorentzian", alias="lorentzian")]
    LorentzianSmearing,
}


impl SmearingMethod {
    pub fn apply_smearing<I1, I2, I3, R>(&self, x: I1, centers: I2, width: f64, scales: Option<I3>) -> R
    where I1: AsRef<[f64]>,
          I2: AsRef<[f64]>,
          I3: AsRef<[f64]>,
          R: FromIterator<f64>,
    {
        //let centers_len = centers.as_ref()
        let x = x.as_ref();
        let centers = centers.as_ref();
        let centers_len = centers.len();
        let scales = scales
            .map(|x| x.as_ref().to_owned())
            .unwrap_or(vec![1.0; centers_len]);

        match self {
            Self::GaussianSmearing   => Self::gaussian(x, centers, width, &scales),
            Self::LorentzianSmearing => Self::lorentzian(x, centers, width, &scales),
        }
    }


    /// f(x,μ,σ) = exp(-(x-μ)^2 / (2*σ^2) / (σ*sqrt(2π))
    fn gaussian<I1, I2, I3, R>(x: I1, mus: I2, sigma: f64, scales: I3) -> R
    where I1: AsRef<[f64]>,
          I2: AsRef<[f64]>,
          I3: AsRef<[f64]>,
          R: FromIterator<f64>,
    {
        let x   = x.as_ref();
        let mus = mus.as_ref();
        let scales = scales.as_ref();

        let inv_two_sgm_sqr = 1.0 / (2.0 * sigma.powi(2));           // 1.0/(2*σ^2)
        let inv_sgm_sqrt2pi = 1.0 / (sigma * (2.0 * PI).sqrt());     // 1.0/(σ*sqrt(2π))

        x.iter().cloned()
            .map(|x| { mus.iter().cloned().zip(scales.iter().cloned())
                .map(|(c, s)| {
                    -((x - c).powi(2) * inv_two_sgm_sqr).exp() * inv_sgm_sqrt2pi * s
                })
                .sum()
            })
            .collect()
    }


    /// lorentz_smearing(x::AbstractArray, x0::Float64, Γ=0.05) = @. Γ/(2π) / ((x-x0)^2 + (Γ/2)^2)
    fn lorentzian<I1, I2, I3, R>(x: I1, x0s: I2, gamma: f64, scales: I3) -> R
    where I1: AsRef<[f64]>,
          I2: AsRef<[f64]>,
          I3: AsRef<[f64]>,
          R: FromIterator<f64>,
    {
        let x   = x.as_ref();
        let x0s = x0s.as_ref();
        let scales = scales.as_ref();
        let gam_div_2pi  = gamma / (2.0 * PI);      // Γ/(2π)
        let gam_half_sqr = (gamma / 2.0).powi(2);   // (Γ/2)^2

        x.iter().cloned()
            .map(|x| { x0s.iter().cloned().zip(scales.iter().cloned())
                .map(|(c, s)| {
                    gam_div_2pi / ((x - c).powi(2) + gam_half_sqr) * s
                })
                .sum()
            })
            .collect()
    }
}


#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
/// Maintain the detailed balance?
///
/// Always: The detailed balance always works, meaning that the upward hoppings would
///     be reduced by the Boltzmann factor P = P0 * exp( -(Ei-Ej) /kBT ), Ei > Ej;
///
/// DependsOnEField: When the external electric field norm |E| > 0, the detailed balance
///     is nolonger maintained since there is external energy input, and the upward
///     hoppings are allowed; When the |E| = 0, the detailed balance works.
///
/// Never: The detailed balances will never be maintained, allowing all upward hoppings.
pub enum DetailedBalance {
    #[serde(alias="always", alias="Always")]
    Always,

    #[serde(alias="depends-on-efield", alias="DependsOnEField")]
    DependsOnEField,

    #[serde(alias="nac-only", alias="NacOnly")]
    NacOnly,

    #[serde(alias="never", alias="Never")]
    Never,
}


impl Default for DetailedBalance {
    fn default() -> Self {
        DetailedBalance::DependsOnEField
    }
}


impl fmt::Display for DetailedBalance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use DetailedBalance as DB;
        let s: &str = match self {
            DB::Always => "Always",
            DB::DependsOnEField => "DependsOnEField",
            DB::NacOnly => "NacOnly",
            DB::Never => "Never",
        };
        f.write_str(s)
    }
}


#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SurfhopConfig {
    hamil_fname: PathBuf,
    namdtime: usize,
    nelm: usize,
    ntraj: usize,
    shmethod: SHMethod,
    outdir: PathBuf,

    #[serde(default = "SurfhopConfig::default_detailed_balance")]
    detailed_balance: DetailedBalance,

    /// Which type of smearing to use:
    ///
    /// - GaussianSmearing: `f(x) = exp(-(x-μ)^2 / (2*σ^2)) / (σ*sqrt(2π))` \n{}
    /// - LorentzianSmearing: `f(x) = Γ/(2π) / ((x-μ)^2 + (Γ/2)^2)`
    #[serde(default = "SurfhopConfig::default_smearing_method")]
    smearing_method: SmearingMethod,

    /// Smearing linewidth, in eV
    #[serde(default = "SurfhopConfig::default_smearing_sigma")]
    smearing_sigma: f64,

    /// Number of points PER eV in photon spectra.
    #[serde(rename = "smearing_npoints_per_eV",
            default = "SurfhopConfig::default_smearing_npoints_per_ev")]
    smearing_npoints_per_ev: usize,

    #[serde(deserialize_with="SurfhopConfig::parse_iniband")]
    iniband: Vec<i32>,
    inisteps: Vec<usize>,
}


impl SurfhopConfig {
    fn default_detailed_balance() -> DetailedBalance { DetailedBalance::default() }
    fn default_smearing_method() -> SmearingMethod { SmearingMethod::LorentzianSmearing }
    fn default_smearing_sigma() -> f64 { 0.01 }
    fn default_smearing_npoints_per_ev() -> usize { 500 }

    fn parse_iniband<'de, D>(deserializer: D) -> std::result::Result<Vec<i32>, D::Error>
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
                        format!("Invalid integer '{}' from `iniband`, 0 is not allowed.", s)
                    ));
                }
                ret.push(s.parse().unwrap())
            } else if re_range.is_match(s) {
                let m = re_range.captures(s).unwrap();
                let start = m.get(1).unwrap().as_str().parse::<i32>().unwrap();
                let end   = m.get(2).unwrap().as_str().parse::<i32>().unwrap();

                if start * end <= 0 {
                    return Err(D::Error::custom(
                        format!("Invalid range '{}' from `iniband`, start and end must have save sign.", s)
                    ));
                }

                let sign = start.signum();
                let to_be_extend = (start.abs() ..= end.abs()).map(|x| x*sign).collect::<Vec<_>>();
                ret.extend(to_be_extend);
            } else {
                return Err(D::Error::custom(
                    format!("Invalid token '{}' from `iniband`, it should be either range (start..end) or integer.", s)
                ));
            }
        }

        // check band occupations, if any larger than 2, return error
        let mut err = "Invalid inibands, some bands have more than 2 electrons: ".to_string();
        let mut ret_dedup = ret.clone();
        ret_dedup.sort();
        let mut fail = false;
        for cnt in ret_dedup.into_iter().dedup_with_count() {
            if cnt.0 > 2 {
                fail = true;
                err.push_str(&format!(" occ({})={}", cnt.1, cnt.0));
            }
        }
        if fail {
            err.push_str(&format!(", please check."));
            return Err(D::Error::custom(err));
        }

        if ret.is_empty() {
            return Err(D::Error::custom("Empty `inibands` is not allowed."));
        }

        return Ok(ret)
    }

    pub fn get_hamil_fname(&self) -> &PathBuf { &self.hamil_fname }
    pub fn get_namdtime(&self) -> usize { self.namdtime }
    pub fn get_nelm(&self) -> usize { self.nelm }
    pub fn get_ntraj(&self) -> usize { self.ntraj }
    pub fn get_shmethod(&self) -> SHMethod { self.shmethod }
    pub fn get_outdir(&self) -> &PathBuf { &self.outdir }
    pub fn get_outdir_mut(&mut self) -> &mut PathBuf { &mut self.outdir }
    pub fn get_detailed_balance(&self) -> DetailedBalance { self.detailed_balance }
    pub fn get_smearing_method(&self) -> SmearingMethod { self.smearing_method }
    pub fn get_smearing_sigma(&self) -> f64 { self.smearing_sigma }
    pub fn get_npoints_per_ev(&self) -> usize { self.smearing_npoints_per_ev }

    pub fn get_iniband(&self) -> &[i32] { &self.iniband }
    pub fn get_inisteps(&self) -> &[usize] { &self.inisteps }

    pub fn print_to_log(&self) {
        let config_print = format!("{}", self);
        let hashtag_line = "#".repeat(120);
        log::info!("SurfhopConfig file loaded. The formatted config is:\n\n{hashtag_line}\n{}{hashtag_line}\n\n", config_print);
    }

    pub fn write_inistep_py<P>(fname: P) -> Result<()>
    where P: AsRef<Path> {
        write_script(fname, INISTEPS_PY, true)
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
            detailed_balance: DetailedBalance::default(),
            smearing_method: SmearingMethod::LorentzianSmearing,
            smearing_sigma: 0.01,
            smearing_npoints_per_ev: 500,

            iniband: vec![0],
            inisteps: vec![1, 2, 3],
        }
    }
}


impl fmt::Display for SurfhopConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "####             NAMD-lmi config for surface-hopping calculation            ####")?;
        writeln!(f, "#### YOUT NEED TO CHANGE THE PARAMETERS IN THE FOLLOWING TO FIT YOUR SYSTEM ####")?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "hamil_fname", self.hamil_fname)?;
        writeln!(f, " {:>20} = {:?}", "namdtime", self.namdtime)?;
        writeln!(f, " {:>20} = {:?}", "nelm", self.nelm)?;
        writeln!(f, " {:>20} = {:?}", "ntraj", self.ntraj)?;
        writeln!(f, " {:>20} = {:#}", "shmethod", self.shmethod)?;
        writeln!(f, " {:>20} = {:?}", "outdir", self.outdir)?;
        writeln!(f, " {:>20} = \"{:?}\"", "detailed_balance", self.detailed_balance)?;
        writeln!(f, " {:>20} = \"{:?}\"", "smearing_method", self.smearing_method)?;
        writeln!(f, " {:>20} = {:?}", "smearing_sigma", self.smearing_sigma)?;
        writeln!(f, " {:>20} = {:?}", "smearing_npoints_per_eV", self.smearing_npoints_per_ev)?;
        writeln!(f)?;

        writeln!(f, " {:>20} = {:?}", "iniband", self.iniband)?;
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
        detailed_balance = "never"
        smearing_method = "gaussian"
        smearing_sigma = 0.05
        smearing_npoints_per_eV = 1000

        iniband = "-1..-4 1..4"
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
            detailed_balance: DetailedBalance::Never,
            smearing_method: SmearingMethod::GaussianSmearing,
            smearing_sigma: 0.05,
            smearing_npoints_per_ev: 1000,

            iniband: vec![-1, -2, -3, -4, 1, 2, 3, 4],
            inisteps: vec![114, 514],
        };

        assert_eq!(expect_cfg, actual_cfg);

        println!("{}", &expect_cfg);
    }
}
