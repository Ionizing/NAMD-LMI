use std::fs;
use std::path::Path;

pub use anyhow::{
    self,
    Result,
    Context,
    bail,
};
pub use ndarray::{
    self,
    Array1,
    Array2,
    Array3,
};
pub use ndarray_linalg;
pub use regex::Regex;
pub use clap::Parser;
pub use tracing::{
    self,
    warn,
    info,
    debug,
    error,
};

#[allow(non_camel_case_types)]
pub type c64 = ndrustfft::Complex<f64>;
#[allow(non_camel_case_types)]
pub type c32 = ndrustfft::Complex<f32>;

// Convenient types
pub type Vector<T> = Array1<T>;  // Define this type to use broadcast operations.
pub type Matrix<T> = Array2<T>;
pub type Cube<T>   = Array3<T>;
pub type MatX3<T> = Vec<[T;3]>;  // Nx3 matrix
pub type Mat33<T> = [[T;3];3];   // 3x3 matrix

pub trait OptProcess : Parser {
    fn process(&self) -> Result<()>;
}

/// Index array containing negative indices => Index array full of positive indices.
/// `-1` means the last index,
/// If `v` contains `0`, selecting the total indices,  `1..=len` is returned.
pub fn index_transform(v: Vec<i32>, len: usize) -> Vec<usize> {
    if v.contains(&0) {
        (1 ..= len).collect()
    } else {
        v.into_iter()
         .map(|i| {
            if i < 0 {
                i.rem_euclid(len as i32) as usize + 1
            } else {
                i as usize
            }
         })
        .collect()
    }
}


/// Parse string containing range and integers into `Vec<i32>`
/// 
/// Valid strings can be `"1..5 12 -1..3 0"` will be parsed as
/// `vec![1, 2, 3, 4, 5, 12, -1, 1, 2, 3]`, `0` is filtered out.
/// If only `0` is supplied, `Ok(vec![0])` is return
pub fn range_parse(input: &str) -> Result<Vec<i32>> {
    if input.trim() == "0" {
        return Ok(vec![0]);
    }

    let mut ret = vec![];

    let re_range = Regex::new(r"^(-?\d+)\.\.(-?\d+)$").unwrap();
    let re_digit = Regex::new(r"^-?\d+$").unwrap();

    for s in input.split_ascii_whitespace() {
        if re_digit.is_match(s) {
            ret.push(s.parse::<i32>().unwrap());
        } else if re_range.is_match(s) {
            let m = re_range.captures(s).unwrap();
            let start = m.get(1).unwrap().as_str().parse::<i32>().unwrap();
            let end   = m.get(2).unwrap().as_str().parse::<i32>().unwrap();
            
            if start > end {
                bail!("[RANGE_PARSE]: start is greater than end in token \'{}\'", s);
            }

            let to_be_extend = (start ..= end).collect::<Vec<_>>();
            ret.extend( to_be_extend );
        } else {
            bail!("[RANGE_PARSE]: token \'{}\' is invalid, cannot be parsed as range or integer", s);
        }
    }

    let ret = ret.into_iter().filter(|x| *x != 0).collect::<Vec<_>>();

    Ok(ret)
}


/// Copies the contents of one file to another. Permission bits are preserved.
///
/// If `to` will be overwritten if it exists. If `to` does not exists, a file with same name of
/// `from` will be created and written.
pub fn copy_file_to<P, Q>(from: P, to: Q) -> Result<u64>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    assert!(from.as_ref().is_file());

    if to.as_ref().is_file() {
        return Ok(fs::copy(from, to)?);
    }

    let fname = from.as_ref().file_name().unwrap();
    let target = to.as_ref().join(&fname);
    return Ok(fs::copy(from, target)?);
}


pub fn link_file_to<P, Q>(from: P, to: Q) -> Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    assert!(from.as_ref().is_file());
    assert!(to.as_ref().is_dir());

    let from = from.as_ref().canonicalize().unwrap();
    let fname = from.file_name().unwrap();
    let target = to.as_ref().join(&fname);

    #[cfg(target_os = "windows")]
    std::os::winodws::fs::symlink_file(from, &target)?;

    #[cfg(not(target_os = "windows"))]
    std::os::unix::fs::symlink(from, &target)?;

    Ok(())
}


#[derive(Clone)]
pub struct Structure {
    pub cell          : Mat33<f64>,
    pub ion_types     : Vec<String>,
    pub ions_per_type : Vec<i32>,
    pub car_pos       : MatX3<f64>,
    pub frac_pos      : MatX3<f64>,
    pub constr        : Option<MatX3<bool>>,
}


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_index_transform() {
        assert_eq!(index_transform(vec![-1], 5), vec![5]);
        assert_eq!(index_transform(vec![-1, 5], 5), vec![5, 5]);
        assert_eq!(index_transform(vec![-1, 0], 5), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_range_parse() {
        assert_eq!(range_parse("1..6").unwrap(), vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(range_parse("9..12").unwrap(), vec![9, 10, 11, 12]);
        assert_eq!(range_parse("-1..5").unwrap(), vec![-1, 1, 2, 3, 4, 5]);
        assert_eq!(range_parse("-10..-5").unwrap(), vec![-10, -9, -8, -7, -6, -5]);
        assert_eq!(range_parse("-10..-9 4    29 \n  -5").unwrap(), vec![-10, -9, 4, 29, -5]);
        assert_eq!(range_parse(" 0   ").unwrap(), vec![0]);
        assert!(range_parse("5..-1").is_err());
        assert!(range_parse("..-1").is_err());
        assert!(range_parse("..10").is_err());
        assert!(range_parse("1 .. 10").is_err());
        assert!(range_parse("1 .. 10").is_err());
        assert!(range_parse("1-2..5").is_err());
    }
}
