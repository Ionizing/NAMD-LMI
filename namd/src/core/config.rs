use std::path::Path;
use std::fmt;
use serde::Deserialize;
use shared::Result;

pub trait NamdConfig: Clone + Sized + Default + fmt::Display + for<'a> Deserialize<'a> {
    fn from_file<P>(fname: P) -> Result<Self>
        where P: AsRef<Path>,
              Self: Sized;
    fn to_file<P>(&self, fname: P) -> Result<()>
        where P: AsRef<Path>;
}
