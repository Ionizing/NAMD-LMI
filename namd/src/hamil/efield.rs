use std::fs::{
    File,
    read_to_string,
    write,
};
use std::fmt::Write;
use std::io::{
    Write as _,
    BufWriter,
};
use std::path::Path;
use std::sync::{
    OnceLock,
    Mutex,
};

use rhai::{
    AST,
    Dynamic,
    Engine,
    FLOAT,
    Scope,
};

use shared::{
    log,
    info,
    warn,
    anyhow::ensure,
    Context,
    Result,
};
use shared::MatX3;


pub struct Efield<'a> {
    raw: String,
    ast: AST,
    engine: Engine,
    scope: Scope<'a>,
}



#[allow(non_snake_case)]
pub fn EFIELD() -> &'static Mutex<Efield<'static>> {
    static EFIELD: OnceLock<Mutex<Efield>> = OnceLock::new();
    EFIELD.get_or_init(|| Mutex::new(Efield::from_str("").unwrap()))
}


impl Efield<'_> {
    pub fn from_str(raw: &str) -> Result<Self> {
        let raw = raw.to_string();

        let engine = Engine::new();
        let ast = engine.compile(&raw).context("Failed eval efield.")?;
        let mut scope = Scope::new();
        scope.push_constant("hbar", 0.6582119569);
        scope.push_constant("pi", std::f64::consts::PI);
        scope.push_constant("e", std::f64::consts::E);

        Ok(Self {
            raw,
            engine,
            ast,
            scope,
        })
    }


    pub fn from_file<P>(fname: P) -> Result<Self>
    where P: AsRef<Path> {
        let raw = read_to_string(fname)?;
        Self::from_str(&raw)
    }


    pub fn get_src(&self) -> &str {
        &self.raw
    }


    pub fn print_to_log(&self) {
        let mut s = String::new();
        let hashtag_line = "#".repeat(120);
        writeln!(s, "{hashtag_line}").unwrap();
        writeln!(s, "## Content of efield file:").unwrap();
        for l in self.raw.lines() {
            writeln!(s, "##     {}", l).unwrap();
        }
        writeln!(s, "{hashtag_line}").unwrap();
        log::info!("{}", s);
    }


    pub fn eval(&mut self, t: f64) -> [f64; 3] {
        let engine = &self.engine;
        let ast    = &self.ast;
        let scope  = &mut self.scope;

        let ret = engine.call_fn::<Dynamic>(scope, ast, "efield", (t,))
            .unwrap()
            .into_typed_array::<FLOAT>()
            .unwrap();
        assert_eq!(3, ret.len());
        
        [ret[0], ret[1], ret[2]]
    }


    pub fn eval_array(&mut self, ts: &[f64]) -> MatX3<f64> {
        let engine = &self.engine;
        let ast    = &self.ast;
        let scope  = &mut self.scope;

        ts.iter().cloned()
            .map(|t| -> [f64; 3] {
                let ret = engine.call_fn::<Dynamic>(scope, ast, "efield", (t,))
                    .unwrap()
                    .into_typed_array::<FLOAT>()
                    .unwrap();
                [ret[0], ret[1], ret[3]]
            })
            .collect()
    }


    pub fn get_eafield_array(&mut self, namdtime: usize, potim: f64, nelm: usize) -> (Vec<f64>, [MatX3<f64>; 2]) {
        let edt = potim / nelm as f64;
        let mut tt = vec![0.0; 0];

        for iion in 0 .. namdtime {
            let t0 = potim * iion as f64;
            for _ielm in 0 .. nelm {
                let t = t0 + edt * nelm as f64;
                tt.push(t);
            }
        }

        let efield = self.eval_array(&tt);
        let mut afield = MatX3::<f64>::new();

        // integrate electric field into vector potential
        // E = - \partial A / \partial t
        // A = - \int dt E
        let mut at = [0.0; 3];
        for e in efield.iter().copied() {
            at = [at[0] - e[0] * edt,
                  at[1] - e[1] * edt,
                  at[2] - e[2] * edt];
            afield.push(at);
        }

        (tt, [efield, afield])
    }


    pub fn print_eafield_tofile<P>(&mut self, dir: P, namdtime: usize, potim: f64, nelm: usize) -> Result<()>
    where P: AsRef<Path> {
        let dir = dir.as_ref().to_owned();
        ensure!(dir.is_dir(), "The parameter dir should be a valid directory.");
        let eafield_fname = dir.with_file_name("EAFIELD.txt");

        let (t, [efield, afield]) = self.get_eafield_array(namdtime, potim, nelm);

        info!("Writing electric field to {:?} ...", &eafield_fname);
        if eafield_fname.is_file() {
            warn!("File {:?} exists, overwriting ...", &eafield_fname);
        }

        let mut f = BufWriter::new(
            File::options().write(true)
                           .create(true)
                           .truncate(true)
                           .open(&eafield_fname)?
            );

        writeln!(f, "# Time(fs) Ex Ey Ez(V/Å) Ax Ay Az(V*fs/Å)")?;
        for i in 0 .. t.len() {
            let e = efield[i];
            let a = afield[i];
            writeln!(f, "{:15.5} {:15.7} {:15.7} {:15.7}    {:15.7} {:15.7} {:15.7}",
                t[i], e[0], e[1], e[2], a[0], a[1], a[2])?;
        }

        Ok(())
    }


    pub fn singleton_from_str(raw: &str) -> Result<&'static Mutex<Self>> {
        let instance = Efield::from_str(raw)?;
        *(EFIELD().lock().unwrap()) = instance;
        Ok(&EFIELD())
    }


    pub fn singleton_from_file<P>(fname: P) -> Result<&'static Mutex<Self>>
    where P: AsRef<Path> {
        let instance = Efield::from_file(fname)?;
        *(EFIELD().lock().unwrap()) = instance;
        Ok(&EFIELD())
    }


    pub fn template_to_file<P>(fname: P) -> Result<()>
    where P: AsRef<Path> {
        write(fname, include_str!("./efield_template.rhai"))?;
        Ok(())
    }
}
