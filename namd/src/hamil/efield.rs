use std::fs::{
    File,
    read_to_string,
};
use std::io::{
    Write as _,
    BufWriter,
};
use std::path::Path;
use std::sync::OnceLock;

use rhai::{
    AST,
    Dynamic,
    Engine,
    FLOAT,
    Scope,
};

use shared::{
    info,
    warn,
    anyhow::ensure,
    Context,
    Result,
};
use shared::MatX3;


const ENGINE: OnceLock<Engine> = <_>::default();
fn get_engine() -> &'static Engine {
    ENGINE.get_or_init(|| Engine::new())
}


#[derive(Clone)]
pub struct Efield<'a> {
    raw: String,
    ast: AST,
    scope: Scope<'a>,
}


impl<'a> Efield<'a> {
    pub fn from_str(raw: &str) -> Result<Self> {
        let raw = raw.to_string();

        let engine = get_engine();
        let ast = engine.compile(&raw).context("Failed eval efield.")?;
        let mut scope = Scope::new();
        scope.push_constant("hbar", 0.6582119569);
        scope.push_constant("pi", std::f64::consts::PI);
        scope.push_constant("e", std::f64::consts::E);

        Ok(Self {
            raw,
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


    pub fn eval(&mut self, t: f64) -> [f64; 3] {
        let engine = get_engine();
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
        let engine = get_engine();
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
        let mut efield = MatX3::<f64>::new();
        let mut afield = MatX3::<f64>::new();

        for iion in 0 .. namdtime {
            let t0 = potim * iion as f64;
            for ielm in 0 .. nelm {
                let t = t0 + edt * nelm as f64;
                tt.push(t);
            }
        }

        efield = self.eval_array(&tt);

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


    pub fn print_efield_tofile<P>(&self, dir: P, namdtime: usize, potim: f64, nelm: usize) -> Result<()>
    where P: AsRef<Path> {
        let dir = dir.as_ref();
        ensure!(dir.is_dir(), "The parameter dir should be a valid directory.");
        let efield_fname = dir.with_file_name("EFIELD.txt");
        let afield_fname = dir.with_file_name("AFIELD.txt");

        let (t, [efield, afield]) = self.get_eafield_array(namdtime, potim, nelm);


        {
            info!("Writing electric field to {:?} ...", efield_fname);
            if efield_fname.is_file() {
                warn!("File {:?} exists, overwriting ...", efield_fname);
            }

            let mut f = BufWriter::new(
                File::options().write(true)
                               .create(true)
                               .truncate(true)
                               .open(efield_fname)?
                );

            writeln!(f, "# Time(fs) Ex Ey Ez(V/Å)")?;
            for (i, e) in efield.iter().enumerate() {
                writeln!(f, "{:15.5} {:15.7} {:15.7} {:15.7}", t[i], e[0], e[1], e[2])?;
            }
        }


        {
            info!("Writing integrated vector potential to {:?} ...", afield_fname);
            if efield_fname.is_file() {
                warn!("File {:?} exists, overwriting ...", afield_fname);
            }

            let mut f = BufWriter::new(
                File::options().write(true)
                               .create(true)
                               .truncate(true)
                               .open(afield_fname)?
                );

            writeln!(f, "# Time(fs) Ax Ay Az(V*fs/Å)")?;
            for (i, a) in afield.iter().enumerate() {
                writeln!(f, "{:15.5} {:15.7} {:15.7} {:15.7}", t[i], a[0], a[1], a[2])?;
            }
        }

        Ok(())
    }
}
