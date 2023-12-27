use std::fs::File;
use std::io::{
    Write as _,
    BufWriter,
};
use std::path::{
    Path,
    PathBuf,
};

use shared::{
    info,
    warn,
    MatX3,
    Result,
};
use rhai::{
    AST,
    Dynamic,
    Engine,
    FLOAT,
    Scope,
};


#[derive(Clone)]
pub struct Efield {
    ast: AST,
    namdtime: usize,
    nelm: usize,
    dt: f64,
    t: Vec<f64>,
    efield: MatX3<f64>,
    afield: MatX3<f64>,
}


impl Efield {
    pub fn from_file(fname: PathBuf) -> Self
    {
        let engine = Engine::new();
        let ast = engine.compile_file(fname).unwrap();

        Self {
            ast,
            namdtime: 0,
            nelm: 0,
            dt: 0.0,
            t: vec![],
            efield: vec![],
            afield: vec![],
        }
    }


    pub fn initialize(&mut self, namdtime: usize, nelm: usize, dt: f64) {
        if self.namdtime == namdtime  &&
           self.nelm == nelm &&
           self.dt == dt   &&
           self.t.len() == namdtime * nelm &&
           self.efield.len() == namdtime * nelm &&
           self.afield.len() == namdtime * nelm {
            return;
        }

        self.namdtime  = namdtime;
        self.nelm = nelm;
        self.dt   = dt;

        let engine = Engine::new();
        let mut scope = Scope::new();
        scope.push_constant("hbar", 0.6582119569);

        let mut f = |t: f64| {
            let ret = engine.call_fn::<Dynamic>(&mut scope, &self.ast, "efield", (t,))
                .unwrap()
                .into_typed_array::<FLOAT>()
                .unwrap();
            assert_eq!(3, ret.len());
            
            [ret[0], ret[1], ret[2]]
        };

        let mut tt = vec![0.0; 0];
        let mut efield = vec![[0.0; 3]; 0];
        let mut afield = vec![[0.0; 3]; 0];

        let edt = dt / nelm as f64;

        // get electric field
        for iion in 0 .. namdtime {
            let t0 = dt * iion as f64;
            for ielm in 0 .. nelm {
                let t = t0 + edt * ielm as f64;

                tt.push(t);
                efield.push(f(t));
            }
        }

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

        self.t      = tt;
        self.efield = efield;
        self.afield = afield;
    }


    /// Get electric field E, in V/Å
    pub fn get_efield(&self, iion: usize, ielm: usize) -> [f64; 3] {
        let i = iion * self.nelm + ielm;
        self.efield[i]
    }


    /// Get vector potential A, in V*fs/Å
    pub fn get_afield(&self, iion: usize, ielm: usize) -> [f64; 3] {
        let i = iion * self.nelm + ielm;
        self.afield[i]
    }


    pub fn print_efield_to_file<P>(&self, fname: &P) -> Result<()>
    where
        P: AsRef<Path> + ?Sized,
    {
        info!("Writing electric field to {:?} ...", fname.as_ref());
        if fname.as_ref().is_file() {
            warn!("File {:?} already exists, overwriting ...", fname.as_ref());
        }

        let mut f = BufWriter::new(
            File::options().write(true)
                           .create(true)
                           .truncate(true)
                           .open(fname)?
            );

        writeln!(f, "# Time(fs) Ex Ey Ez(V/Å)")?;
        
        for (i, e) in self.efield.iter().enumerate() {
            writeln!(f, "{:15.5} {:15.7} {:15.7} {:15.7}", self.t[i], e[0], e[1], e[2])?;
        }

        Ok(())
    }


    pub fn print_afield_to_file<P>(&self, fname: &P) -> Result<()>
    where
        P: AsRef<Path> + ?Sized,
    {
        info!("Writing integrated vector potential to {:?} ...", fname.as_ref());
        if fname.as_ref().is_file() {
            warn!("File {:?} already exists, overwriting ...", fname.as_ref());
        }

        let mut f = BufWriter::new(
            File::options().write(true)
                           .create(true)
                           .truncate(true)
                           .open(fname)?
            );

        writeln!(f, "# Time(fs) Ax Ay Az(V*fs/Å)")?;
        
        for (i, a) in self.afield.iter().enumerate() {
            writeln!(f, "{:15.5} {:15.7} {:15.7} {:15.7}", self.t[i], a[0], a[1], a[2])?;
        }

        Ok(())
    }
}
