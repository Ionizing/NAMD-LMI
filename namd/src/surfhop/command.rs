use std::path::PathBuf;
use std::fs::create_dir_all;

use rayon::prelude::*;

use clap::{Parser, ValueEnum};
use shared::{
    Result,
    bail,
    log,
    copy_file_to,
    link_file_to,
};
use crate::OptProcess;

use crate::core::{
    SurfaceHopping,
    NamdConfig,
};
use crate::surfhop;


#[derive(Debug, Parser)]
/// Perform the surface-hopping process with given Hamiltonian file and config file.
#[command(arg_required_else_help(true))]
pub struct SurfhopCommand {
    #[arg(short='n', long, default_value_t=0)]
    /// Number of threads for parallel calculation.
    /// 
    /// If 0 is set, it will fall back to the number of logic CPU cores of you machine.
    nthreads: usize,

    #[arg(short='c', long, default_value="surfhop_config.toml", aliases=["cfg", "conf"])]
    /// Config file name.
    ///
    /// Aliases: "cfg", "conf".
    config: PathBuf,

    #[arg(long, value_enum, alias="gen")]
    /// Generate auxiliary files for the calculation and analysis.
    ///
    /// The surface-hopping will not run if this flag is set.
    ///
    /// Alias: "gen".
    generate: Option<TemplateGenerator>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum TemplateGenerator {
    #[value(aliases=["config", "cfg", "conf"])]
    /// Generate config template for Hamiltonian generation.
    /// Aliases: "config", "cfg" and "conf".
    ConfigTemplate,

    #[value(name="inistep.py", aliases=["inistep", "ini"])]
    /// Generate Python script to help append `inistep` field.
    /// Aliases: "inistep" and "ini"
    Inistep,

    #[value(aliases=["post-process", "postprocess", "pp"])]
    /// Generate post-process scripts for surface-hopping analysis.
    /// Aliases: "post-process", "postprocess", "pp".
    PostprocessTemplate,
}


impl OptProcess for SurfhopCommand {
    fn process(&self) -> Result<()> {
        use TemplateGenerator::*;

        if let Some(g) = self.generate {
            return match g {
                ConfigTemplate => {
                    {
                        log::info!("writing `surfhop_config_template.toml` ...");
                        surfhop::SurfhopConfig::default().to_file("surfhop_config_template.toml")
                    }.and_then(|_| {
                        log::info!("Writing `inisteps.py` ...");
                        surfhop::SurfhopConfig::write_inistep_py("inisteps.py")
                    })
                },
                Inistep => {
                    log::info!("Writing `inisteps.py` ...");
                    surfhop::SurfhopConfig::write_inistep_py("inisteps.py")
                },
                PostprocessTemplate => todo!(),
            }
        }

        // Start running surface-hopping method.
        log::info!("Prepare to run surface-hopping method ...");
        rayon::ThreadPoolBuilder::new().num_threads(self.nthreads).build_global().unwrap();

        let mut cfg = surfhop::SurfhopConfig::from_file(&self.config)?;
        create_outputdir(cfg.get_outdir_mut())?;
        let cfg = cfg;      // cancel mutability

        crate::logging::logger_redirect(&cfg.get_outdir())?;
        log::info!("Got Surface Hopping config:\n{}", &cfg);
        let sh = surfhop::Surfhop::from_config(&cfg)?;

        sh.into_par_iter()
            .map(|mut v| v.run())
            .collect::<Result<Vec<()>>>()?;

        Ok(())
    }
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
            log::warn!("The outdir {:?} is already exists and will be switched to {:?} for this run.", dir, newdir);
            *dir = newdir;
        } else {
            bail!("Existed outdir reached maximum homonymy outdirs: {:?}", tmpdir);
        }
    }

    log::info!("Log and output files will be stored in {:?} .", dir);
    create_dir_all(dir)?;

    Ok(())
}
