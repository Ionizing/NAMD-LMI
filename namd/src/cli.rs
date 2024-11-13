use std::path::Path;
use std::fs;
use std::io::Write as _;
use std::sync::OnceLock;

use clap::{
    Parser,
    builder::styling::{
        AnsiColor,
        Effects,
        Styles,
    },
};

use shared::{log, Result};
use crate::version::Version;
use crate::logging::logger_init;


pub fn get_style() -> Styles {
    static INSTANCE: OnceLock<Styles> = OnceLock::new();
    INSTANCE.get_or_init(|| {
        Styles::styled()
            .header(AnsiColor::Yellow.on_default() | Effects::BOLD)
            .usage(AnsiColor::Green.on_default()   | Effects::BOLD)
            .literal(AnsiColor::Green.on_default() | Effects::BOLD)
            .placeholder(AnsiColor::BrightBlue.on_default())
            .error(AnsiColor::BrightRed.on_default())
            .valid(AnsiColor::BrightYellow.on_default())
    }).to_owned()
}


pub trait OptProcess : Parser {
    fn process(&self) -> Result<()>;
}


#[derive(Debug, Parser)]
#[command(name = "NAMD-LMI",
          about = Version::new().to_string(),
          long_about = format!("{:#}", Version::new()),
          version,
          author = "@Ionizing github.com/Ionizing/NAMD-LMI",
          styles = get_style())]
enum Opt {
    Nac(crate::nac::NacCommand),
    Hamil(crate::hamil::HamilCommand),
    Surfhop(crate::surfhop::SurfhopCommand),
}


impl OptProcess for Opt {
    fn process(&self) -> Result<()> {
        use Opt::*;
        
        logger_init();
        log::info!("Global logger initialized with targets being stderr and \"./globalrun.log\"");

        match self {
            Nac(cmd) => cmd.process(),
            Hamil(cmd) => cmd.process(),
            Surfhop(cmd) => cmd.process(),
        }
    }
}


pub fn run() -> Result<()> {
    Opt::parse().process()
}


pub fn write_script<P>(fname: P, content: &str, with_exe_permission: bool) -> Result<()>
where P: AsRef<Path> {
    let mut f = fs::File::create(fname)?;
    f.write(content.as_bytes())?;

    #[cfg(unix)]
    if with_exe_permission {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = f.metadata()?.permissions();
        perms.set_mode(0o755);
        f.set_permissions(perms)?;
    }

    Ok(())
}
