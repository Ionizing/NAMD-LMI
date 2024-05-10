use std::fs;
use std::path::Path;
use shared::{
    Result,
    info,
};


pub const EFIELD_EXAMPLE:  &'static str = include_str!("./efield.rhai");
pub const INIBANDS_PY:     &'static str = include_str!("./inibands.py");
pub const POST_PROCESS_PY: &'static str = include_str!("./post-process.py");


fn write_file_to_dir<P>(dir: P, fname: &str, content: &str) -> Result<()>
where
    P: AsRef<Path>,
{
    assert!(dir.as_ref().is_dir());

    let target = dir.as_ref().join(fname);
    if target.is_file() {
        info!("Script {:?} exists, skip writing.", target);
        return Ok(());
    }

    fs::write(target, content)?;
    Ok(())
}


pub enum Scripts {
    EField,
    IniBands,
    PostProcess,
}


pub fn write_script<P>(dir: P, script: Scripts) -> Result<()>
where
    P: AsRef<Path>,
{
    use Scripts::*;

    let (fname, content) = match script {
        EField      => ("efield.rhai",     EFIELD_EXAMPLE),
        IniBands    => ("inibands.py",     INIBANDS_PY),
        PostProcess => ("post-process.py", POST_PROCESS_PY),
    };

    write_file_to_dir(dir, fname, content)
}
