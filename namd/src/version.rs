use std::fmt;
use once_cell::sync::OnceCell;

pub mod built_info {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}


fn built_time() -> &'static str {
    static INSTANCE: OnceCell<String> = OnceCell::new();
    INSTANCE.get_or_init(|| {
        built::util::strptime(built_info::BUILT_TIME_UTC)
            .with_timezone(&built::chrono::offset::Local)
            .to_string()
    })
    .as_str()
}


#[derive(Debug)]
pub struct Version<'a> {
    // https://www.asciiart.eu/text-to-ascii-art
    name:           &'a str,
    logo:           &'a str,
    version_str:    &'a str,
    authors:        &'a str,
    built_time:     &'a str,
    git_hash_long:  Option<&'a str>,
    git_hash_short: Option<&'a str>,
    git_dirty:      Option<bool>,

    host:           &'a str,
    opt_level:      &'a str,
    profile:        &'a str,
    rustc:          &'a str,
    target:         &'a str,
}


impl<'a> Version<'a> {
    pub fn new() -> Self {
        Self {
            name:           built_info::PKG_NAME,
            logo:           include_str!("./logo.txt"),
            version_str:    built_info::PKG_VERSION,
            authors:        built_info::PKG_AUTHORS,
            built_time:     built_time(),
            git_hash_long:  built_info::GIT_COMMIT_HASH,
            git_hash_short: built_info::GIT_COMMIT_HASH_SHORT,
            git_dirty:      built_info::GIT_DIRTY,
            host:           built_info::HOST,
            opt_level:      built_info::OPT_LEVEL,
            profile:        built_info::PROFILE,
            rustc:          built_info::RUSTC_VERSION,
            target:         built_info::TARGET,
        }
    }
}


impl<'a> fmt::Display for Version<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.logo)?;
        writeln!(f, "Welcome to use {}!", self.name)?;
        writeln!(f, "    current version:    {}", self.version_str)?;
        writeln!(f, "    git hash:           {}", self.git_hash_short.unwrap_or("NO GIT INFO"))?;
        writeln!(f, "    authors:            {}", self.authors)?;
        writeln!(f, "    host:               {}", self.host)?;
        writeln!(f, "    build time:         {}", self.built_time)?;

        if f.alternate() {
            writeln!(f, "        git_hash_long:  {}", self.git_hash_long.unwrap_or("NO GIT INFO"))?;
            if self.git_hash_long.is_some() {
                writeln!(f, "        is git dirty?   {:#}", self.git_dirty.unwrap())?;
            }

            writeln!(f, "        opt level:      {}", self.opt_level)?;
            writeln!(f, "        build profile:  {}", self.profile)?;
            writeln!(f, "        rustc version:  {}", self.rustc)?;
            writeln!(f, "        build target:   {}", self.target)?;
        }
        Ok(())
    }
}
