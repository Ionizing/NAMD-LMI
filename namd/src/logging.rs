use std::path::Path;
use once_cell::sync::Lazy;
use std::sync::Mutex;

use shared::{
    Result,
    log::LevelFilter,
};
use log4rs::{
    append::{
        console::{
            ConsoleAppender,
            Target,
        },
        file::FileAppender,
    },
    config::{
        Appender,
        Config,
        Root,
    },
    encode::pattern::PatternEncoder,
    init_config,
    Handle,
};


pub static HANDLE: Lazy<Mutex<Handle>> = Lazy::new(|| {
    let config = gen_logger_config(Option::<&str>::None).unwrap();
    let handle = init_config(config).unwrap();
    Mutex::new(handle)
});


fn gen_logger_config(path: Option<impl AsRef<Path>>) -> Result<Config> {
    const ENCODE_STR: &'static str = "{d(%Y-%m-%d %H:%M:%S)} [{h({l:>5})}] {m}{n}";

    let level = LevelFilter::Info;

    let stderr = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new(ENCODE_STR)))
        .target(Target::Stderr)
        .build();

    let global_log = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new(ENCODE_STR)))
        .build("./globalrun.log")
        .unwrap();

    let logfile = path.map(|p| p.as_ref().join("run.log"))
        .map(|file_path| {
            FileAppender::builder()
                .encoder(Box::new(PatternEncoder::new(ENCODE_STR)))
                .build(file_path)
                .unwrap()
        });

    let config = if let Some(logfile) = logfile {
        Config::builder()
            .appender(Appender::builder().build("stderr", Box::new(stderr)))
            .appender(Appender::builder().build("global_log", Box::new(global_log)))
            .appender(Appender::builder().build("logfile", Box::new(logfile)))
            .build(
                Root::builder()
                    .appender("stderr")
                    .appender("global_log")
                    .appender("logfile")
                    .build(level)
            )?
    } else {
        Config::builder()
            .appender(Appender::builder().build("stderr", Box::new(stderr)))
            .appender(Appender::builder().build("global_log", Box::new(global_log)))
            .build(
                Root::builder()
                    .appender("stderr")
                    .appender("global_log")
                    .build(level)
            )?
    };

    return Ok(config);
}


pub fn logger_init() {
    Lazy::force(&HANDLE);
}


pub fn logger_redirect(path: impl AsRef<Path>) -> Result<()> {
    let config =gen_logger_config(Some(path))?;
    HANDLE.lock().unwrap().set_config(config);
    Ok(())
}
