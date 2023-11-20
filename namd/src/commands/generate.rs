use clap::Parser;


#[derive(Parser)]
pub struct Generate {
    #[arg(default_value_t = String::from("input-example.toml") )]
    /// Where to write the example input file.
    filename: String,
}
