use thiserror::Error;

#[derive(Error, Debug)]
pub enum VcfParserError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(#[from] vcf::VCFError),
    #[error("Invalid argument: {0}")]
    InvalidArgument(String),
}