use thiserror::Error;

#[derive(Error, Debug)]
pub enum VcfParserError<'a> {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(#[from] vcf::VCFError),
    #[error("Invalid argument: {0}")]
    InvalidArgument(String),
    #[error("Invalid filter expression: {0}")]
    InvalidFilter(String),
    #[error("variant parse error: {0}")]
    VariantParseError(&'a str),
}