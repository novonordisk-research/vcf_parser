// Read from stdin or a .vcf[.gz] file
use calm_io::stdoutln;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str;
use vcf::{VCFReader, VCFRecord};
use crate::variant::Variant;

pub mod variant;

#[derive(Parser)]
#[command(version = "0.1.0", about = "Read a .vcf[.gz] and produce json objects per line. CSQ-aware, multiallelic-aware", long_about = None)]
#[command(styles=get_styles())]
struct Args {
    /// input .vcf[.gz] file, or ignore to read from stdin
    #[arg(short, long, value_parser = vcf_extension_validator, default_value = "-")]
    input: String,

    /// threads to use, default to use all available
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// filter yaml file to use
    #[arg(short, long, default_value = "-")]
    filter: String,
}





fn parse_csq_header(header: &str) -> Vec<String> {
    // parse the csq header to produce a list of fields
    header.split(": ").collect::<Vec<&str>>()[1]
        .split("|")
        .map(|x| x.trim().to_string())
        .collect::<Vec<String>>()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    // if thread is not 0, set the number of threads
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }
    let reader: Box<dyn BufRead + Send + Sync> = if args.input == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else if args.input.ends_with(".vcf.gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(File::open(args.input)?)))
    } else {
        Box::new(BufReader::new(File::open(args.input)?))
    };

    let reader = VCFReader::new(reader)?;
    // get info and vep headers
    let mut info_headers: Vec<&str> = Vec::new();
    let mut csq_headers: HashMap<String, Vec<String>> = HashMap::new();
    for info in reader.header().info_list() {
        let info_str = str::from_utf8(&info)?;
        info_headers.push(&info_str);
        // if in the description it has 'Format: Allele|' then presume it is a csq header
        let desc = str::from_utf8(reader.header().info(info).unwrap().description).unwrap();
        if desc.contains(": Allele") {
            csq_headers.insert(info_str.to_string(), parse_csq_header(desc));
        }
    }
    let header = reader.header().clone();
    reader.reader.lines().par_bridge().for_each(|line| {
        let line = line.unwrap();
        if line.starts_with("#") {
            return;
        }
        let line = line.as_bytes() as &[u8];
        let vcf_record = VCFRecord::from_bytes(line, 1, header.clone()).unwrap();
        let variant = Variant::new(&vcf_record, header.samples(), &csq_headers);
        let j = serde_json::to_string(&variant).unwrap();
        match stdoutln!("{}", j) {
            Ok(_) => Ok(()),
            Err(e) => match e.kind() {
                std::io::ErrorKind::BrokenPipe => std::process::exit(0),
                _ => Err(e),
            },
        }
        .unwrap();
    });
    Ok(())
}

fn vcf_extension_validator(fname: &str) -> Result<String, String> {
    if fname == "-" {
        return Ok(format!("{fname}"));
    }
    if !fname.contains(".") {
        return Err(format!(
            "Input file has to be a .vcf or a .gz file, or use - to read from stdin"
        ));
    }
    let file_extension = Path::new(&fname).extension().unwrap().to_str().unwrap();
    if ["vcf", "gz"].contains(&file_extension) {
        // if return a &str from the input it won't compile
        if !Path::new(fname).exists() {
            return Err(format!("input file {fname} does not exist"));
        }
        Ok(format!("{fname}"))
    } else {
        Err(format!(
            "Input file has to be a .vcf or a .gz file, or use - to read from stdin"
        ))
    }
}

fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}
