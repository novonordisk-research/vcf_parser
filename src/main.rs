// Read from stdin or a .vcf[.gz] file
use calm_io::stdoutln;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Number, Value};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str;
use vcf::{VCFReader, VCFRecord};

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
}

#[derive(Serialize, Deserialize)]
struct Variant {
    chromosome: String,
    position: u64,
    id: Vec<String>,
    reference: String,
    alternative: Vec<String>,
    qual: Option<f64>,
    filter: Vec<String>,
    info: Map<String, Value>,
    genotype: Map<String, Value>,
}

fn try_parse_number(input: &str) -> Value {
    if input.contains('.') || input.contains('e') || input.contains('E') {
        // Try to parse as f64
        match input.parse::<f64>() {
            Ok(num) => Value::Number(Number::from_f64(num).unwrap()),
            Err(_) => Value::String(input.to_string()),
        }
    } else if input == "" {
        Value::Null
    } else {
        // Try to parse as i64
        match input.parse::<i64>() {
            Ok(num) => Value::Number(Number::from(num)),
            Err(_) => Value::String(input.to_string()),
        }
    }
}

impl Variant {
    fn new(
        vcf_record: &VCFRecord,
        samples: &[Vec<u8>],
        csq_headers: &HashMap<String, Vec<String>>,
    ) -> Self {
        // parse genotype
        let mut genotype = Map::new();
        for sample in samples {
            let mut sample_genotype = Map::new();
            for key in &vcf_record.format {
                let genotype_raw = vcf_record.genotype(&sample, &key);
                let f_k = str::from_utf8(key).unwrap();
                let val = str::from_utf8(&genotype_raw.unwrap()[0]).unwrap();
                match vcf_record.header().format(key).unwrap().value_type {
                    vcf::ValueType::Integer => sample_genotype.insert(
                        f_k.to_string(),
                        Value::Number(Number::from(val.parse::<i64>().unwrap())),
                    ),
                    vcf::ValueType::Float => sample_genotype.insert(
                        f_k.to_string(),
                        Value::Number(Number::from_f64(val.parse::<f64>().unwrap()).unwrap()),
                    ),
                    _ => sample_genotype.insert(f_k.to_string(), Value::String(val.to_string())),
                };
            }
            genotype.insert(
                str::from_utf8(&sample).unwrap().to_string(),
                Value::Object(sample_genotype),
            );
        }
        // parse info
        let mut info = Map::new();
        for id in vcf_record.header().info_list() {
            let field = vcf_record.header().info(id).unwrap();
            let field_str = str::from_utf8(&field.id).unwrap();
            let val = match vcf_record.info(id) {
                Some(dat) => {
                    if *field.value_type == vcf::ValueType::Flag {
                        // flag type
                        Value::Bool(true)
                    } else if *field.value_type == vcf::ValueType::Integer {
                        if *field.number == vcf::Number::Allele {
                            Value::Array(
                                dat.iter()
                                    .map(|x| {
                                        Value::from(Number::from(
                                            str::from_utf8(x).unwrap().parse::<i64>().unwrap(),
                                        ))
                                    })
                                    .collect::<Vec<Value>>(),
                            )
                        } else {
                            Value::from(str::from_utf8(&dat[0]).unwrap().parse::<i64>().unwrap())
                        }
                    } else if *field.value_type == vcf::ValueType::Float {
                        if *field.number == vcf::Number::Allele {
                            Value::Array(
                                dat.iter()
                                    .map(|x| {
                                        Value::from(Number::from_f64(
                                            str::from_utf8(x).unwrap().parse::<f64>().unwrap(),
                                        ))
                                    })
                                    .collect::<Vec<Value>>(),
                            )
                        } else {
                            Value::from(str::from_utf8(&dat[0]).unwrap().parse::<f64>().unwrap())
                        }
                    } else if csq_headers.contains_key(field_str) {
                        dat.iter()
                            .map(|csq_field| {
                                let csq_vec = csq_headers
                                    .get(field_str)
                                    .unwrap()
                                    .iter()
                                    .zip(str::from_utf8(csq_field).unwrap().split("|"));
                                let mut csq = Map::new();
                                for (k, v) in csq_vec {
                                    csq.insert(k.to_string(), try_parse_number(v));
                                }
                                Value::Object(csq)
                            })
                            .collect::<Value>()
                    } else {
                        if *field.number == vcf::Number::Allele {
                            Value::Array(
                                dat.iter()
                                    .map(|x| Value::String(str::from_utf8(x).unwrap().to_string()))
                                    .collect::<Vec<Value>>(),
                            )
                        } else {
                            Value::String(str::from_utf8(&dat[0]).unwrap().to_string())
                        }
                    }
                }
                None => {
                    if *field.value_type == vcf::ValueType::Flag {
                        // flag type
                        Value::Bool(false)
                    } else {
                        continue;
                    }
                }
            };
            info.insert(str::from_utf8(&field.id).unwrap().to_string(), val);
        }
        Variant {
            chromosome: str::from_utf8(&vcf_record.chromosome).unwrap().to_string(),
            position: vcf_record.position,
            id: vcf_record
                .id
                .iter()
                .map(|x| str::from_utf8(&x).unwrap().to_string())
                .collect::<Vec<String>>(),
            reference: str::from_utf8(&vcf_record.reference).unwrap().to_string(),
            alternative: vcf_record
                .alternative
                .iter()
                .map(|x| str::from_utf8(&x).unwrap().to_string())
                .collect::<Vec<String>>(),
            qual: vcf_record.qual,
            filter: vcf_record
                .filter
                .iter()
                .map(|x| str::from_utf8(&x).unwrap().to_string())
                .collect::<Vec<String>>(),
            info,
            genotype,
        }
    }
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
