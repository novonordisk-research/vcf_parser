// Read from stdin or a .vcf[.gz] file

use flate2::read::MultiGzDecoder;
use clap::Parser;
use rayon::prelude::*;
use std::{fs::File, collections::HashMap, error::Error, str, path::Path};
use std::io::{self, BufRead, BufReader};
use vcf::{VCFReader, VCFRecord};
use crate::variant::Variant;
use serde_json::{Map, Value};
use serde::Serialize;

mod variant;
mod utils;

#[derive(Parser)]
#[command(version = "0.2.2", about = "Read a (normalised) .vcf[.gz] and output tsv/json. CSQ-aware.", long_about = None)]
#[command(styles=get_styles())]
pub struct Args {
    /// input .vcf[.gz] file, or ignore to read from stdin
    #[arg(short, long, value_parser = vcf_extension_validator, default_value = "-")]
    input: String,

    /// threads to use, default to use all available
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// filter yaml file to use
    #[arg(short, long, default_value = "-")]
    filter: String,

    /// list of columns available for query/output
    #[arg(short, long, default_value_t = false)]
    list: bool,

    /// nested fields with | separator to parse, such as CSQ
    #[arg(long, value_parser, value_delimiter = ',', default_value = "CSQ")]
    fields: Vec<String>,

    /// nested fields join together on e.g. transcript id, such as Feature
    #[arg(long, value_parser, value_delimiter = ',', default_value = "Feature")]
    fields_join: Vec<String>,

    /// specify output columns
    #[arg(short, long, value_delimiter = ',')]
    columns: Option<Vec<String>>,

    /// output format.
    #[arg(long, default_value_t, value_enum)]
    output_format: OutputFormat,
}

pub fn run(args:Args) -> Result<(), Box<dyn Error>> {
    // if args.fields is greater than 1, then its length should be the same with args.fields_join
    // if args.fields is 1, then args.fields_join
    if args.fields.len() >1 && args.fields.len() != args.fields_join.len() {
        return Err("Number of fields should be equal to the number of fields_join".into());
    }
    // if thread is not 0, set the number of threads
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }
    // if thread is not 0, set the number of threads
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }
    // read filter if given
    let filters: serde_json::Value = if args.filter == "-" {
        serde_json::Value::Null
    } else {
        let filter_file = File::open(args.filter)?;
        serde_yaml::from_reader(filter_file)?
    };

    // read input.
    // if "-", read from stdin.
    // if *.vcf.gz, use gzdecoder
    // otherwise just plain text read
    let reader: Box<dyn BufRead + Send + Sync> = match args.input.as_str() {
        "-" => Box::new(BufReader::new(io::stdin())),
        inp if inp.ends_with(".vcf.gz") => Box::new(BufReader::new(MultiGzDecoder::new(File::open(args.input)?))),
        _ => Box::new(BufReader::new(File::open(args.input)?))
    };

    /* 
    vcfrs reads the vcf header from the stream, then reads a variant/site at a time to couple with the header to make a VCFRecord.
    To enable multithread, I would need to make a VCFRecord in each thread when it reads a line, and couple it with the header.
    This would require the access of a private method (VCFRecord::reader) in the vcfrs repo, which isn't possible. Therefore I forked that repo and exposed the field..
    */

    let reader = VCFReader::new(reader)?;
    // get info and csq-like headers
    let mut info_headers: Vec<&str> = Vec::new();
    let mut csq_headers: HashMap<String, Vec<String>> = HashMap::new();
    for info in reader.header().info_list() {
        let info_str = str::from_utf8(&info)?;
        info_headers.push(&info_str);
        let desc = str::from_utf8(reader.header().info(info).unwrap().description).unwrap();
        if args.fields.contains(&info_str.to_string()) {
            csq_headers.insert(info_str.to_string(), utils::parse_csq_header(desc));
        }
    }
    let header = reader.header().clone();
    

    // info_fields are fields with 'info.' prefix, so CSQ becomes info.CSQ
    let info_fields = args.fields.iter().map(|x| format!("info.{}", x)).collect::<Vec<String>>();
    // Feature becomes info.CSQ.Feature
    let info_fields_join = args.fields_join.iter().enumerate().map(|(ind, x)| format!("{}.{}", info_fields[ind], x)).collect::<Vec<String>>();

    // if tsv, write the header
    let tsv_header = utils::get_output_header(info_headers, &csq_headers, &args.columns);

    // flush header and quit if --l
    if args.list {
        utils::print_line_to_stdout(&tsv_header.join("\n"))?;
        return Ok(());
    }

    // write tsv header to stdout
    if args.output_format == OutputFormat::T {
        
        utils::print_line_to_stdout(&tsv_header.join("\t"))?;
    }
    
    // parallel processing each variant/site
    reader.reader.lines().par_bridge().for_each(|line| {
        let line = line.unwrap();
        if line.starts_with("#") {
            return;
        }
        let line = line.as_bytes() as &[u8];
        let vcf_record = VCFRecord::from_bytes(line, 1, header.clone()).unwrap();
        let variant = Variant::new(&vcf_record, header.samples(), &csq_headers);
        let explodeds = info_fields.iter().map(|x| utils::explode_data(serde_json::to_value(&variant).unwrap(), x, &info_fields)).collect::<Vec<Vec<Map<String, Value>>>>();
        let joined = utils::outer_join(explodeds, &info_fields_join).unwrap();
        joined.iter().filter(|x| utils::filter_record(x, &filters)).for_each(|x| {
            match args.output_format {
                OutputFormat::T => {
                    let row = utils::get_row(x.clone(), &tsv_header);
                    utils::print_line_to_stdout(&row.join("\t")).unwrap();
                },
                OutputFormat::J => {
                    let j = serde_json::to_string(&x).unwrap();
                    utils::print_line_to_stdout(&j).unwrap();
                },
                OutputFormat::V => {
                    unimplemented!();
                }
            }
        });
        
        
    });
    
    Ok(())
}



#[derive(
    clap::ValueEnum, Clone, Default, Debug, Serialize,
)]
#[derive(PartialEq)]
#[serde(rename_all = "kebab-case")]
enum OutputFormat {
    /// tsv
    #[default]
    T,
    /// json
    J,
    /// VCF
    V,
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





#[cfg(test)]
mod tests {
    use super::*;
    pub fn prepare_test(fields:&Vec<String>)-> Result<(VCFReader<Box<dyn BufRead + Send + Sync>>, HashMap<String, Vec<String>>, Value), Box<dyn Error>> {
        let filter_file = File::open("test/filter.yml")?;
        let filter = serde_yaml::from_reader(filter_file)?;
        let reader: Box<dyn BufRead + Send + Sync> = 
            Box::new(BufReader::new(File::open("test/test.vcf")?));
    
        let reader = VCFReader::new(reader)?;
        // get info and vep headers
        let mut info_headers: Vec<&str> = Vec::new();
        let mut csq_headers: HashMap<String, Vec<String>> = HashMap::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(&info_str);
            // if in the description it has 'Format: Allele|' then presume it is a csq header
            let desc = str::from_utf8(reader.header().info(info).unwrap().description).unwrap();
            //if desc.contains(": Allele") {
            if fields.contains(&info_str.to_string()) {
                csq_headers.insert(info_str.to_string(), utils::parse_csq_header(desc));
            }
        }
        Ok((reader, csq_headers, filter))
    }

    #[test]
    fn test_variants() -> Result<(), Box<dyn Error>> {
        let (mut reader, csq_headers, _filter) = prepare_test(&vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut vcf_record = reader.empty_record();
        let mut variants: Vec<Variant> = Vec::new();
        
        while reader.next_record(&mut vcf_record).unwrap() {
            let variant = Variant::new(&vcf_record, reader.header().samples(), &csq_headers);
            variants.push(variant);
            
        }
        assert!(variants.len() == 5);
        
        Ok(())
    }

    #[test]
    fn test_get_output_header() -> Result<(), Box<dyn Error>> {
        let (reader, csq_headers, _filter) = prepare_test(&vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut info_headers: Vec<&str> = Vec::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(&info_str);
        }
        let header = utils::get_output_header(info_headers, &csq_headers);
        let expected = vec!["chromosome", "position", "id", "reference", "alternative", "qual", "filter", "info.AC", "info.AF", "info.CADD_PHRED", "info.CADD_RAW", "info.CSQ.Allele", "info.CSQ.CANONICAL", "info.CSQ.Consequence", "info.CSQ.Feature", "info.CSQ.Feature_type", "info.CSQ.Gene", "info.CSQ.IMPACT", "info.CSQ.SYMBOL", "info.Pangolin.pangolin_gene", "info.Pangolin.pangolin_max_score", "info.Pangolin.pangolin_transcript", "info.tag", "info.what", "info.who"];
        assert_eq!(header, expected);
        Ok(())
    }
    #[test]
    fn test_explode_data() -> Result<(), Box<dyn Error>> {
        let data = serde_json::from_str(r#"{"a":1, "b":2, "c":[{"foo":"A","bar":"B"},{"foo":"a", "bar":"b"}]}"#)?;
        let exploded = utils::explode_data(data, "c", &vec!["b".to_string(),"c".to_string()]);
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"a":1, "c.foo":"A", "c.bar":"B"}"#).unwrap(),
            serde_json::from_str(r#"{"a":1, "c.foo":"a", "c.bar":"b"}"#).unwrap(),
        ];
        assert_eq!(exploded, expected);
        Ok(())
    }

    #[test]
    fn test_filter() -> Result<(), Box<dyn Error>> {
        let (mut reader, csq_headers, filter) = prepare_test(&vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut vcf_record = reader.empty_record();
        let fields = vec!["info.CSQ".to_string(), "info.Pangolin".to_string()];
        let fields_join = vec!["info.CSQ.Feature".to_string(), "info.Pangolin.pangolin_transcript".to_string()];
        let mut results: Vec<Map<String,Value>> = Vec::new();
        while reader.next_record(&mut vcf_record).unwrap() {
            let variant = Variant::new(&vcf_record, reader.header().samples(), &csq_headers);
            //let val = serde_json::to_value(&variant)?;
            //println!("{:?}", serde_json::to_string(&variant)?);
            let explodeds = fields.iter().map(|x| utils::explode_data(serde_json::to_value(&variant).unwrap(), x, &fields)).collect::<Vec<Vec<Map<String, Value>>>>();
            let joined = utils::outer_join(explodeds, &fields_join)?;
            //println!("====================");
            //println!("{}", serde_json::to_string_pretty(&joined)?);
            let filtered_record: Vec<Map<String, Value>> = joined.into_iter().filter(|x| utils::filter_record(x, &filter)).collect();
            results.extend(filtered_record.into_iter());
        }
        assert_eq!(results.len(), 6);
        Ok(())
    }

    #[test]
    fn test_outer_join1() -> Result<(), Box<dyn Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key3": "value6"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key3": "value8"}"#).unwrap(),
        ];
        let joined_table = utils::outer_join(vec![table1, table2], &vec!["key1".to_string(), "key1".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2", "key3": "value6"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key2": null, "key3": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4", "key3": null}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
    #[test]
    fn test_outer_join2() -> Result<(), Box<dyn Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key3": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key3": "value8"}"#).unwrap(),
        ];
        let joined_table = utils::outer_join(vec![table1, table2], &vec!["key1".to_string(), "key1".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2", "key3": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key2": null, "key3": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4", "key3": null}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
    #[test]
    fn test_outer_join3() -> Result<(), Box<dyn Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key11": "value1.8", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key11": "value7", "key2": "value8"}"#).unwrap(),
        ];
        let joined_table = utils::outer_join(vec![table1, table2], &vec!["key1".to_string(), "key11".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key11": "value1.8", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": null, "key11": "value7", "key2": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key11": null, "key2": "value4"}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
}