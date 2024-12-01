// Read from stdin or a .vcf[.gz] file

use flate2::read::MultiGzDecoder;
use clap::Parser;
use rayon::prelude::*;
use vcfparser::VcfParser;
use std::{fs::File, error::Error, str, path::Path};
use std::io::{self, BufRead, BufReader};
pub use vcf::VCFRecord;
use crate::variant::Variant;
use serde_json::{Map, Value};
use serde::Serialize;
use anyhow::Result;

pub mod vcfparser;
pub mod variant;
pub mod utils;
pub mod error;
pub mod parser;

#[derive(Parser)]
#[command(version = "0.2.4", about = "Read a (normalised) .vcf[.gz] and output tsv/json. CSQ-aware.", long_about = None)]
#[command(styles=get_styles())]
pub struct Args {
    /// input .vcf[.gz] file, or ignore to read from stdin
    #[arg(short, long, value_parser = vcf_extension_validator)]
    input: Option<String>,

    /// threads to use, default to use all available
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// filter expression (experimental) or yaml file to use
    #[arg(short, long)]
    filter: Option<String>,

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



pub fn run(args:Args)-> Result<(), Box<dyn Error>> {
    // read filter if given
    // if space is present, treat it as a logic expression
    // otherwise, treat it as a file
    let filters: serde_json::Value = if let Some(filter_string) = args.filter  {
        if filter_string.contains(" ") {
            parser::parse_logic_expr(&filter_string).map_err(|e| error::VcfParserError::InvalidFilter(e.to_string()))?
        } else {
            let filter_file = File::open(filter_string)?;
            serde_yaml::from_reader(filter_file)?
        }
    } else {
        serde_json::Value::Null
    };
    
    let reader: Box<dyn BufRead + Send + Sync> = match args.input {
        None => Box::new(BufReader::new(io::stdin())),
        Some(inp) if inp.ends_with(".vcf.gz") => Box::new(BufReader::new(MultiGzDecoder::new(File::open(inp)?))),
        Some(rest) => Box::new(BufReader::new(File::open(rest)?))
    };
    let vcf_parser = VcfParser::new(filters, args.fields, args.fields_join, args.columns, args.output_format, reader)?;

    // if --list, print the headers and quit
    if args.list {
        let header = vcf_parser.tsv_headers;
        utils::print_line_to_stdout(&header.join("\n"))?;
        return Ok(());
    }
    // write tsv header to stdout
    let tsv_header = vcf_parser.tsv_headers;
    if vcf_parser.output_format == OutputFormat::T {
        
        utils::print_line_to_stdout(&tsv_header.join("\t"))?;
    }
    
    // parallel processing each variant/site
    
    vcf_parser.reader.reader.lines().par_bridge().for_each(|line| {
        let line = line.unwrap();
        if line.starts_with("#") {
            return;
        }
        let line = line.as_bytes() as &[u8];
        let vcf_record = VCFRecord::from_bytes(line, 1, (*vcf_parser.header).clone()).unwrap();
        let variant = Variant::new(&vcf_record, vcf_parser.header.samples(), &vcf_parser.csq_headers);
        let explodeds = vcf_parser.info_fields.iter().map(|x| utils::explode_data(serde_json::to_value(&variant).unwrap(), x, &vcf_parser.info_fields)).collect::<Vec<Vec<Map<String, Value>>>>();
        let joined = utils::outer_join(explodeds, &vcf_parser.fields_join).unwrap();
        joined.iter().filter(|x| utils::filter_record(x, &vcf_parser.filters)).for_each(|x| {
            match vcf_parser.output_format {
                OutputFormat::T => {
                    let row = utils::get_row(&x, &tsv_header);
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
pub enum OutputFormat {
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
    use std::collections::HashMap;
    use vcf::VCFReader;
    pub fn prepare_test(vcf_file: Option<&str>, fields:&Vec<String>)-> Result<(VCFReader<Box<dyn BufRead + Send + Sync>>, HashMap<String, Vec<String>>, Value), Box<dyn Error>> {
        let filter_file = File::open("test/filter.yml")?;
        let filter = serde_yaml::from_reader(filter_file)?;
        let vcf_file = vcf_file.unwrap_or("test/test.vcf");
        let reader: Box<dyn BufRead + Send + Sync> = 
            Box::new(BufReader::new(File::open(vcf_file)?));
    
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
        let (mut reader, csq_headers, _filter) = prepare_test(None, &vec!["CSQ".to_string(), "Pangolin".to_string()])?;
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
        let (reader, csq_headers, _filter) = prepare_test(None, &vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut info_headers: Vec<String> = Vec::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(info_str.to_string());
        }
        let header = utils::get_output_header(&info_headers, &csq_headers, reader.header().samples(), &None);
        let expected = vec!["chromosome", "position", "id", "reference", "alternative", "qual", "filter", "info.AC", "info.AF", "info.CADD_PHRED", "info.CADD_RAW", "info.CSQ.Allele", "info.CSQ.CANONICAL", "info.CSQ.Consequence", "info.CSQ.Feature", "info.CSQ.Feature_type", "info.CSQ.Gene", "info.CSQ.IMPACT", "info.CSQ.SYMBOL", "info.Pangolin.pangolin_gene", "info.Pangolin.pangolin_max_score", "info.Pangolin.pangolin_transcript", "info.tag", "info.what", "info.who"];
        assert_eq!(header, expected);

        let header = utils::get_output_header(&info_headers, &csq_headers, reader.header().samples(), &Some(vec!["info.CSQ.Consequence".to_string(), "reference".to_string()]));
        let expected = vec!["info.CSQ.Consequence", "reference"];
        assert_eq!(header, expected);

        Ok(())
    }
    #[test]
    fn test_get_output_header_with_samples() -> Result<(), Box<dyn Error>> {
        let (reader, csq_headers, _filter) = prepare_test(Some("test/test_samples.vcf"), &vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut info_headers: Vec<String> = Vec::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(info_str.to_string());
        }
        let header = utils::get_output_header(&info_headers, &csq_headers, reader.header().samples(), &None);
        let expected = vec!["chromosome", "position", "id", "reference", "alternative", "qual", "filter", "info.AC", "info.AF", "info.CADD_PHRED", "info.CADD_RAW", "info.CSQ.Allele", "info.CSQ.CANONICAL", "info.CSQ.Consequence", "info.CSQ.Feature", "info.CSQ.Feature_type", "info.CSQ.Gene", "info.CSQ.IMPACT", "info.CSQ.SYMBOL", "info.Pangolin.pangolin_gene", "info.Pangolin.pangolin_max_score", "info.Pangolin.pangolin_transcript", "info.tag", "info.what", "info.who", "S1", "S2"];
        assert_eq!(header, expected);

        let header = utils::get_output_header(&info_headers, &csq_headers, reader.header().samples(), &Some(vec!["info.CSQ.Consequence".to_string(), "reference".to_string()]));
        let expected = vec!["info.CSQ.Consequence", "reference"];
        assert_eq!(header, expected);

        Ok(())
    }
    #[test]
    #[should_panic]
    fn test_get_output_header_panic() {
        let (reader, csq_headers, _filter) = prepare_test(None, &vec!["CSQ".to_string(), "Pangolin".to_string()]).unwrap();
        let mut info_headers: Vec<String> = Vec::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info).unwrap();
            info_headers.push(info_str.to_string());
        }
        utils::get_output_header(&info_headers, &csq_headers, reader.header().samples(), &Some(vec!["info.CSQ.Consequence".to_string(), "doesnotexist".to_string()]));
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
        let (mut reader, csq_headers, filter) = prepare_test(None, &vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut vcf_record = reader.empty_record();
        let fields = vec!["info.CSQ".to_string(), "info.Pangolin".to_string()];
        let fields_join = vec!["info.CSQ.Feature".to_string(), "info.Pangolin.pangolin_transcript".to_string()];
        let mut results: Vec<Map<String,Value>> = Vec::new();
        while reader.next_record(&mut vcf_record).unwrap() {
            let variant = Variant::new(&vcf_record, reader.header().samples(), &csq_headers);
            let explodeds = fields.iter().map(|x| utils::explode_data(serde_json::to_value(&variant).unwrap(), x, &fields)).collect::<Vec<Vec<Map<String, Value>>>>();
            let joined = utils::outer_join(explodeds, &fields_join)?;
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

    #[test]
    fn test_logic() -> Result<(), Box<dyn Error>> {

        let exprs = [
            (r#"foo = bar AND baz > 10"#, r#"{"AND":[{"name":"foo","op":"=","value":"bar"},{"name":"baz","op":">","value":10.0}]}"#),
            (r#"(foo == bar or baz > 10) AND baz <= 5"#, r#"{"AND":[{"OR":[{"name":"foo","op":"==","value":"bar"},{"name":"baz","op":">","value":10.0}]},{"name":"baz","op":"<=","value":5.0}]}"#),
        ];
        for (expr, expected) in exprs.iter() {
            let result = parser::parse_logic_expr(expr)?;
            assert_eq!(result, serde_json::from_str::<serde_json::Value>(expected)?);
        }
        Ok(())
    }
}