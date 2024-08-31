// Read from stdin or a .vcf[.gz] file
use calm_io::stdoutln;
use flate2::read::MultiGzDecoder;
use clap::Parser;
use rayon::prelude::*;
use core::panic;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str;
use vcf::{VCFReader, VCFRecord};
use crate::variant::Variant;
use serde_json::{Map, Value};

mod variant;
#[derive(Parser)]
#[command(version = "0.1.1", about = "Read a .vcf[.gz] and produce json objects per line. CSQ-aware, multiallelic-aware", long_about = None)]
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

    /// nested fields with | separator to parse, such as CSQ
    #[arg(long, value_parser, value_delimiter = ',', default_value = "CSQ")]
    fields: Vec<String>,

    /// nested fields join together on e.g. transcript id, such as Feature
    #[arg(long, value_parser, value_delimiter = ',', default_value = "Feature")]
    fields_join: Vec<String>,
}
fn parse_csq_header(header: &str) -> Vec<String> {
    // parse the csq header to produce a list of fields
    header.split(": ").collect::<Vec<&str>>()[1]
        .split("|")
        .map(|x| x.trim().to_string())
        .collect::<Vec<String>>()
}

fn filter_record(record: &Map<String,Value>, filters: &Value) -> bool {
    // filter the variant based on the filters. Filter is like:
    /*
    {  
        "AND":[
            {
                "AND":[
                    {"name":"gnomAD_exome_V4.0_AF","op":"le","value":0.05},
                    {"name":"gnomAD_genome_V4.0_AF","op":"le","value":0.05}
                ]
            },
            {
                "OR":[
                    {"name":"CADD_PHRED","op":"ge","value":20},
                    {"name":"Lof","op":"ne","value":null}
                ]
            }
        ]
    }
    This function is recursive, so it can handle nested AND/OR
    */
    match filters {
        Value::Object(map) => {
            for (k, v) in map {
                if k == "AND" {
                    match v {
                        Value::Array(vv) => {
                            if vv.into_iter().all(|x| filter_record(record, x)) {
                                return true;
                            } else {
                                return false;
                            }
                        }
                        _ => panic!("AND should be a list of filters"),
                    }
                } else if k == "OR" {
                    match v {
                        Value::Array(vv) => {
                            if vv.into_iter().any(|x| filter_record(record, x)) {
                                return true;
                            } else {
                                return false;
                            }
                        }
                        _ => panic!("OR should be a list of filters"),
                    }
                } else {
                    let name = map["name"].as_str().unwrap();
                    let op = map["op"].as_str().unwrap();
                    let value = &map["value"];
                    let val = record.get(name).unwrap_or(&Value::Null);
                    match op {
                        "eq" => {
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.contains(value) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if *val == *value {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        "ne" => {
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.iter().any(|x| x != value) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if *val != *value {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        "gt" => {
                            // if val is null, returning false. essentially treat null as 0, and assume value is positive.
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.iter().any(|x| {
                                        if x.is_null() {
                                            if 0. > value.as_f64().unwrap() {
                                                return true;
                                            } else {
                                                return false;
                                            }
                                        }
                                        x.as_f64().unwrap() > value.as_f64().unwrap()
                                    }) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if val.is_null() {
                                        if 0. > value.as_f64().unwrap() {
                                            return true;
                                        } else {
                                            return false;
                                        }
                                    }
                                    if val.as_f64().unwrap() > value.as_f64().unwrap() {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        "ge" => {
                            // if val is null, treat it as 0.
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.iter().any(|x| {
                                        if x.is_null() {
                                            if 0. >= value.as_f64().unwrap() {
                                                return true;
                                            } else {
                                                return false;
                                            }
                                        }
                                        x.as_f64().unwrap() >= value.as_f64().unwrap()
                                    }) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if val.is_null() {
                                        if 0. >= value.as_f64().unwrap() {
                                            return true;
                                        } else {
                                            return false;
                                        }
                                    }
                                    if val.as_f64().unwrap() >= value.as_f64().unwrap() {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        "lt" => {
                            // if val is null, returning true. essentially treat null as 0, and assume value is positive.
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.iter().any(|x| {
                                        if x.is_null() {
                                            if 0. < value.as_f64().unwrap() {
                                                return true;
                                            } else {
                                                return false;
                                            }
                                        }
                                        x.as_f64().unwrap() < value.as_f64().unwrap()
                                    }) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if val.is_null() {
                                        if 0. < value.as_f64().unwrap() {
                                            return true;
                                        } else {
                                            return false;
                                        }
                                    }
                                    if val.as_f64().unwrap() < value.as_f64().unwrap() {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        "le" => {
                            // if val is null, returning true. essentially treat null as 0, and assume value is positive.
                            match val {
                                serde_json::Value::Array(arr) => {
                                    if arr.iter().any(|x| {
                                        if x.is_null() {
                                            if 0. <= value.as_f64().unwrap() {
                                                return true;
                                            } else {
                                                return false;
                                            }
                                        }
                                        x.as_f64().unwrap() <= value.as_f64().unwrap()
                                    }) {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                                _ => {
                                    if val.is_null() {
                                        if 0. <= value.as_f64().unwrap() {
                                            return true;
                                        } else {
                                            return false;
                                        }
                                    }
                                    if val.as_f64().unwrap() <= value.as_f64().unwrap() {
                                        return true;
                                    } else {
                                        return false;
                                    }
                                }
                            }
                        }
                        _ => panic!("Unknown operator"),
                    }
                }
            };
            false
        },
        serde_json::Value::Null => true,
        _ => false
    }
}
pub fn run(args:Args) -> Result<(), Box<dyn std::error::Error>> {
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
    // read filter if given
    let filters: serde_json::Value = if args.filter == "-" {
        serde_json::Value::Null
    } else {
        let filter_file = File::open(args.filter)?;
        serde_yaml::from_reader(filter_file)?
    };
    // read input
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
        //if desc.contains(": Allele") {
        if args.fields.contains(&info_str.to_string()) {
            csq_headers.insert(info_str.to_string(), parse_csq_header(desc));
        }
    }
    let header = reader.header().clone();
    //let csv_header = 

    // info_fields are fields with 'info.' prefix, so CSQ becomes info.CSQ
    let info_fields = args.fields.iter().map(|x| format!("info.{}", x)).collect::<Vec<String>>();
    // Feature becomes info.CSQ.Feature
    let info_fields_join = args.fields_join.iter().enumerate().map(|(ind, x)| format!("{}.{}", info_fields[ind], x)).collect::<Vec<String>>();
    
    reader.reader.lines().par_bridge().for_each(|line| {
        let line = line.unwrap();
        if line.starts_with("#") {
            return;
        }
        let line = line.as_bytes() as &[u8];
        let vcf_record = VCFRecord::from_bytes(line, 1, header.clone()).unwrap();
        let variant = Variant::new(&vcf_record, header.samples(), &csq_headers);
        let explodeds = info_fields.iter().map(|x| explode_data(serde_json::to_value(&variant).unwrap(), x, &info_fields)).collect::<Vec<Vec<Map<String, Value>>>>();
        let joined = outer_join(explodeds.clone(), &info_fields_join).unwrap();
        joined.iter().filter(|x| filter_record(x, &filters)).for_each(|x| {
            let j = serde_json::to_string(&x).unwrap();
            match stdoutln!("{}", j) {
                Ok(_) => Ok(()),
                Err(e) => match e.kind() {
                    std::io::ErrorKind::BrokenPipe => std::process::exit(0),
                    _ => Err(e),
                },
            }
            .unwrap();
        });
        
        
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



pub fn outer_join(mut tables: Vec<Vec<Map<String, Value>>>, keys: &Vec<String>) -> Result<Vec<Map<String, Value>>, Box<dyn std::error::Error>> {
    if tables.len() == 0 {
        return Ok(vec![]);
    }
    if tables.len() != keys.len() {
        return Err("Number of tables should be equal to the number of keys".into());
    }
    // Create a map to hold the result of the outer join
    // right_table would be the joint table
    let mut keys = keys.iter().map(|x| x.to_string()).collect::<Vec<String>>();
    let mut right_table: Vec<Map<String, Value>> = tables.pop().unwrap();
    let right_key = keys.pop().unwrap();
    while let Some(mut left_table) = tables.pop() {
        // get all keys from the left table
        let left_firstrow = left_table[0].clone();
        let left_keys = left_firstrow.keys().collect::<Vec<&String>>();
        let right_firstrow = right_table[0].clone();
        let right_keys = right_firstrow.keys().collect::<Vec<&String>>();
        let left_key = keys.pop().unwrap();
        // Iterate over each entry in the right table, and remove  entries in the left table that are joined
        right_table.iter_mut().for_each(|right_entry| {
            let right_value = match right_entry.get(&right_key).unwrap_or(&Value::Null){
                Value::Null => "__MISSING__",
                x => x.as_str().unwrap(),
            }.split(".").collect::<Vec<&str>>();
            // fight the matching left entry, and delete the left entry from the left table
            if let Some(left_index) = left_table.iter().position(|left_entry| {
                let left_value = match left_entry.get(&left_key).unwrap_or(&Value::Null){
                    Value::Null => "__MISSING__",
                    x => x.as_str().unwrap(),
                }.split(".").collect::<Vec<&str>>();
                left_value.first().unwrap() == right_value.first().unwrap()
            }) {
                let left_entry = left_table.remove(left_index);
                for (k, v) in left_entry {
                    // skip if the key is already in the entry, and it is not Null
                    if right_entry.get(&k).unwrap_or(&Value::Null) == &Value::Null {
                        right_entry.insert(k.clone(), v.clone());
                    }
                }
            } else {
                // fill in null to all the keys in the left table that are not in the right table
                for k in left_keys.clone() {
                    if right_entry.get(k).unwrap_or(&Value::Null) == &Value::Null {
                        right_entry.insert(k.clone(), Value::Null);
                    }
                }
            }
        });
        // add the left table (whatever left) to the right table. fill missing right keys with null
        for left_entry in left_table {
            let mut new_entry = left_entry.clone();
            for k in right_keys.clone() {
                if new_entry.get(k).unwrap_or(&Value::Null) == &Value::Null {
                    new_entry.insert(k.clone(), Value::Null);
                }
            }
            right_table.push(new_entry);
        }
    }
    Ok(right_table)
}

pub fn explode_data(data:Value, key: &str, drops: &Vec<String>) -> Vec<Map<String, Value>> {
    // explode on key, but drop the keys in drops
    let mut result: Vec<Map<String, Value>> = vec![];
    match data.get(key).unwrap_or(&Value::Null) {
        Value::Array(arr) => {
            for a in arr {
                let mut new_record = data.as_object().unwrap().clone();
                for drop in drops {
                    new_record.remove(drop);
                }
                match a {
                    Value::Object(map) => {
                        for (k, v) in map {
                            let k = format!("{}.{}", key, k);
                            new_record.insert(k, v.clone());
                        }
                        result.push(new_record.clone());
                    },
                    _ => panic!("Array should contain objects"),
                }
            }
        },
        _ => {
            println!("{:?}", data);
            println!("{:?}", key);
            panic!("Data should be an array")
        },
    }
    result
}




#[cfg(test)]
mod tests {
    use super::*;
    pub fn prepare_test(fields:&Vec<String>)-> Result<(VCFReader<Box<dyn BufRead + Send + Sync>>, HashMap<String, Vec<String>>, Value), Box<dyn std::error::Error>> {
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
                csq_headers.insert(info_str.to_string(), parse_csq_header(desc));
            }
        }
        Ok((reader, csq_headers, filter))
    }
    #[test]
    fn args() {
        let args = Args::parse_from(&["test", "-i", "test/test.vcf", "-t", "1", "-f", "test/filter.yml", "--fields", "CSQ"]);
        assert_eq!(args.input, "test/test.vcf");
        assert_eq!(args.threads, 1);
        assert_eq!(args.filter, "test/filter.yml");
        assert_eq!(args.fields, vec!["CSQ".to_string()]);
    }

    #[test]
    fn test_variants() -> Result<(), Box<dyn std::error::Error>> {
        let (mut reader, csq_headers, _filter) = prepare_test(&vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut vcf_record = reader.empty_record();
        let mut variants:Vec<Variant> = Vec::new();
        while reader.next_record(&mut vcf_record).unwrap() {
            let variant = Variant::new(&vcf_record, reader.header().clone().samples(), &csq_headers);
            variants.push(variant);
            
        }
        assert!(variants.len() == 5);
        Ok(())
    }

    #[test]
    fn test_explode_data() -> Result<(), Box<dyn std::error::Error>> {
        let data = serde_json::from_str(r#"{"a":1, "b":2, "c":[{"foo":"A","bar":"B"},{"foo":"a", "bar":"b"}]}"#)?;
        let exploded = explode_data(data, "c", &vec!["b".to_string(),"c".to_string()]);
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"a":1, "c.foo":"A", "c.bar":"B"}"#).unwrap(),
            serde_json::from_str(r#"{"a":1, "c.foo":"a", "c.bar":"b"}"#).unwrap(),
        ];
        assert_eq!(exploded, expected);
        Ok(())
    }

    #[test]
    fn test_test() -> Result<(), Box<dyn std::error::Error>> {
        let (mut reader, csq_headers, filter) = prepare_test(&vec!["CSQ".to_string(), "Pangolin".to_string()])?;
        let mut vcf_record = reader.empty_record();
        let fields = vec!["info.CSQ".to_string(), "info.Pangolin".to_string()];
        let fields_join = vec!["info.CSQ.Feature".to_string(), "info.Pangolin.pangolin_transcript".to_string()];
        while reader.next_record(&mut vcf_record).unwrap() {
            let variant = Variant::new(&vcf_record, reader.header().clone().samples(), &csq_headers);
            //let val = serde_json::to_value(&variant)?;
            //println!("{:?}", serde_json::to_string(&variant)?);
            let explodeds = fields.iter().map(|x| explode_data(serde_json::to_value(&variant).unwrap(), x, &fields)).collect::<Vec<Vec<Map<String, Value>>>>();
            let joined = outer_join(explodeds.clone(), &fields_join)?;
            println!("====================");
            //println!("{}", serde_json::to_string_pretty(&joined)?);
            let filtered_record = joined.iter().filter(|x| filter_record(x, &filter)).collect::<Vec<&Map<String, Value>>>();
            println!("{}", serde_json::to_string_pretty(&filtered_record)?);
        }
        //println!("{:?}", serde_json::to_string(&variants)?);
        Ok(())
    }
    #[test]
    fn test_outer_join1() -> Result<(), Box<dyn std::error::Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key3": "value6"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key3": "value8"}"#).unwrap(),
        ];
        let joined_table = outer_join(vec![table1, table2], &vec!["key1".to_string(), "key1".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2", "key3": "value6"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key2": null, "key3": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4", "key3": null}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
    #[test]
    fn test_outer_join2() -> Result<(), Box<dyn std::error::Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key3": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key3": "value8"}"#).unwrap(),
        ];
        let joined_table = outer_join(vec![table1, table2], &vec!["key1".to_string(), "key1".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": "value2", "key3": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value7", "key2": null, "key3": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4", "key3": null}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
    #[test]
    fn test_outer_join3() -> Result<(), Box<dyn std::error::Error>> {
        let table1: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key2": null}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key2": "value4"}"#).unwrap(),
        ];
        let table2: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key11": "value1.8", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key11": "value7", "key2": "value8"}"#).unwrap(),
        ];
        let joined_table = outer_join(vec![table1, table2], &vec!["key1".to_string(), "key11".to_string()])?;
        let expected: Vec<Map<String, Value>> = vec![
            serde_json::from_str(r#"{"key1": "value1", "key11": "value1.8", "key2": "value2"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": null, "key11": "value7", "key2": "value8"}"#).unwrap(),
            serde_json::from_str(r#"{"key1": "value3", "key11": null, "key2": "value4"}"#).unwrap(),
        ];
        assert_eq!(joined_table, expected);
        Ok(())
    }
}