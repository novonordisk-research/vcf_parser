use std::error::Error;
use std::collections::HashMap;
use calm_io::stdoutln;
use serde_json::{Map, Value, Number};

pub fn print_line_to_stdout(line: &str) -> Result<(), Box<dyn Error>> {
    // output line to stdout.
    // also captures broken pipe that might be caused by | head, for example.
    match stdoutln!("{}", line) {
        Ok(_) => Ok(()),
        Err(e) => match e.kind() {
            std::io::ErrorKind::BrokenPipe => std::process::exit(0),
            _ => Err(Box::new(e)),
        },
    }
}

pub fn parse_csq_header(header: &str) -> Vec<String> {
    // parse the csq header to produce a list of fields
    header.split(": ").collect::<Vec<&str>>()[1]
        .split("|")
        .map(|x| x.trim().to_string())
        .collect::<Vec<String>>()
}

pub fn filter_record(record: &Map<String,Value>, filters: &Value) -> bool {
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

    One can potentially convert this into logic expression and parse through coolrule, but it would be many times slower!
    */
    match filters {
        Value::Object(map) => {
            for (k, v) in map {
                if k.eq_ignore_ascii_case("AND") {
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
                } else if k.eq_ignore_ascii_case("OR") {
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
                        "eq" | "=" | "==" => {
                            if *val == *value {
                                return true;
                            } 
                            return false;    
                        }
                        "ne" | "!=" | "≠" => {
                            if *val != *value {
                                return true;
                            }
                            return false;
                        }
                        "gt" | ">" => {
                            if val.is_null() {
                                return false;
                            }
                            if val.as_f64().unwrap() > value.as_f64().unwrap() {
                                return true;
                            }
                            return false;
                        }
                        "ge" | ">=" | "≥" => {
                            if val.is_null() {
                                return false;
                            }
                            if val.as_f64().unwrap() >= value.as_f64().unwrap() {
                                return true;
                            } 
                            return false;
                        }
                        "lt" | "<" => {
                            if val.is_null() {
                                return false;
                            }
                            if val.as_f64().unwrap() < value.as_f64().unwrap() {
                                return true;
                            }
                            return false;
                        }
                        "le" | "<=" | "≤" => {
                            if val.is_null() {
                                return false;
                            }
                            if val.as_f64().unwrap() <= value.as_f64().unwrap() {
                                return true;
                            } 
                            return false;
                        }
                        "in" | "∈" => {
                            match value {
                                serde_json::Value::Array(arr) => {
                                    return arr.contains(val)
                                }
                                _ => {
                                    if val == value {
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

pub fn outer_join(mut tables: Vec<Vec<Map<String, Value>>>, keys: &Vec<String>) -> Result<Vec<Map<String, Value>>, Box<dyn Error>> {
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
                        right_entry.insert(k.to_string(), v);
                    }
                }
            } else {
                // fill in null to all the keys in the left table that are not in the right table
                for k in left_keys.clone() {
                    if right_entry.get(k).unwrap_or(&Value::Null) == &Value::Null {
                        right_entry.insert(k.to_string(), Value::Null);
                    }
                }
            }
        });
        // add the left table (whatever left) to the right table. fill missing right keys with null
        for left_entry in left_table {
            let mut new_entry = left_entry;
            for k in right_keys.clone() {
                if new_entry.get(k).unwrap_or(&Value::Null) == &Value::Null {
                    new_entry.insert(k.to_string(), Value::Null);
                }
            }
            right_table.push(new_entry);
        }
    }
    Ok(right_table)
}

pub fn explode_data(data:Value, key: &str, drops: &Vec<String>) -> Vec<Map<String, Value>> {
    // explode on key, but drop the columns in `drops`.
    // dropping columns is desirable if one row is going to be exploded on multiple keys, and it will result in duplicated columns
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
                        result.push(new_record);
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

pub fn get_row(data:&Map<String, Value>, header:&Vec<String>) -> Vec<String> {
    // given a data and a header, return a row
    // for tsv output
    header.iter().map(|x| {
        match data.get(x).unwrap_or(&Value::Null) {
            Value::Null => "".to_string(),
            Value::Number(n) => n.to_string(),
            Value::String(s) => s.as_str().to_string(),
            x => x.to_string(),
        }
    }).collect::<Vec<String>>()
}

pub fn get_output_header<'a>(info_header: &Vec<String>, csq_header: &HashMap<String, Vec<String>>, samples: &[Vec<u8>], user_columns: &Option<Vec<String>>) -> Vec<String> {
    // get the header for csv output
    // essential columns are in the front. All info columns are in the back, sorted alphabetically.
    let mut header: Vec<String> = Vec::new();
    for h in info_header {
        if csq_header.contains_key(h) {
            for csq in csq_header.get(h).unwrap() {
                header.push(format!("info.{}.{}", h, csq));
            }
        } else {
            header.push(format!("info.{}", h));
        }
    }
    // sort header
    header.sort();

    // turn samples into strings
    let samples = samples.iter().map(|x| std::str::from_utf8(x).unwrap().to_string()).collect::<Vec<String>>();

    // add chromosome, position, id, ref, alt, qual, filter
    let header = vec![
        "chromosome".to_string(),
        "position".to_string(),
        "id".to_string(),
        "reference".to_string(),
        "alternative".to_string(),
        "qual".to_string(),
        "filter".to_string(),
    ].iter()
    .chain(header.iter())
    .chain(samples.iter())
    .map(|x| x.to_string()).collect::<Vec<String>>();

    // if user_columns is None, then all columns are selected
    // validate if user_columns are a subset of header
    match user_columns {
        Some(user_columns) => {
            for user_column in user_columns {
                if !header.contains(&user_column) {
                    panic!("Column {} is not in the header", user_column);
                }
            }
            user_columns.clone()
        },
        None => header,
    }
}

pub fn try_parse_number(input: &str) -> Value {
    // Can't just rely on VCF header to parse info fields if there are nested ones like CSQ-like fields, 
    // as their data types are not necessarily exposed in the VCF header.
    match input {
        inp if inp.contains('.') || input.contains('e') || input.contains('E') => {
            // Try to parse as f64
            match input.parse::<f64>() {
                Ok(num) => Value::Number(Number::from_f64(num).unwrap()),
                Err(_) => Value::String(input.to_string()),
            }
        },
        "" => Value::Null,
        _ =>  match input.parse::<i64>() {
            Ok(num) => Value::Number(Number::from(num)),
            Err(_) => Value::String(input.to_string()),
        }
    }
}