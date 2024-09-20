use std::error::Error;
use calm_io::stdoutln;
use serde_json::{json, Map, Value};
use regex::Regex;

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

pub fn logical_expression_to_json(expression: &str) -> Value {
    let re = Regex::new(r"(\w+)\s*(==|<=|>=|<|>|in)\s*(\([\w.,\s]+\)|[\w.]+)").unwrap();

    let mut json_obj = json!({});

    if expression.contains("and") {
        let parts: Vec<&str> = expression.split("and").collect();
        let mut and_arr = vec![];
        for part in parts {
            and_arr.push(logical_expression_to_json(part.trim()));
        }
        json_obj["AND"] = Value::Array(and_arr);
    } else if expression.contains("or") {
        let parts: Vec<&str> = expression.split("or").collect();
        let mut or_arr = vec![];
        for part in parts {
            or_arr.push(logical_expression_to_json(part.trim()));
        }
        json_obj["OR"] = Value::Array(or_arr);
    } else {
        if let Some(captures) = re.captures(expression) {
            let name = captures.get(1).unwrap().as_str();
            let op = captures.get(2).unwrap().as_str();
            let value = captures.get(3).unwrap().as_str();
            if op == "in" {
                let value_list: Vec<&str> = value
                    .trim_matches(|c| c == '(' || c == ')')
                    .split(',')
                    .map(|s| s.trim())
                    .collect();
                json_obj = json!({
                    "name": name,
                    "op": op,
                    "value": value_list
                });
            } else {
                json_obj = json!({
                    "name": name,
                    "op": op,
                    "value": value
                });
            }
        }
    }

    json_obj
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