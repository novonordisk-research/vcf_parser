use serde::{Deserialize, Serialize};
use serde_json::{Map, Number, Value};
use vcf::VCFRecord;
use std::collections::HashMap;
use std::str;

#[derive(Serialize, Deserialize)]
pub struct Variant {
    pub chromosome: String,
    pub position: u64,
    pub id: Vec<String>,
    pub reference: String,
    pub alternative: Vec<String>,
    pub qual: Option<f64>,
    pub filter: Vec<String>,
    pub info: Map<String, Value>,
    pub genotype: Map<String, Value>,
}

impl Variant {
    pub fn new(
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
                        if val == "." {
                            Value::Null
                        } else {
                            Value::Number(Number::from(val.parse::<i64>().unwrap_or_else( |_| panic!("parse error {f_k}, {val}"))))
                        },
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
                            if str::from_utf8(&dat[0]).unwrap() == "." {
                                Value::Null
                            } else {
                                Value::from(str::from_utf8(&dat[0]).unwrap().parse::<f64>().unwrap_or_else(|_| panic!("parse error {}, {}", field_str, str::from_utf8(&dat[0]).unwrap())))
                            }
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
    pub fn get_value(self: &Self, key: &str) -> Value{
        // Helper function to recursively search for the key
        fn search<'a>(value: &'a Value, key: &str, results: &mut Vec<&'a Value>) {
            match value {
                Value::Object(map) => {
                    for (k, v) in map {
                        if k == key {
                            results.push(v);
                        } else {
                            search(v, key, results);
                        }
                    }
                }
                Value::Array(arr) => {
                    for v in arr {
                        search(v, key, results);
                    }
                }
                _ => {}
            }
        }
        let mut results = Vec::new();
        let bind = Value::Object(self.info.clone());
        search(&bind, key, &mut results);
        // If only one result, return it directly, otherwise return an array of results
        if results.len() == 1 {
            results[0].clone()
        } else {
            Value::Array(results.into_iter().cloned().collect())
        }
    }
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
