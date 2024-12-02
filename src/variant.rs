use serde::{Deserialize, Serialize};
use serde_json::{Map, Number, Value};
use vcf::VCFRecord;
use std::collections::HashMap;
use std::str;
use crate::utils;
use crate::error::VcfParserError;
use crate::Result;

serde_with::with_prefix!(prefix_info "info.");

#[derive(Serialize, Deserialize)]
pub struct Variant {
    pub chromosome: String,
    pub position: u64,
    pub id: String,
    pub reference: String,
    pub alternative: String,
    pub qual: Option<f64>,
    pub filter: String,
    #[serde(flatten, with="prefix_info")]
    pub info: Map<String, Value>,
    pub genotype: Map<String, Value>,
}

impl Variant {
    pub fn new(
        vcf_record: &VCFRecord,
        samples: &[Vec<u8>],
        csq_headers: &HashMap<String, Vec<String>>,
    ) -> Result<Self> {
        // parse genotype
        let mut genotype = Map::new();
        for sample in samples {
            let mut sample_genotype = Map::new();
            for key in &vcf_record.format {
                let genotype_raw = vcf_record.genotype(&sample, &key);
                let f_k = str::from_utf8(key)?;
                let val = str::from_utf8(&genotype_raw.ok_or(VcfParserError::VariantParseError("genotype parse error"))?[0])?;
                match vcf_record.header().format(key).ok_or(VcfParserError::VariantParseError("cannot find header"))?.value_type {
                    vcf::ValueType::Integer => sample_genotype.insert(
                        f_k.to_string(),
                        if val == "." {
                            Value::Null
                        } else {
                            Value::Number(Number::from(val.parse::<i64>()?))
                        },
                    ),
                    vcf::ValueType::Float => sample_genotype.insert(
                        f_k.to_string(),
                        Value::Number(Number::from_f64(val.parse::<f64>()?).ok_or(VcfParserError::VariantParseError("parse error"))?),
                    ),
                    _ => sample_genotype.insert(f_k.to_string(), Value::String(val.to_string())),
                };
            }
            genotype.insert(
                str::from_utf8(&sample)?.to_string(),
                Value::Object(sample_genotype),
            );
        }
        // parse info
        let mut info = Map::new();
        for id in vcf_record.header().info_list() {
            let field = vcf_record.header().info(id).ok_or_else(|| VcfParserError::VariantParseError("field not found"))?;
            let field_str = str::from_utf8(&field.id)?;
            let val = match vcf_record.info(id) {
                Some(dat) => {
                    if *field.value_type == vcf::ValueType::Flag {
                        // flag type
                        Value::Bool(true)
                    } else if *field.value_type == vcf::ValueType::Integer {
                        if *field.number == vcf::Number::Allele {
                            // assume input is normalised vcf. not care about the number of alleles.
                            // just take the first element.
                            Value::from(
                                dat.iter()
                                    .map(|x| {
                                        if str::from_utf8(x)? == "." {
                                            Ok(Value::Null)
                                        } else {
                                            Ok(Value::from(Number::from(
                                                str::from_utf8(x)?.parse::<i64>()?,
                                            )))
                                        }
                                    })
                                    .collect::<Result<Vec<Value>>>()?[0].to_owned(),
                            )
                        } else {
                            if dat.len() == 0 || str::from_utf8(&dat[0])? == "." {
                                Value::Null
                            } else {
                                Value::from(str::from_utf8(&dat[0])?.parse::<i64>()?)
                            }
                        }
                    } else if *field.value_type == vcf::ValueType::Float {
                        if *field.number == vcf::Number::Allele {
                            Value::from(
                                dat.iter()
                                    .map(|x| {
                                        if str::from_utf8(x)? == "." {
                                            Ok(Value::Null)
                                        } else {
                                            Ok(Value::from(Number::from_f64(
                                                str::from_utf8(x)?.parse::<f64>()?,
                                            )))
                                        }
                                    })
                                    .collect::<Result<Vec<Value>>>()?[0].to_owned(),
                            )
                        } else {
                            if dat.len() == 0 || str::from_utf8(&dat[0])? == "." {
                                Value::Null
                            } else {
                                Value::from(str::from_utf8(&dat[0])?.parse::<f64>()?)
                            }
                        }
                    } else if csq_headers.contains_key(field_str) {
                        dat.iter()
                            .map(|csq_field| {
                                let csq_vec = csq_headers
                                    .get(field_str)
                                    .ok_or(VcfParserError::VariantParseError("field not found"))?
                                    .iter()
                                    .zip(str::from_utf8(csq_field)?.split("|"));
                                let mut csq = Map::new();
                                for (k, v) in csq_vec {
                                    csq.insert(k.to_string(), utils::try_parse_number(v));
                                }
                                Ok(Value::Object(csq))
                            })
                            .collect::<Result<Value>>()?
                    } else {
                        if *field.number == vcf::Number::Allele {
                            Value::from(
                                dat.iter()
                                    .map(|x| Ok(Value::String(str::from_utf8(x)?.to_string())))
                                    .collect::<Result<Vec<Value>>>()?[0].clone(),
                            )
                        } else {
                            if dat.len() == 0 {
                                Value::Null
                            } else {
                                Value::String(str::from_utf8(&dat[0])?.to_string())
                            }
                        }
                    }
                }
                None => {
                    match *field.value_type {
                        // flag type
                        vcf::ValueType::Flag => Value::Bool(false),
                        // if csq_header contains the field, parse it as csq. Values are all null.
                        _ if csq_headers.contains_key(field_str) => Value::Array(vec![csq_headers
                            .get(field_str)
                            .ok_or(VcfParserError::VariantParseError("field not found"))?
                            .iter()
                            .map(|k| (k.to_string(), Value::Null))
                            .collect::<Map<String, Value>>()
                            .into()]),
                        _ => Value::Null
                    }
                }
            };
            info.insert(str::from_utf8(&field.id)?.to_string(), val);
        }
        Ok(Variant {
            chromosome: str::from_utf8(&vcf_record.chromosome)?.to_string(),
            position: vcf_record.position,
            // flatten id / filter / alternative
            id: vcf_record
                .id
                .iter()
                .map(|x| Ok(str::from_utf8(&x)?.to_string()))
                .collect::<Result<Vec<String>>>()?
                .join(";"),
            alternative: vcf_record
                .alternative
                .iter()
                .map(|x| Ok(str::from_utf8(&x)?.to_string()))
                .collect::<Result<Vec<String>>>()?
                .join(","),
            filter: vcf_record
                .filter
                .iter()
                .map(|x| Ok(str::from_utf8(&x)?.to_string()))
                .collect::<Result<Vec<String>>>()?
                .join(","),
            reference: str::from_utf8(&vcf_record.reference)?.to_string(),
            qual: vcf_record.qual,
            info,
            genotype: genotype,
        })
    }
}

