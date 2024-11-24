use super::{OutputFormat, error::VcfParserError};
use std::collections::HashMap;
use std::str;
use anyhow::Result;
use std::io::BufRead;
use vcf::VCFReader;
use serde_json;
use crate::utils;
pub struct VcfParser {
    /// filter to use for filtering variants
    pub filters: serde_json::Value,
    /// fields to parse, such as CSQ,VEP
    fields: Vec<String>,
    /// info fields have `info.` prefix, such as info.CSQ, info.VEP
    pub info_fields: Vec<String>,
    /// fields to join on, such as Feature,Transcript_id. Version numbers will be ignored.
    pub fields_join: Vec<String>,
    /// columns to output
    columns: Option<Vec<String>>,
    /// output format, tsv, json, vcf(coming soon)
    pub output_format: OutputFormat,
    /// reader to read from
    pub reader: Option<VCFReader<Box<dyn BufRead + Send + Sync>>>,
    /// vcf header
    header: Option<vcf::VCFHeader>,
    /// info headers
    info_headers: Option<Vec<String>>,
    /// CSQ headers
    csq_headers: Option<HashMap<String, Vec<String>>>,
    /// tsv headers
    pub tsv_headers: Option<Vec<String>>,
}
impl VcfParser{
    pub fn new(
        filters: serde_json::Value,
        fields: Vec<String>,
        fields_join: Vec<String>,
        columns: Option<Vec<String>>,
        output_format: OutputFormat,
    ) -> Result<Self> {
        if fields.len() >1 && fields.len() != fields_join.len() {
            return Err(VcfParserError::InvalidArgument("Number of fields should be equal to the number of fields_join".into()).into());
        }
        let info_fields = fields.iter().map(|x| format!("info.{}", x)).collect::<Vec<String>>();
        let fields_join = fields_join.iter().enumerate().map(|(ind, x)| format!("{}.{}", info_fields[ind], x)).collect::<Vec<String>>();
        Ok(VcfParser {
            filters,
            fields,
            info_fields,
            fields_join,
            columns,
            output_format,
            reader: None,
            header: None,
            info_headers: None,
            csq_headers: None,
            tsv_headers: None,
        })
    }
    /// asign reader to self.reader
    pub fn with_reader(mut self, reader: Box<dyn BufRead + Send + Sync>) -> Result<Self> {
        
        self.reader = Some(VCFReader::new(reader)?);
        Ok(self)
    }
    /// get headers from reader
    pub fn with_headers(mut self) -> Result<Self> {
        let reader = self.reader.as_mut().unwrap();
        let mut info_headers: Vec<String> = Vec::new();
        let mut csq_headers: HashMap<String, Vec<String>> = HashMap::new();
        for info in reader.header().info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(info_str.to_string());
            let desc = str::from_utf8(reader.header().info(info).unwrap().description).unwrap();
            if self.fields.contains(&info_str.to_string()) {
                csq_headers.insert(info_str.to_string(), utils::parse_csq_header(desc));
            }
        }
        self.header = Some(reader.header().clone());
        self.info_headers = Some(info_headers.clone());
        self.csq_headers = Some(csq_headers.clone());
        let tsv_headers = utils::get_output_header(info_headers, &csq_headers, &self.columns);
        self.tsv_headers = Some(tsv_headers);
        Ok(self)
    }

    pub fn into_reader(& mut self) -> VCFReader<Box<dyn BufRead + Send + Sync>> {
        self.reader.take().unwrap()
    }

    pub fn header(&self) -> vcf::VCFHeader {
        self.header.as_ref().unwrap().to_owned()
    }

    pub fn csq_headers(&self) -> HashMap<String, Vec<String>> {
        self.csq_headers.as_ref().unwrap().to_owned()
    }
}