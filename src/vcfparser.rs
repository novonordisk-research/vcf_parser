use super::{OutputFormat, error::VcfParserError};
use std::collections::HashMap;
use std::str;
use anyhow::Result;
use std::io::BufRead;
use vcf::VCFReader;
use std::sync::Arc;
use serde_json;
use crate::utils;
pub struct VcfParser<T>
where T: BufRead + Send + Sync,
{
    /// filter to use for filtering variants
    pub filters: serde_json::Value,
    /// info fields have `info.` prefix, such as info.CSQ, info.VEP
    pub info_fields: Vec<String>,
    /// fields to join on, such as Feature,Transcript_id. Version numbers will be ignored.
    pub fields_join: Vec<String>,
    /// output format, tsv, json, vcf(coming soon)
    pub output_format: OutputFormat,
    /// reader to read from
    pub reader: VCFReader<T>,
    /// vcf header
    pub header: Arc<vcf::VCFHeader>,
    /// info headers
    /// CSQ headers
    pub csq_headers: Arc<HashMap<String, Vec<String>>>,
    /// tsv headers
    pub tsv_headers: Vec<String>,
}
impl <T> VcfParser <T>
where T: BufRead + Send + Sync,
{
    pub fn new(
        filters: serde_json::Value,
        fields: Vec<String>,
        fields_join: Vec<String>,
        columns: Option<Vec<String>>,
        output_format: OutputFormat,
        reader: T,
    ) -> Result<Self> {
        if fields.len() >1 && fields.len() != fields_join.len() {
            return Err(VcfParserError::InvalidArgument("Number of fields should be equal to the number of fields_join".into()).into());
        }
        let info_fields = fields.iter().map(|x| format!("info.{}", x)).collect::<Vec<String>>();
        let fields_join = fields_join.iter().enumerate().map(|(ind, x)| format!("{}.{}", info_fields[ind], x)).collect::<Vec<String>>();
        let mut info_headers: Vec<String> = Vec::new();
        let mut csq_headers: HashMap<String, Vec<String>> = HashMap::new();
        let reader = VCFReader::new(reader)?;
        let header = Arc::new(reader.header().clone());
        for info in header.info_list() {
            let info_str = str::from_utf8(&info)?;
            info_headers.push(info_str.to_string());
            let desc = str::from_utf8(reader.header().info(info).unwrap().description).unwrap();
            if fields.contains(&info_str.to_string()) {
                csq_headers.insert(info_str.to_string(), utils::parse_csq_header(desc));
            }
        }
        let info_headers = Arc::new(info_headers);
        let csq_headers = Arc::new(csq_headers);
        let tsv_headers = utils::get_output_header(&info_headers, &csq_headers, &columns);
        Ok(VcfParser {
            filters,
            info_fields,
            fields_join,
            output_format,
            reader,
            csq_headers,
            tsv_headers,
            header,
        })
    }

}