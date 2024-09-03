A VCF Parser written in Rust
============================
It reads a VCF file (plain or gzipped) and parse the output to a JSON/TSV format.
The input VCF needs to be normalised so that there is no more than one alt per line. By default it uses all cores available to process the input.

Change log:
------------
### 0.2.0
* Breaking change: Flatten the output json, so there is no longer nested structure in the output.
* Breaking change: Transcript aware filtering. If you want canonical trancripts, you only get canonical transcripts.
* Provide TSV as an option (also set as default).
* Join different CSQ-like fields with keys. For instance if you have two CSQ-like fields, they can be joined by transcript id. vcf_parser ignores the version number when it joins records.



Installation
------------
Assume you have `cargo` installed, you can install it with:
```bash
cargo install --path .
```
There's a binary on `marjorie` already in `/nfs_home/projects/departments/nnrco/genetic_department/bin`

Usage
------
```bash
vcf_parser -i test/test.vcf >output.tsv
```
Or from stdin
```bash
cat test/test.vcf | vcf_parser >output.tsv
```

Parse multiple CSQ-like fields and join by transcript_id, and filter variants as defined in `filter.yml`, output json format
```bash
vcf_parser -i test/test.vcf --fields CSQ,Pangolin --fields-join Feature,pangolin_transcript -f filter.yml --output-format j >output.json
```

### Options
```
-h #help 
-i <input.vcf[.gz]>  
-f <filter.yaml>
--output-format <j|t> #j for json, t for tsv
--fields #fields to explode. default to CSQ
--fields-join #keys to join fields, in the same order
```
