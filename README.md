A VCF Parser written in Rust
============================
It reads a VCF file (plain or gzipped) and parse the output to a JSON/TSV format.
The input VCF needs to be normalised so that there is no more than one alt per line. By default it uses all cores available to process the input.

Features
------------
* Unnest INFO field
* Explode CSQ-like fields
* Join exploded fields by transcript_id
* Accept arbitrarily sophisticated filters defined in a yaml file.
* Multithreaded
* Output in JSON or TSV format
* Support logic expression for filters (experimental)

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

An example of `filter.yml` can be found in the `test/` folder. You can replace `eq` with `=` or `==`, `le` with `<=` or `≤`, etc.
Equivalently, you can now pass a logic expression as a string to the `-f` option. For example:
```bash
vcf_parser -i test/test.vcf -f "(info.AF <= 0.01 AND info.CSQ.IMPACT in (HIGH,MODERATE)) AND (info.CADD_PHRED >=20 OR info.Pangolin.pangolin_max_score >= 0.5 or info.Pangolin.pangolin_max_score <= -0.5)" --fields CSQ,Pangolin --fields-join Feature,pangolin_transcript
```
**NOTE**: As in the current implementation, "AND" and "OR" have the same precedence. Please use parentheses to make the logic expression unambiguous.

### Options
```
-h #help 
-i <input.vcf[.gz]>  
-f <filter.yaml or expression>
-t <thread number>
-l # to list columns and exit
-c <columns to output>
--output-format <j|t> #j for json, t for tsv
--fields #fields to explode. default to CSQ
--fields-join #keys to join fields, in the same order
```
