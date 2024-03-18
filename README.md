A VCF Parser written in Rust
============================
It reads a VCF file (plain or gzipped) and parse the output to a JSON format.
It is multiallelic aware and CSQ aware, and use threads to speed up the process.

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
vcf_parser -i <vcf_file>
```
Or from stdin
```bash
zcat <vcf.gz> | vcf_parser
```

By default it uses all threads available, but you can specify the number of threads to use with the `-t` option.

```bash
# Use 4 threads
vcf_parser -i <vcf_file> -t 4
```

Parsing a densily annotated VCF file in the CSQ field takes about 40 seconds for 100K lines (~2500 lines per second), which is about twice as fast as the python version.
~20 seconds for 100K lines with 2 threads.
~10 seconds for 100K lines with 4 threads.
~8 seconds for 100K lines with 8 threads.

Output
------
It output a json object per VCF line. This enables one to use tools like `jq` to filter the output.

```bash
vcf_parser -i <vcf_file> | jq -f filter.jq
```
And in a filter.jq file:
```jq
select (
  (
    (any(.info.CSQ[]; ."gnomAD.genomes_4.0_AF" // 0 | . <= 0.04)) and
    (any(.info.CSQ[]; ."gnomAD.exomes_4.0_AF" // 0 | . <= 0.04))
  ) and (
    (any(.info.CSQ[]; .ESM1b // 0 | . <= -7.5)) or
    (any(.info.CSQ[]; .CADD_PHRED // 0 | . >= 20)) or
    (any(.info.CSQ[]; .LoF // "" | . == "HC" )) or
    (any(.info.CSQ[]; .LoF // "" | . == "LC" ))
  )
) | .info.CSQ |= map(with_entries(select(.value != null)))
```
The `.info.CSQ |= map(with_entries(select(.value != null)))` line is used to remove null values from the CSQ field.
An example of the output
[!jsonexample.png](./asset/json_example.png)
