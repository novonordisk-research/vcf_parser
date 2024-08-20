// Read from stdin or a .vcf[.gz] file
use clap::Parser;
use vcf_parser::Args;
fn main()  {
    let args = Args::parse();
    if let Err(e) = vcf_parser::run(args) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}