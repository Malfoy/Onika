#[macro_use]
extern crate clap;
extern crate flate2;
extern crate rayon;
extern crate zstd;
extern crate f128;
extern crate num_traits;
extern crate needletail;
extern crate byteorder;

mod onika_index;
mod utils;

use clap::{App, Arg, ArgGroup};
use onika_index::{Index, IndexBuilder, SketchMode};
use std::time::Instant;
use std::fs::File;
use std::io::{BufReader, BufRead};
use needletail::parse_fastx_file;

/// Estimates the average document size by sampling the first 100 documents or sequences.
fn estimate_genome_size(matches: &clap::ArgMatches) -> u32 {
    const SAMPLE_SIZE: usize = 100;
    let mut total_len: u64 = 0;
    let mut doc_count: u64 = 0;

    println!("Sampling up to {} documents/sequences to estimate average size...", SAMPLE_SIZE);

    if let Some(fof_path) = matches.value_of("index_fof") {
        let file = File::open(fof_path).expect("Could not open file of files for estimation");
        let reader = BufReader::new(file);

        'outer: for filepath_res in reader.lines() {
            if let Ok(filepath) = filepath_res {
                if let Ok(mut file_reader) = parse_fastx_file(&filepath) {
                    let mut current_doc_size: u64 = 0;
                    while let Some(Ok(record)) = file_reader.next() {
                        current_doc_size += record.seq().len() as u64;
                    }

                    if current_doc_size > 0 {
                        total_len += current_doc_size;
                        doc_count += 1;
                        if doc_count >= SAMPLE_SIZE as u64 {
                            break 'outer;
                        }
                    }
                }
            }
        }
    } else if let Some(single_file_path) = matches.value_of("index_by_line") {
        if let Ok(mut reader) = parse_fastx_file(single_file_path) {
            while let Some(Ok(record)) = reader.next() {
                if doc_count >= SAMPLE_SIZE as u64 {
                    break;
                }
                total_len += record.seq().len() as u64;
                doc_count += 1;
            }
        }
    }

    if doc_count > 0 {
        (total_len / doc_count) as u32
    } else {
        5_000_000
    }
}


fn main() {
    let matches = App::new("Onika-rs")
        .version("1.0")
        .author("Rust Translation")
        .about("A Rust implementation of the Onika indexing tool")
        .arg(Arg::with_name("index_fof")
            .short("I")
            .long("index-fof")
            .value_name("FILE")
            .help("Index a file of files, where each file is a document.")
            .takes_value(true))
        .arg(Arg::with_name("index_by_line")
            .short("i")
            .long("index-by-line")
            .value_name("FILE")
            .help("Index a file where each sequence is a document.")
            .takes_value(true))
        .group(ArgGroup::with_name("index_mode")
            .args(&["index_fof", "index_by_line"]))
        .arg(Arg::with_name("dump")
            .short("d")
            .long("dump")
            .value_name("FILE")
            .help("File to dump the index to (default: onika_index.bin)")
            .takes_value(true))
        .arg(Arg::with_name("load")
            .short("l")
            .long("load")
            .value_name("FILE")
            .help("Load index from a binary file")
            .takes_value(true))
        .arg(Arg::with_name("query_fof")
            .short("Q")
            .long("query-fof")
            .value_name("FILE")
            .help("Query with a file of files, where each file is a query.")
            .takes_value(true))
        .arg(Arg::with_name("query_by_line")
            .short("q")
            .long("query-by-line")
            .value_name("FILE")
            .help("Query with a file where each line is a query.")
            .takes_value(true))
        .group(ArgGroup::with_name("query_mode")
            .args(&["query_fof", "query_by_line"]))
        .arg(Arg::with_name("sketch_mode")
            .long("sketch-mode")
            .value_name("MODE")
            .help("Sets the sketching mode [default|perfect].")
            .default_value("default")
            .takes_value(true))
        .arg(Arg::with_name("matrix")
            .long("matrix")
            .help("Output results in matrix format instead of the default sparse format."))
        .arg(Arg::with_name("zstd_level")
            .long("zstd-level")
            .value_name("LEVEL")
            .help("Compression level for zstd (0 for none, default 1).")
            .takes_value(true))
        .arg(Arg::with_name("print_stats")
            .long("print-stats")
            .help("Prints statistics about the index after construction."))
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .value_name("FILE")
            .help("File where the results are written (default: output.txt)")
            .takes_value(true))
        .arg(Arg::with_name("w")
            .short("w")
            .long("w_size")
            .value_name("INT")
            .help("Fingerprint size (default: 12)")
            .takes_value(true))
        .arg(Arg::with_name("s")
            .short("s")
            .long("s_size")
            .value_name("INT")
            .help("Set sketch size to 2^S (default: 15)")
            .takes_value(true))
        .arg(Arg::with_name("k")
            .short("k")
            .long("k_size")
            .value_name("INT")
            .help("K-mer size (default: 31)")
            .takes_value(true))
        .arg(Arg::with_name("e")
            .short("e")
            .long("e_size")
            .value_name("INT")
            .help("Expected genome size. If not provided in perfect mode, it will be estimated.")
            .takes_value(true))
        .arg(Arg::with_name("threads")
            .long("threads")
            .value_name("INT")
            .help("Number of threads used (default: 32)")
            .takes_value(true))
        .arg(Arg::with_name("threshold")
            .long("threshold")
            .value_name("FLOAT")
            .help("Threshold for query results (default: 0.5)")
            .takes_value(true))
        .arg(Arg::with_name("sub_div_factor")
            .long("sub-div-factor")
            .value_name("B")
            .help("Sub-division factor for memory optimization (power of 2).")
            .default_value("4")
            .takes_value(true))
        .get_matches();

    if !matches.is_present("load") && !matches.is_present("index_mode") {
        println!("Error: An indexing mode (-I or -i) is required when not loading an index with -l.");
        return;
    }

    let w = value_t!(matches, "w", u32).unwrap_or(16);
    let s = value_t!(matches, "s", u32).unwrap_or(16);
    let k = value_t!(matches, "k", u32).unwrap_or(31);
    let threads = value_t!(matches, "threads", u16).unwrap_or(32);
    let threshold = value_t!(matches, "threshold", f64).unwrap_or(0.0);
    let output = matches.value_of("output").unwrap_or("output.txt");
    let is_matrix = matches.is_present("matrix");
    let zstd_level = value_t!(matches, "zstd_level", i32).unwrap_or(1);
    
    let sketch_mode = match matches.value_of("sketch_mode").unwrap() {
        "perfect" => SketchMode::Perfect,
        _ => SketchMode::Default,
    };

    let sub_div_factor = value_t!(matches, "sub_div_factor", u32).unwrap_or(1);

    let e = if sketch_mode == SketchMode::Perfect && !matches.is_present("e") {
        println!("Perfect fingerprinting mode enabled without -e flag. Estimating genome size...");
        let estimated_e = estimate_genome_size(&matches);
        println!("Estimated genome size (e): {}", estimated_e);
        estimated_e
    } else {
        value_t!(matches, "e", u32).unwrap_or(5_000_000)
    };

    rayon::ThreadPoolBuilder::new().num_threads(threads as usize).build_global().unwrap();

    let index = if let Some(load_file) = matches.value_of("load") {
        let start_loading = Instant::now();
        println!("ONIKA INDEX LOADING from {}", load_file);
        let index = Index::from_file(load_file).expect("Failed to load index");
        println!("Loading took {} seconds.", start_loading.elapsed().as_secs());
        index
    } else {
        let start_indexing = Instant::now();
        let builder = IndexBuilder::new(s, k, w, e, sketch_mode, sub_div_factor,1024);
        
        if let Some(fof) = matches.value_of("index_fof") {
            println!("Indexing file of files: {}", fof);
            builder.index_file_of_files(fof);
        } else if let Some(line_file) = matches.value_of("index_by_line") {
            println!("Indexing file by line: {}", line_file);
            builder.index_file_line_by_line(line_file);
        }

        let final_index = builder.into_final_index();
        println!("Indexing took {} seconds.", start_indexing.elapsed().as_secs());

        if let Some(dump_file) = matches.value_of("dump") {
             println!("Dumping index to file: {}", dump_file);
             if let Err(e) = final_index.dump_index_disk(dump_file, zstd_level) {
                 eprintln!("Error dumping index: {}", e);
             }
        }
        final_index
    };

    if matches.is_present("print_stats") {
        index.print_stats();
    }

    if let Some(fof) = matches.value_of("query_fof") {
        println!("\nONIKA QUERY with file of files: {}", fof);
        let start_query = Instant::now();
        index.query_file_of_file(fof, threshold, is_matrix, zstd_level, output);
        println!("Querying took {} seconds.", start_query.elapsed().as_secs());
    } else if let Some(line_file) = matches.value_of("query_by_line") {
        println!("\nONIKA QUERY by line: {}", line_file);
        let start_query = Instant::now();
        index.query_file_line_by_line(line_file, threshold, is_matrix, zstd_level, output);
        println!("Querying took {} seconds.", start_query.elapsed().as_secs());
    }
}
