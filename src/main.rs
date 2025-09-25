#[macro_use]
extern crate clap;
extern crate flate2;
extern crate rayon;
extern crate zstd;
extern crate f128;
extern crate num_traits;
extern crate needletail;
extern crate byteorder;
extern crate ahash;
extern crate tempfile;


mod onika_index;
mod utils;

use clap::{App, SubCommand, Arg, ArgGroup, ArgMatches};
use onika_index::{Index, IndexBuilder, SketchMode, open_compressed_file};
use std::time::Instant;
use std::io::{ BufRead};

/// Estimates the average document size by sampling the first 100 documents or sequences.
fn estimate_genome_size(matches: &ArgMatches) -> u32 {
    const SAMPLE_SIZE: usize = 100;
    let mut total_len: u64 = 0;
    let mut doc_count: u64 = 0;

    println!("Sampling up to {} documents/sequences to estimate average size...", SAMPLE_SIZE);
    
    let fof_path = matches.value_of("input_fof").or(matches.value_of("ref_fof"));
    let fasta_path = matches.value_of("input_fasta").or(matches.value_of("ref_fasta"));

    if let Some(fof_path) = fof_path {
        if let Ok(reader) = open_compressed_file(fof_path) {
            'outer: for filepath_res in reader.lines() {
                if let Ok(filepath) = filepath_res {
                    if let Ok(reader) = open_compressed_file(&filepath) {
                        if let Ok(mut file_reader) = needletail::parse_fastx_reader(reader) {
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
            }
        }
    } else if let Some(single_file_path) = fasta_path {
        if let Ok(reader) = open_compressed_file(single_file_path) {
            if let Ok(mut fastx_reader) = needletail::parse_fastx_reader(reader) {
                while let Some(Ok(record)) = fastx_reader.next() {
                    if doc_count >= SAMPLE_SIZE as u64 {
                        break;
                    }
                    total_len += record.seq().len() as u64;
                    doc_count += 1;
                }
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
        .about("A tool for MinHash sketching and all-versus-all comparison.")
        .subcommand(SubCommand::with_name("sketch")
            .about("Builds sketches from a dataset and saves them to a file.")
            .arg(Arg::with_name("input_fof").long("input-fof").value_name("FILE").help("Input dataset: a file of FASTA/Q file paths.").takes_value(true))
            .arg(Arg::with_name("input_fasta").long("input-fasta").value_name("FILE").help("Input dataset: a single FASTA/Q file where each record is a document.").takes_value(true))
            .group(ArgGroup::with_name("input_mode").args(&["input_fof", "input_fasta"]).required(true))
            .arg(Arg::with_name("output").short("o").long("output").value_name("FILE").help("Output file to save the sketch index (.bin).").required(true).takes_value(true))
            .arg(Arg::with_name("no_compress").long("no-compress").help("Skips sorting and compressing the index; faster build, larger index."))
            .arg(Arg::with_name("k").short("k").long("k_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("s").short("s").long("s_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("w").short("w").long("w_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("e").short("e").long("e_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("threads").long("threads").value_name("INT").takes_value(true))
            .arg(Arg::with_name("sketch_mode").long("sketch-mode").value_name("MODE").default_value("default").takes_value(true))
            .arg(Arg::with_name("zstd_level").long("zstd-level").value_name("LEVEL").default_value("1").takes_value(true))
        )
        .subcommand(SubCommand::with_name("compare")
            .about("Compares two datasets, which can be raw files or pre-computed sketches.")
            .arg(Arg::with_name("ref_fof").long("ref-fof").value_name("FILE").takes_value(true))
            .arg(Arg::with_name("ref_fasta").long("ref-fasta").value_name("FILE").takes_value(true))
            .arg(Arg::with_name("ref_sketch").long("ref-sketch").value_name("FILE").takes_value(true))
            .group(ArgGroup::with_name("reference").args(&["ref_fof", "ref_fasta", "ref_sketch"]).required(true))
            .arg(Arg::with_name("query_fof").long("query-fof").value_name("FILE").takes_value(true))
            .arg(Arg::with_name("query_fasta").long("query-fasta").value_name("FILE").takes_value(true))
            .arg(Arg::with_name("query_sketch").long("query-sketch").value_name("FILE").takes_value(true))
            .group(ArgGroup::with_name("query").args(&["query_fof", "query_fasta", "query_sketch"]).required(true))
            .arg(Arg::with_name("output").short("o").long("output").value_name("FILE").default_value("output.txt").takes_value(true))
            .arg(Arg::with_name("threshold").long("threshold").value_name("FLOAT").default_value("0.0").takes_value(true))
            .arg(Arg::with_name("matrix").long("matrix"))
            .arg(Arg::with_name("print_stats").long("print-stats"))
            .arg(Arg::with_name("no_compress").long("no-compress").help("For on-the-fly sketching, skips sorting and compression."))
            .arg(Arg::with_name("k").short("k").long("k_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("s").short("s").long("s_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("w").short("w").long("w_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("e").short("e").long("e_size").value_name("INT").takes_value(true))
            .arg(Arg::with_name("threads").long("threads").value_name("INT").takes_value(true))
            .arg(Arg::with_name("sketch_mode").long("sketch-mode").value_name("MODE").default_value("default").takes_value(true))
            .arg(Arg::with_name("zstd_level").long("zstd-level").value_name("LEVEL").default_value("1").takes_value(true))
            .arg(Arg::with_name("io_threads").long("io-threads").value_name("INT").default_value("16").takes_value(true))
            .arg(Arg::with_name("io_buffer").long("io-buffer").value_name("INT").default_value("4").takes_value(true))
        )
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        let w = value_t!(matches, "w", u32).unwrap_or(16);
        let s = value_t!(matches, "s", u32).unwrap_or(16);
        let k = value_t!(matches, "k", u32).unwrap_or(31);
        let threads = value_t!(matches, "threads", usize).unwrap_or(32);
        let sketch_mode = match matches.value_of("sketch_mode").unwrap() {
            "perfect" => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let zstd_level = value_t!(matches, "zstd_level", i32).unwrap_or(1);
        let output_file = matches.value_of("output").unwrap();
        let compress = !matches.is_present("no_compress");
        
        let e = if sketch_mode == SketchMode::Perfect && !matches.is_present("e") {
            let estimated_e = estimate_genome_size(matches);
            println!("Estimated genome size (e): {}", estimated_e);
            estimated_e
        } else {
            value_t!(matches, "e", u32).unwrap_or(5_000_000)
        };

        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
        
        println!("--- Building Sketch Index ---");
        let start_time = Instant::now();
        let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);

        if let Some(fof) = matches.value_of("input_fof") {
            builder.index_file_of_files(fof);
        } else if let Some(fasta_file) = matches.value_of("input_fasta") {
            builder.index_fasta_file(fasta_file);
        }
        
        let index = builder.into_final_index("tmp_onika",compress, zstd_level).expect("Failed to finalize index.");
        println!("Sketching took {} seconds.", start_time.elapsed().as_secs());
        
        println!("Dumping sketch index to file: {}", output_file);
        if let Err(e) = index.dump_index_disk(output_file) {
            eprintln!("Error dumping index: {}", e);
        }

    } else if let Some(matches) = matches.subcommand_matches("compare") {
        let threads = value_t!(matches, "threads", usize).unwrap_or(32);
        let sketch_mode = match matches.value_of("sketch_mode").unwrap() {
            "perfect" => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let zstd_level = value_t!(matches, "zstd_level", i32).unwrap_or(1);
        let output_file = matches.value_of("output").unwrap();
        let threshold = value_t!(matches, "threshold", f64).unwrap_or(0.0);
        let is_matrix = matches.is_present("matrix");
        let compress = !matches.is_present("no_compress");
        let io_threads = value_t!(matches, "io_threads", usize).unwrap_or(16);
        let io_buffer = value_t!(matches, "io_buffer", usize).unwrap_or(16);
        
        let w = value_t!(matches, "w", u32).unwrap_or(16);
        let s = value_t!(matches, "s", u32).unwrap_or(16);
        let k = value_t!(matches, "k", u32).unwrap_or(31);
        let e = if sketch_mode == SketchMode::Perfect && !matches.is_present("e") {
            let estimated_e = estimate_genome_size(matches);
            println!("Estimated genome size (e): {}", estimated_e);
            estimated_e
        } else {
            value_t!(matches, "e", u32).unwrap_or(5_000_000)
        };

        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
        
        let ref_index = if let Some(sketch_file) = matches.value_of("ref_sketch") {
            Index::from_file(sketch_file).expect("Failed to load reference sketch file")
        } else {
            println!("--- Building Reference Index On-the-fly ---");
            let start = Instant::now();
            let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);
            if let Some(fof) = matches.value_of("ref_fof") {
                builder.index_file_of_files(fof);
            } else if let Some(fasta_file) = matches.value_of("ref_fasta") {
                builder.index_fasta_file(fasta_file);
            }
            let index = builder.into_final_index("tmp_onika",compress, zstd_level).expect("Failed to finalize reference index.");
            println!("Reference indexing took {} seconds.", start.elapsed().as_secs());
            index
        };

        let query_index = if let Some(sketch_file) = matches.value_of("query_sketch") {
            Index::from_file(sketch_file).expect("Failed to load query sketch file")
        } else {
            println!("\n--- Building Query Index On-the-fly ---");
            let start = Instant::now();
            let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);
            if let Some(fof) = matches.value_of("query_fof") {
                builder.index_file_of_files(fof);
            } else if let Some(fasta_file) = matches.value_of("query_fasta") {
                builder.index_fasta_file(fasta_file);
            }
            let index = builder.into_final_index("tmp_onika",compress, zstd_level).expect("Failed to finalize query index.");
            println!("Query indexing took {} seconds.", start.elapsed().as_secs());
            index
        };
        
        if matches.is_present("print_stats") {
            println!("\n--- Reference Index Stats ---");
            ref_index.print_stats();
            println!("\n--- Query Index Stats ---");
            query_index.print_stats();
        }
        
        ref_index.all_vs_all_comparison(&query_index, threshold, is_matrix, zstd_level, output_file, io_threads, io_buffer,64);
    }
}

