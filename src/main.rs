extern crate byteorder;
extern crate clap;
extern crate f128;
extern crate flate2;
extern crate needletail;
extern crate num_cpus;
extern crate num_traits;
extern crate rayon;
extern crate tempfile;
extern crate zstd;
use jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

mod onika_index;
mod utils;

use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command};
use onika_index::{open_compressed_file, Index, IndexBuilder, SketchMode};
use std::io::BufRead;
use std::time::Instant;

/// Estimates the average document size by sampling the first 100 documents or sequences.
fn estimate_genome_size(matches: &ArgMatches) -> u32 {
    const SAMPLE_SIZE: usize = 100;
    let mut total_len: u64 = 0;
    let mut doc_count: u64 = 0;

    println!(
        "Sampling up to {} documents/sequences to estimate average size...",
        SAMPLE_SIZE
    );

    let fof_path: Option<&str> = matches
        .get_one::<String>("input_fof")
        .map(|s| s.as_str())
        .or_else(|| matches.get_one::<String>("ref_fof").map(|s| s.as_str()));
    let fasta_path: Option<&str> = matches
        .get_one::<String>("input_fasta")
        .map(|s| s.as_str())
        .or_else(|| matches.get_one::<String>("ref_fasta").map(|s| s.as_str()));

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
    let cli_matches = Command::new("Onika-rs")
        .version("1.0")
        .author("Rust Translation")
        .about("A tool for MinHash sketching and all-versus-all comparison.")
        .arg(
            Arg::new("threads")
                .long("threads")
                .global(true)
                .help("Set the number of threads for parallel operations.")
                .value_name("INT")
                .value_parser(clap::value_parser!(usize)),
        )
        .subcommand(
            Command::new("sketch")
                .about("Builds sketches from a dataset and saves them to a file.")
                .arg(
                    Arg::new("input_fof")
                        .long("input-fof")
                        .value_name("FILE")
                        .help("Input dataset: a file of FASTA/Q file paths."),
                )
                .arg(
                    Arg::new("input_fasta")
                        .long("input-fasta")
                        .value_name("FILE")
                        .help("Input dataset: a single FASTA/Q file where each record is a document."),
                )
                .group(
                    ArgGroup::new("input_mode")
                        .args(["input_fof", "input_fasta"])
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_name("FILE")
                        .help("Output file to save the sketch index (.bin).")
                        .required(true),
                )
                .arg(
                    Arg::new("no_compress")
                        .long("no-compress")
                        .help("Skips stream-vbyte compression (still delta encodes); faster build, larger index.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("reorder_similarity")
                        .long("reorder-similarity")
                        .help("After building the sketch, greedily reorder sketches by self-similarity to improve locality.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("reorder_threshold")
                        .long("reorder-threshold")
                        .value_name("FLOAT")
                        .default_value("0.5")
                        .help("Minimum estimated similarity (0-1) a pair must reach to be considered during greedy reordering.")
                        .value_parser(clap::value_parser!(f64)),
                )
                .arg(
                    Arg::new("reorder_sample")
                        .long("reorder-sample")
                        .value_name("INT")
                        .default_value("100")
                        .help("Use only the first N partitions (sketch positions) when estimating similarity during reordering (0 = full sketch).")
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    Arg::new("k")
                        .short('k')
                        .long("k_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("s")
                        .short('s')
                        .long("s_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("w")
                        .short('w')
                        .long("w_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("e")
                        .short('e')
                        .long("e_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("sketch_mode")
                        .long("sketch-mode")
                        .value_name("MODE")
                        .default_value("default")
                        .value_parser(["default", "perfect"]),
                )
                .arg(
                    Arg::new("zstd_level")
                        .long("zstd-level")
                        .value_name("LEVEL")
                        .default_value("1")
                        .value_parser(clap::value_parser!(i32)),
                ),
        )
        .subcommand(
            Command::new("compare")
                .about("Compares two datasets, which can be raw files or pre-computed sketches.")
                .arg(
                    Arg::new("ref_fof")
                        .long("ref-fof")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("ref_fasta")
                        .long("ref-fasta")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("ref_sketch")
                        .long("ref-sketch")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("query_fof")
                        .long("query-fof")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("query_fasta")
                        .long("query-fasta")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("query_sketch")
                        .long("query-sketch")
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_name("FILE")
                        .required(true),
                )
                .arg(
                    Arg::new("threshold")
                        .long("threshold")
                        .value_name("FLOAT")
                        .value_parser(clap::value_parser!(f64)),
                )
                .arg(
                    Arg::new("matrix")
                        .long("matrix")
                        .help("Output a full similarity matrix instead of sparse format.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("probabilistic_threshold")
                        .long("prob-threshold-heuristic")
                        .help("Enable probability-based pruning when estimating pairs unlikely to reach the threshold.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("prob_threshold_probability")
                        .long("prob-threshold-probability")
                        .value_name("FLOAT")
                        .default_value("0.001")
                        .help("Tail probability used by the probabilistic pruning heuristic (e.g. 0.001 for 1 in 1000).")
                        .value_parser(clap::value_parser!(f64)),
                )
                .arg(
                    Arg::new("print_stats")
                        .long("print-stats")
                        .help("Print basic statistics for both indices before comparison.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("no_compress")
                        .long("no-compress")
                        .help("Write uncompressed shard files. Useful for debugging or if compression overhead is too high.")
                        .action(ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("k")
                        .short('k')
                        .long("k_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("s")
                        .short('s')
                        .long("s_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("w")
                        .short('w')
                        .long("w_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("e")
                        .short('e')
                        .long("e_size")
                        .value_name("INT")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("sketch_mode")
                        .long("sketch-mode")
                        .value_name("MODE")
                        .default_value("default")
                        .value_parser(["default", "perfect"]),
                )
                .arg(
                    Arg::new("zstd_level")
                        .long("zstd-level")
                        .value_name("LEVEL")
                        .default_value("1")
                        .value_parser(clap::value_parser!(i32)),
                )
                .arg(
                    Arg::new("shard_zstd_level")
                        .long("shard-zstd-level")
                        .value_name("LEVEL")
                        .default_value("1")
                        .help("Zstd compression level for each shard file. Higher values yield better compression but slower writes. Use 0 to disable.")
                        .value_parser(clap::value_parser!(i32)),
                )
                .arg(
                    Arg::new("temp_dir")
                        .long("temp-dir")
                        .value_name("PATH")
                        .help("Directory for temporary files. Defaults to the system temp directory."),
                ),
        )
        .get_matches();

    let threads = cli_matches
        .get_one::<usize>("threads")
        .copied()
        .unwrap_or_else(num_cpus::get);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    if let Some(matches) = cli_matches.subcommand_matches("sketch") {
        let w = matches.get_one::<u32>("w").copied().unwrap_or(16);
        let s = matches.get_one::<u32>("s").copied().unwrap_or(10);
        let k = matches.get_one::<u32>("k").copied().unwrap_or(31);
        let sketch_mode = match matches
            .get_one::<String>("sketch_mode")
            .map(|s| s.as_str())
            .unwrap_or("default")
        {
            "perfect" => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let zstd_level = matches.get_one::<i32>("zstd_level").copied().unwrap_or(1);
        let output_file = matches.get_one::<String>("output").unwrap();
        let compress = !matches.get_flag("no_compress");

        let e_arg = matches.get_one::<u32>("e").copied();
        let e = if sketch_mode == SketchMode::Perfect && e_arg.is_none() {
            let estimated_e = estimate_genome_size(matches);
            println!("Estimated genome size (e): {}", estimated_e);
            estimated_e
        } else {
            e_arg.unwrap_or(5_000_000)
        };

        println!("--- Building Sketch Index ---");
        let start_time = Instant::now();
        let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);

        if let Some(fof) = matches.get_one::<String>("input_fof") {
            builder.index_file_of_files(fof);
        } else if let Some(fasta_file) = matches.get_one::<String>("input_fasta") {
            builder.index_fasta_file(fasta_file);
        }

        if matches.get_flag("reorder_similarity") {
            let reorder_threshold = matches
                .get_one::<f64>("reorder_threshold")
                .copied()
                .unwrap_or(0.5)
                .clamp(0.0, 1.0);
            let reorder_sample = matches
                .get_one::<usize>("reorder_sample")
                .copied()
                .unwrap_or(100);
            println!(
                "Reordering sketches by greedy self-similarity chain (min sim {:.3}, sample {})...",
                reorder_threshold,
                if reorder_sample == 0 {
                    "all".to_string()
                } else {
                    reorder_sample.to_string()
                }
            );
            let reorder_start = Instant::now();
            match builder.reorder_by_similarity(reorder_threshold, reorder_sample, threads) {
                Ok(order) => {
                    let preview_len = order.len().min(8);
                    if preview_len > 0 {
                        let preview: Vec<String> = order
                            .iter()
                            .take(preview_len)
                            .map(|idx| idx.to_string())
                            .collect();
                        println!(
                            "Reorder complete in {} seconds. Preview (new -> original ids): [{}{}]",
                            reorder_start.elapsed().as_secs(),
                            preview.join(", "),
                            if order.len() > preview_len {
                                ", ..."
                            } else {
                                ""
                            }
                        );
                    } else {
                        println!(
                            "Reorder complete in {} seconds. No sketches to reorder.",
                            reorder_start.elapsed().as_secs()
                        );
                    }
                }
                Err(err) => {
                    eprintln!(
                        "Similarity-based reorder failed after {} seconds: {}",
                        reorder_start.elapsed().as_secs(),
                        err
                    );
                }
            }
        }

        let index = builder
            .into_final_index(output_file, compress, zstd_level)
            .expect("Failed to finalize index.");
        println!("Sketching took {} seconds.", start_time.elapsed().as_secs());

        println!("Dumping sketch index to file: {}", output_file);
        if let Err(e) = index.dump_index_disk(output_file) {
            eprintln!("Error dumping index: {}", e);
        }
    } else if let Some(matches) = cli_matches.subcommand_matches("compare") {
        let sketch_mode = match matches
            .get_one::<String>("sketch_mode")
            .map(|s| s.as_str())
            .unwrap_or("default")
        {
            "perfect" => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let zstd_level = matches.get_one::<i32>("zstd_level").copied().unwrap_or(1);
        let shard_zstd_level = matches
            .get_one::<i32>("shard_zstd_level")
            .copied()
            .unwrap_or(2);
        let temp_dir_path = matches.get_one::<String>("temp_dir").map(|s| s.as_str());
        let output_file = matches.get_one::<String>("output").unwrap();
        let threshold = matches.get_one::<f64>("threshold").copied().unwrap_or(0.0);
        let is_matrix = matches.get_flag("matrix");
        let compress = !matches.get_flag("no_compress");
        let use_probabilistic_pruning = matches.get_flag("probabilistic_threshold");
        let prob_threshold_probability = matches
            .get_one::<f64>("prob_threshold_probability")
            .copied()
            .unwrap_or(1.0 / 10_000.0);

        let w = matches.get_one::<u32>("w").copied().unwrap_or(16);
        let s = matches.get_one::<u32>("s").copied().unwrap_or(16);
        let k = matches.get_one::<u32>("k").copied().unwrap_or(31);
        let e_arg = matches.get_one::<u32>("e").copied();
        let e = if sketch_mode == SketchMode::Perfect && e_arg.is_none() {
            let estimated_e = estimate_genome_size(matches);
            println!("Estimated genome size (e): {}", estimated_e);
            estimated_e
        } else {
            e_arg.unwrap_or(5_000_000)
        };

        let ref_index = if let Some(sketch_file) = matches.get_one::<String>("ref_sketch") {
            Index::from_file(sketch_file).expect("Failed to load reference sketch file")
        } else {
            println!("--- Building Reference Index On-the-fly ---");
            let start = Instant::now();
            let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);
            if let Some(fof) = matches.get_one::<String>("ref_fof") {
                builder.index_file_of_files(fof);
            } else if let Some(fasta_file) = matches.get_one::<String>("ref_fasta") {
                builder.index_fasta_file(fasta_file);
            }
            let index = builder
                .into_final_index(output_file, compress, zstd_level)
                .expect("Failed to finalize reference index.");
            println!(
                "Reference indexing took {} seconds.",
                start.elapsed().as_secs()
            );
            index
        };

        let query_index = if let Some(sketch_file) = matches.get_one::<String>("query_sketch") {
            Index::from_file(sketch_file).expect("Failed to load query sketch file")
        } else {
            println!("\n--- Building Query Index On-the-fly ---");
            let start = Instant::now();
            let mut builder = IndexBuilder::new(s, k, w, e, sketch_mode);
            if let Some(fof) = matches.get_one::<String>("query_fof") {
                builder.index_file_of_files(fof);
            } else if let Some(fasta_file) = matches.get_one::<String>("query_fasta") {
                builder.index_fasta_file(fasta_file);
            }
            let index = builder
                .into_final_index(output_file, compress, zstd_level)
                .expect("Failed to finalize query index.");
            println!("Query indexing took {} seconds.", start.elapsed().as_secs());
            index
        };

        if matches.get_flag("print_stats") {
            println!("\n--- Reference Index Stats ---");
            ref_index.print_stats();
            println!("\n--- Query Index Stats ---");
            query_index.print_stats();
        }

        ref_index.all_vs_all_comparison(
            &query_index,
            threshold,
            is_matrix,
            zstd_level,
            output_file,
            threads,
            shard_zstd_level,
            temp_dir_path,
            use_probabilistic_pruning,
            prob_threshold_probability,
        );
    }
}
