use crate::utils::*;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc, Mutex,
};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use f128::f128;
use flate2::read::GzDecoder;
use needletail::parse_fastx_file;
use nthash::NtHashIterator;
use num_traits::{Float, ToPrimitive};
use rayon::prelude::*;
use rlimit;
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use stream_vbyte::{decode::decode, scalar::Scalar};
use zstd::stream::write::Encoder as ZstdEncoder;
use io::BufWriter;

type Gid = u32;

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(u8)]
pub enum SketchMode {
    Default = 0,
    Perfect = 1,
}

/// Input for the compute_sketch function, allowing either a single sequence or a file path.
pub enum SketchInput<'a> {
    Sequence(&'a [u8]),
    FilePath(&'a str),
}

/// The final, immutable, lock-free index used for querying.
/// GID lists can be stored in a compressed or raw format.
pub struct Index {
    ls: u32,
    w: u32,
    k: u32,
    e: u32,
    s: u64,
    sketch_size: u64,
    fingerprint_range: u64,
    mi: u64,
    b: u32, // Sub-division factor
    // Stores (original_gid_count, data_bytes)
    // data_bytes can be compressed GIDs or raw little-endian u32 GIDs
    data: Vec<(u32, Vec<u8>)>,
    genome_numbers: u64,
    sketch_mode: SketchMode,
    is_compressed: bool,
}

/// A temporary, mutable struct for building the index in memory.
pub struct IndexBuilder {
    ls: u32,
    w: u32,
    k: u32,
    e: u32,
    s: u64,
    sketch_size: u64,
    fingerprint_range: u64,
    mi: u64,
    b: u32, // Sub-division factor
    b_mask: u32,
    // Direct, thread-safe version of the final index structure.
    intermediate_data: Vec<Mutex<Vec<u32>>>,
    genome_numbers: AtomicU64,
    empty_sketches_count: AtomicU64,
    sketch_mode: SketchMode,
}

/// Opens a file, transparently handling .gz and .zst compression.
pub(crate) fn open_compressed_file(path: &str) -> io::Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    match Path::new(path).extension().and_then(OsStr::to_str) {
        Some("gz") => {
            let decoder = GzDecoder::new(file);
            Ok(Box::new(BufReader::new(decoder)))
        }
        Some("zst") => {
            let decoder = zstd::stream::read::Decoder::new(file)?;
            Ok(Box::new(BufReader::new(decoder)))
        }
        _ => Ok(Box::new(BufReader::new(file))),
    }
}

impl IndexBuilder {
    /// Creates a new IndexBuilder for in-memory construction.
    pub fn new(
        ils: u32,
        ik: u32,
        iw: u32,
        ie: u32,
        mode: SketchMode,
        sub_div_factor: u32,
    ) -> Self {
        let s = 1u64 << ils;
        let fingerprint_range = 1u64 << iw;
        
        let total_size = ((fingerprint_range * s) >> sub_div_factor) as usize;
        
        println!("Allocating in-memory index structure...");
        let intermediate_data = (0..total_size)
            .map(|_| Mutex::new(Vec::new()))
            .collect();

        Self {
            ls: ils,
            w: iw,
            k: ik,
            e: ie,
            s,
            sketch_size: s,
            fingerprint_range,
            mi: u64::MAX,
            b: sub_div_factor,
            b_mask: (1u32 << sub_div_factor) - 1,
            intermediate_data,
            genome_numbers: AtomicU64::new(0),
            empty_sketches_count: AtomicU64::new(0),
            sketch_mode: mode,
        }
    }

    /// Finalizes the index construction.
    /// # Arguments
    /// * `compress` - If true, GID lists will be sorted, delta-encoded, and compressed.
    pub fn into_final_index(self, compress: bool) -> Index {
        let data: Vec<(u32, Vec<u8>)> = if compress {
            println!("Sorting, delta-encoding, and compressing GID lists...");
            self.intermediate_data
                .into_par_iter()
                .map(|mutex| {
                    let mut gids = mutex.into_inner().unwrap();
                    if gids.is_empty() { return (0, Vec::new()); }
                    
                    gids.sort_unstable();
                    let original_count = gids.len() as u32;

                    let mut prev = gids[0];
                    for i in 1..gids.len() {
                        let current = gids[i];
                        gids[i] = current.wrapping_sub(prev);
                        prev = current;
                    }

                    let mut encoded_data = vec![0; 5 * gids.len()];
                    let encoded_len =
                        stream_vbyte::encode::encode::<Scalar>(&gids, &mut encoded_data);
                    encoded_data.truncate(encoded_len);
                    (original_count, encoded_data)
                })
                .collect()
        } else {
            println!("Storing raw, unsorted GID lists...");
            self.intermediate_data
                .into_par_iter()
                .map(|mutex| {
                    let gids = mutex.into_inner().unwrap();
                    if gids.is_empty() { return (0, Vec::new()); }

                    let original_count = gids.len() as u32;
                    let mut gids_as_bytes = Vec::with_capacity(gids.len() * 4);
                    for gid in gids {
                        gids_as_bytes.write_u32::<LittleEndian>(gid).unwrap();
                    }
                    (original_count, gids_as_bytes)
                })
                .collect()
        };

        let empty_count = self.empty_sketches_count.load(Ordering::Relaxed);
        if empty_count > 0 {
            eprintln!("\nWarning: {} input documents resulted in empty sketches and were not indexed.", empty_count);
        }
        println!("Index construction complete.");

        Index {
            ls: self.ls,
            w: self.w,
            k: self.k,
            e: self.e,
            s: self.s,
            sketch_size: self.sketch_size,
            fingerprint_range: self.fingerprint_range,
            mi: self.mi,
            b: self.b,
            data,
            genome_numbers: self.genome_numbers.load(Ordering::Relaxed),
            sketch_mode: self.sketch_mode,
            is_compressed: compress,
        }
    }

    pub fn index_file_line_by_line(&self, filestr: &str) {
        let reader = open_compressed_file(filestr)
            .expect("Could not open file with sequences line by line");
        let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();
        self.genome_numbers
            .store(lines.len() as u64, Ordering::SeqCst);

        lines.into_par_iter().enumerate().for_each(|(i, line)| {
            let seq_id = i as u32;
            let mut sketch = vec![self.mi; self.s as usize];
            self.compute_sketch(SketchInput::Sequence(line.as_bytes()), &mut sketch);

            if sketch.iter().all(|&x| x == self.mi) {
                self.empty_sketches_count.fetch_add(1, Ordering::Relaxed);
                return;
            }

            for (i, &fp) in sketch.iter().enumerate() {
                if fp < self.fingerprint_range {
                    let flat_index = (fp * self.sketch_size) + i as u64;
                    let new_index = (flat_index >> self.b) as usize;
                    let offset = (flat_index as u32) & self.b_mask;

                    if seq_id < (1 << (32 - self.b)) {
                        let encoded_gid = (seq_id << self.b) | offset;
                        self.intermediate_data[new_index].lock().unwrap().push(encoded_gid);
                    }
                }
            }
        });
    }

    pub fn index_file_of_files(&self, fof_path: &str) {
        let reader =
            open_compressed_file(fof_path).expect("Could not open file of files");
        let files_to_process: Vec<String> = reader.lines().filter_map(Result::ok).collect();
        self.genome_numbers
            .store(files_to_process.len() as u64, Ordering::SeqCst);

        files_to_process
            .into_par_iter()
            .enumerate()
            .for_each(|(i, filepath)| {
                let seq_id = i as u32;
                let mut sketch = vec![self.mi; self.s as usize];
                self.compute_sketch(SketchInput::FilePath(&filepath), &mut sketch);
                
                if sketch.iter().all(|&x| x == self.mi) {
                    self.empty_sketches_count.fetch_add(1, Ordering::Relaxed);
                    return;
                }

                for (i, &fp) in sketch.iter().enumerate() {
                    if fp < self.fingerprint_range {
                        let flat_index = (fp * self.sketch_size) + i as u64;
                        let new_index = (flat_index >> self.b) as usize;
                        let offset = (flat_index as u32) & self.b_mask;
    
                        if seq_id < (1 << (32 - self.b)) {
                            let encoded_gid = (seq_id << self.b) | offset;
                            self.intermediate_data[new_index].lock().unwrap().push(encoded_gid);
                        }
                    }
                }
            });
    }

    fn compute_sketch(&self, input: SketchInput, sketch: &mut Vec<u64>) {
        sketch.fill(self.mi);
        let mut empty_cell = self.s;
        let is_valid_dna =
            |c: &u8| matches!(*c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');

        let mut update_sketch = |sequence: &[u8]| {
            for subseq in sequence.split(|c| !is_valid_dna(c)) {
                if subseq.len() >= self.k as usize {
                    let h_iter = NtHashIterator::new(subseq, self.k as usize).unwrap();
                    for canon_hash in h_iter {
                        let fp = revhash64(canon_hash);
                        let bucket_id = (unrevhash64(canon_hash) >> (64 - self.ls)) as usize;
                        if bucket_id < sketch.len() {
                            if sketch[bucket_id] == self.mi {
                                empty_cell -= 1;
                                sketch[bucket_id] = fp;
                            } else if sketch[bucket_id] > fp {
                                sketch[bucket_id] = fp;
                            }
                        }
                    }
                }
            }
        };

        match input {
            SketchInput::Sequence(seq) => update_sketch(seq),
            SketchInput::FilePath(path) => {
                match open_compressed_file(path) {
                    Ok(reader) => {
                        if let Ok(mut file_reader) = needletail::parse_fastx_reader(reader) {
                            while let Some(Ok(record)) = file_reader.next() {
                                update_sketch(&record.seq());
                            }
                        }
                    }
                    Err(e) => eprintln!("Warning: could not open file {}: {}", path, e),
                }
            }
        }

        if empty_cell < self.s {
            self.sketch_densification(sketch, empty_cell as u32);
            for val in sketch.iter_mut() {
                if *val != self.mi {
                    *val = match self.sketch_mode {
                        SketchMode::Default => *val & (self.fingerprint_range - 1),
                        SketchMode::Perfect => self.get_perfect_fingerprint(*val),
                    };
                }
            }
        }
    }

    fn sketch_densification(&self, sketch: &mut Vec<u64>, mut empty_cell: u32) {
        let mut step = 0;
        let size = sketch.len();
        while empty_cell != 0 {
            let mut changes = Vec::new();
            for i in 0..size {
                if sketch[i] != self.mi {
                    let hash_pos = (hash_family(sketch[i], step) % size as u64) as usize;
                    if hash_pos < sketch.len() && sketch[hash_pos] == self.mi {
                        changes.push((hash_pos, sketch[i]));
                    }
                }
            }
            if changes.is_empty() && empty_cell > 0 {
                break;
            }

            for (pos, val) in changes {
                if sketch[pos] == self.mi {
                    sketch[pos] = val;
                    empty_cell -= 1;
                    if empty_cell == 0 {
                        return;
                    }
                }
            }
            step += 1;
        }
    }

    fn get_perfect_fingerprint(&self, hashed: u64) -> u64 {
        if hashed == self.mi {
            return self.mi;
        }

        let b = f128::from(hashed >> 32);
        let twopowern = f128::from(1u128 << 32);
        let ratio = f128::from(self.e) / f128::from(self.s);

        let top = (twopowern - b).powf(ratio);
        let bottom = twopowern.powf(ratio);

        let scaled_top = f128::from(self.fingerprint_range) * top;
        let result_f = f128::from(self.fingerprint_range) - (scaled_top / bottom);

        if result_f < f128::from(0.0) {
            0
        } else {
            result_f.floor().to_u64().unwrap_or(0)
        }
    }
}

impl Index {
    pub fn from_file(filestr: &str) -> io::Result<Self> {
        println!("Loading index from disk: {}...", filestr);
        let file = File::open(filestr)?;
        let mut zstd_reader = zstd::stream::read::Decoder::new(file)?;

        let ls = zstd_reader.read_u32::<LittleEndian>()?;
        let k = zstd_reader.read_u32::<LittleEndian>()?;
        let w = zstd_reader.read_u32::<LittleEndian>()?;
        let e = zstd_reader.read_u32::<LittleEndian>()?;
        let fingerprint_range = zstd_reader.read_u64::<LittleEndian>()?;
        let genome_numbers = zstd_reader.read_u64::<LittleEndian>()?;
        let sketch_mode = match zstd_reader.read_u8()? {
            1 => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let b = zstd_reader.read_u32::<LittleEndian>()?;
        let is_compressed = zstd_reader.read_u8()? != 0;

        let s = 1u64 << ls;
        let total_size = ((fingerprint_range * s) >> b) as usize;
        let mut data = Vec::with_capacity(total_size);

        for _ in 0..total_size {
            let original_count = zstd_reader.read_u32::<LittleEndian>()?;
            let gids_len = zstd_reader.read_u32::<LittleEndian>()?;
            let mut gids_bytes = vec![0u8; gids_len as usize];
            zstd_reader.read_exact(&mut gids_bytes)?;
            data.push((original_count, gids_bytes));
        }

        println!("Finished loading index.");
        Ok(Self {
            ls,
            w,
            k,
            e,
            s,
            sketch_size: s,
            fingerprint_range,
            mi: u64::MAX,
            b,
            data,
            genome_numbers,
            sketch_mode,
            is_compressed,
        })
    }

    pub fn dump_index_disk(&self, filestr: &str, zstd_level: i32) -> io::Result<()> {
        if !self.is_compressed {
            eprintln!("Warning: Dumping an uncompressed index is not recommended for portability.");
        }
        let file = File::create(filestr)?;
        let mut writer = ZstdEncoder::new(file, zstd_level)?.auto_finish();

        writer.write_u32::<LittleEndian>(self.ls)?;
        writer.write_u32::<LittleEndian>(self.k)?;
        writer.write_u32::<LittleEndian>(self.w)?;
        writer.write_u32::<LittleEndian>(self.e)?;
        writer.write_u64::<LittleEndian>(self.fingerprint_range)?;
        writer.write_u64::<LittleEndian>(self.genome_numbers)?;
        writer.write_u8(self.sketch_mode as u8)?;
        writer.write_u32::<LittleEndian>(self.b)?;
        writer.write_u8(self.is_compressed as u8)?;

        for (original_count, gids_bytes) in &self.data {
            writer.write_u32::<LittleEndian>(*original_count)?;
            writer.write_u32::<LittleEndian>(gids_bytes.len() as u32)?;
            writer.write_all(gids_bytes)?;
        }

        Ok(())
    }

    pub fn print_stats(&self) {
        println!("\nINDEX STATISTICS");
        println!("=========================================================");

        let sizes: Vec<u64> = self
            .data
            .par_iter()
            .map(|(count, _)| *count as u64)
            .filter(|&size| size > 0)
            .collect();

        if sizes.is_empty() {
            println!("Index is empty. No stats to compute.");
            return;
        }

        let num_lists = sizes.len() as f64;
        let sum: u64 = sizes.iter().sum();
        let mean = sum as f64 / num_lists;

        let min_size = *sizes.iter().min().unwrap_or(&0);
        let max_size = *sizes.iter().max().unwrap_or(&0);

        println!("Number of non-empty lists: {}", sizes.len());
        println!("Total elements in index: {}", sum);
        println!("Min list size: {}", min_size);
        println!("Max list size: {}", max_size);
        println!("Mean list size: {:.4}", mean);
        println!("Index is compressed: {}", self.is_compressed);
    }
      /// Performs an all-versus-all comparison between this index (reference) and another (query).
    pub fn all_vs_all_comparison(
        &self,
        query_index: &Index,
        threshold: f64,
        is_matrix: bool,
        zstd_level: i32,
        output_file: &str,
    ) {
        println!("\n--- Starting All-Versus-All Comparison ---");
        let start_time = std::time::Instant::now();

        let num_queries = query_index.genome_numbers as usize;
        let num_refs = self.genome_numbers as usize;

        if self.data.len() != query_index.data.len() {
            panic!("Incompatible indices for comparison: different total sizes.");
        }
        
        // OPTIMIZATION: Use a sharded HashMap to avoid the reduce bottleneck.
        // The number of shards is equal to the number of threads Rayon will use.
        let num_shards = rayon::current_num_threads();
        let scores_shards: Vec<Mutex<HashMap<(u32, u32), u32>>> =
            (0..num_shards).map(|_| Mutex::new(HashMap::new())).collect();

        (0..self.data.len()).into_par_iter().for_each(|flat_index| {
            let (ref_count, ref_bytes) = &self.data[flat_index];
            if *ref_count > 0 {
                let (query_count, query_bytes) = &query_index.data[flat_index];
                if *query_count > 0 {
                    let ref_gids_encoded =
                        decode_gid_list(*ref_count, ref_bytes, self.is_compressed);
                    let query_gids_encoded =
                        decode_gid_list(*query_count, query_bytes, query_index.is_compressed);

                    let b_mask = (1u32 << self.b) - 1;
                    let b_shift = self.b;

                    let mut ref_pairs: Vec<(u32, u32)> = ref_gids_encoded
                        .iter()
                        .map(|&gid| (gid & b_mask, gid >> b_shift))
                        .collect();
                    ref_pairs.sort_unstable();

                    let mut query_pairs: Vec<(u32, u32)> = query_gids_encoded
                        .iter()
                        .map(|&gid| (gid & b_mask, gid >> b_shift))
                        .collect();
                    query_pairs.sort_unstable();

                    let mut i = 0;
                    let mut j = 0;
                    while i < ref_pairs.len() && j < query_pairs.len() {
                        let ref_offset = ref_pairs[i].0;
                        let query_offset = query_pairs[j].0;

                        if ref_offset < query_offset {
                            i += 1;
                        } else if query_offset < ref_offset {
                            j += 1;
                        } else {
                            let mut i_end = i;
                            while i_end < ref_pairs.len() && ref_pairs[i_end].0 == ref_offset {
                                i_end += 1;
                            }

                            let mut j_end = j;
                            while j_end < query_pairs.len() && query_pairs[j_end].0 == query_offset {
                                j_end += 1;
                            }

                            for query_idx in j..j_end {
                                for ref_idx in i..i_end {
                                    let query_gid = query_pairs[query_idx].1;
                                    let ref_gid = ref_pairs[ref_idx].1;
                                    
                                    // Hash the pair to determine which shard to lock
                                    let key = (query_gid, ref_gid);
                                    let shard_index = key.0 as usize % num_shards;
                                    let mut shard = scores_shards[shard_index].lock().unwrap();
                                    *shard.entry(key).or_insert(0) += 1;
                                }
                            }

                            i = i_end;
                            j = j_end;
                        }
                    }
                }
            }
        });
        
        println!(
            "Score calculation finished in {} seconds. Aggregating and formatting results...",
            start_time.elapsed().as_secs()
        );

        let mut scores_by_query: Vec<Vec<(Gid, u32)>> = vec![Vec::new(); num_queries];
        for shard_mutex in scores_shards {
            let shard_map = shard_mutex.into_inner().unwrap();
            for ((query_gid, ref_gid), count) in shard_map {
                if (query_gid as usize) < num_queries {
                    scores_by_query[query_gid as usize].push((ref_gid, count));
                }
            }
        }
        
        println!("Parallel formatting of results...");
        let output_lines: Vec<String> = (0..num_queries)
            .into_par_iter()
            .map(|query_gid| {
                let mut line = String::with_capacity(256); // Pre-allocate string
                write!(&mut line, "query_{}\t", query_gid).unwrap();
                
                let res = &scores_by_query[query_gid];

                if is_matrix {
                    let mut full_row = vec![0.0; num_refs];
                    for &(ref_gid, count) in res {
                        if (ref_gid as usize) < num_refs {
                            full_row[ref_gid as usize] = count as f64 / self.s as f64;
                        }
                    }
                    let similarities: Vec<String> = full_row
                        .iter()
                        .map(|&similarity| {
                            if similarity >= threshold {
                                format!("{:.4}", similarity)
                            } else {
                                "0.0000".to_string()
                            }
                        })
                        .collect();
                    line.push_str(&similarities.join(","));
                } else { // Sparse format
                    let mut hits: Vec<(Gid, f64)> = res
                        .iter()
                        .map(|&(ref_gid, count)| (ref_gid, count as f64 / self.s as f64))
                        .filter(|&(_, sim)| sim >= threshold)
                        .collect();

                    hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

                    line.push_str(
                        &hits
                            .iter()
                            .map(|(gid, sim)| format!("{}:{:.4}", gid, sim))
                            .collect::<Vec<String>>()
                            .join(","),
                    );
                }
                line
            })
            .collect();
        
        println!("Writing results to disk...");
        let mut writer = self.create_writer(output_file, zstd_level);
        for line in output_lines {
            writeln!(writer, "{}", line).unwrap();
        }

        println!(
            "All-vs-all comparison complete. Total time: {} seconds.",
            start_time.elapsed().as_secs()
        );
    }


    fn create_writer<'a>(&self, output_file: &'a str, zstd_level: i32) -> Box<dyn Write + 'a> {
        let file = File::create(output_file).expect("Cannot create output file");
        if zstd_level > 0 {
            Box::new(ZstdEncoder::new(file, zstd_level).unwrap().auto_finish())
        } else {
            Box::new(BufWriter::new(file))
        }
    }
}


/// Helper to decode a GID list from bytes depending on whether it's compressed.
fn decode_gid_list(count: u32, bytes: &[u8], is_compressed: bool) -> Vec<u32> {
    if count == 0 {
        return Vec::new();
    }

    if is_compressed {
        let mut gids = vec![0u32; count as usize]; // Pre-allocate with length for the decoder
        decode::<Scalar>(bytes, count as usize, &mut gids);
        for j in 1..gids.len() {
            gids[j] = gids[j].wrapping_add(gids[j - 1]);
        }
        gids
    } else {
        let mut gids = Vec::with_capacity(count as usize); // Allocate capacity and push
        let mut cursor = io::Cursor::new(bytes);
        for _ in 0..count {
            gids.push(cursor.read_u32::<LittleEndian>().unwrap());
        }
        gids
    }
}

