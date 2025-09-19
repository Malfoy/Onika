use crate::utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Read};
use std::sync::{atomic::{AtomicU64, Ordering}, Arc, Mutex};
use rayon::prelude::*;
use nthash::NtHashIterator;
use needletail::{parse_fastx_file};
use zstd::stream::write::Encoder as ZstdEncoder;
use f128::f128;
use num_traits::{Float, ToPrimitive};
use byteorder::{ReadBytesExt, WriteBytesExt, LittleEndian};
use tempfile::TempDir;
use rlimit;
use stream_vbyte::{decode::decode, scalar::Scalar};
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;


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
/// GID lists are stored in a compressed format using stream-vbyte.
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
    // Stores (original_gid_count, compressed_delta_encoded_gids)
    data: Vec<(u32, Vec<u8>)>, 
    genome_numbers: u64,
    sketch_mode: SketchMode,
}

/// A temporary, mutable struct for building the index in parallel using temporary files.
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
    temp_dir: TempDir,
    writers: Vec<Arc<Mutex<BufWriter<File>>>>,
    genome_numbers: AtomicU64,
    sketch_mode: SketchMode,
}

impl IndexBuilder {
    /// Creates a new IndexBuilder.
    ///
    /// # Arguments
    /// * `num_temp_files` - The number of temporary files to use for distributing write load.
    ///   A value around 256 or 512 is a reasonable start.
    pub fn new(ils: u32, ik: u32, iw: u32, ie: u32, mode: SketchMode, sub_div_factor: u32, num_temp_files: usize) -> Self {
        // Attempt to raise the soft limit for open file descriptors to the hard limit.
        // This is often necessary when creating a large number of temporary files.
        match rlimit::getrlimit(rlimit::Resource::NOFILE) {
            Ok((soft, hard)) => {
                if soft < hard {
                    if let Err(e) = rlimit::setrlimit(rlimit::Resource::NOFILE, hard, hard) {
                        eprintln!("Warning: Failed to raise file descriptor limit: {}. This might cause issues if num_temp_files is large.", e);
                    }
                }
            },
            Err(e) => {
                 eprintln!("Warning: Could not get the current file descriptor limit: {}.", e);
            }
        }

        let s = 1u64 << ils;
        let fingerprint_range = 1u64 << iw;
        
        let temp_dir = tempfile::Builder::new()
            .prefix("pafsketch-builder-")
            .tempdir()
            .expect("Could not create temporary directory for index construction");

        // Create a pool of writers for the temporary bucket files.
        let writers = (0..num_temp_files)
            .map(|i| {
                let file_path = temp_dir.path().join(i.to_string());
                let file = File::create(file_path).expect("Could not create temporary bucket file. Check file permissions and ulimit -n.");
                Arc::new(Mutex::new(BufWriter::new(file)))
            })
            .collect();

        Self {
            ls: ils, w: iw, k: ik, e: ie, s,
            sketch_size: s,
            fingerprint_range,
            mi: u64::MAX,
            b: sub_div_factor,
            b_mask: (1u32 << sub_div_factor) - 1,
            temp_dir,
            writers,
            genome_numbers: AtomicU64::new(0),
            sketch_mode: mode,
        }
    }

    /// Finalizes the index construction.
    /// This method reads all temporary files, sorts and compresses the GID lists,
    /// and constructs the final, queryable Index struct.
    pub fn into_final_index(self) -> Index {
        println!("Finalizing index: flushing temporary files...");
        // Ensure all writers are flushed before reading from the temporary files.
        for writer_arc in self.writers.iter() {
            writer_arc.lock().unwrap().flush().expect("Failed to flush temporary bucket file.");
        }

        let total_size = ((self.fingerprint_range * self.s) >> self.b) as usize;
        let mut intermediate_data: Vec<Vec<u32>> = vec![Vec::new(); total_size];

        println!("Reading and aggregating data from temporary files...");
        // Iterate over the temporary files and populate the intermediate data vector.
        for i in 0..self.writers.len() {
            let file_path = self.temp_dir.path().join(i.to_string());
            let file = File::open(file_path).expect("Could not open temporary bucket file for reading.");
            let mut reader = BufReader::new(file);

            // Read all (index, gid) pairs from the current bucket file.
            loop {
                match reader.read_u64::<LittleEndian>() {
                    Ok(new_index) => {
                        let encoded_gid = reader.read_u32::<LittleEndian>()
                            .expect("Incomplete GID/position pair in bucket file. File may be corrupt.");
                        if (new_index as usize) < intermediate_data.len() {
                            intermediate_data[new_index as usize].push(encoded_gid);
                        }
                    },
                    Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                        // Reached the end of the current temporary file.
                        break;
                    }
                    Err(e) => {
                        panic!("Error reading from bucket file: {}", e);
                    }
                }
            }
        }
        
        println!("Compressing GID lists...");
        // In parallel, sort, delta-encode, and compress each GID list.
        let data: Vec<(u32, Vec<u8>)> = intermediate_data
            .into_par_iter()
            .map(|mut gids| {
                if gids.is_empty() {
                    return (0, Vec::new());
                }

                gids.sort_unstable();
                let original_count = gids.len() as u32;

                // Delta encode the sorted GIDs to improve compression ratio.
                let mut prev = gids[0];
                for i in 1..gids.len() {
                    let current = gids[i];
                    gids[i] = current.wrapping_sub(prev);
                    prev = current;
                }
                
                // stream-vbyte requires a buffer that is larger than the input.
                // 5 bytes per u32 is a safe upper bound.
                let mut encoded_data = vec![0; 5 * gids.len()];
                let encoded_len = stream_vbyte::encode::encode::<stream_vbyte::scalar::Scalar>(&gids, &mut encoded_data);
                encoded_data.truncate(encoded_len);
                
                (original_count, encoded_data)
            })
            .collect();
        
        println!("Index construction complete.");
        // The temp_dir is automatically removed when `self` is dropped.

        Index {
            ls: self.ls, w: self.w, k: self.k, e: self.e, s: self.s,
            sketch_size: self.sketch_size,
            fingerprint_range: self.fingerprint_range,
            mi: self.mi,
            b: self.b,
            data,
            genome_numbers: self.genome_numbers.load(Ordering::Relaxed),
            sketch_mode: self.sketch_mode,
        }
    }

    pub fn index_file_line_by_line(&self, filestr: &str) {
        let file = File::open(filestr).expect("Could not open file");
        let reader = BufReader::new(file);
        
        let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();
        self.genome_numbers.store(lines.len() as u64, Ordering::SeqCst);

        lines.into_par_iter().enumerate().for_each(|(i, line)| {
            let seq_id = i as u32;
            let mut sketch = vec![self.mi; self.s as usize];
            self.compute_sketch(SketchInput::Sequence(line.as_bytes()), &mut sketch);
            if !sketch.iter().all(|&x| x == self.mi) {
                self.insert_sketch(&sketch, seq_id);
            }
        });
    }

    pub fn index_file_of_files(&self, fof_path: &str) {
        let file = File::open(fof_path).expect("Could not open file of files");
        let reader = BufReader::new(file);
        
        let files_to_process: Vec<String> = reader.lines().filter_map(Result::ok).collect();
        self.genome_numbers.store(files_to_process.len() as u64, Ordering::SeqCst);

        files_to_process.into_par_iter().enumerate().for_each(|(i, filepath)| {
            let seq_id = i as u32;
            let mut sketch = vec![self.mi; self.s as usize];
            self.compute_sketch(SketchInput::FilePath(&filepath), &mut sketch);
            if !sketch.iter().all(|&x| x == self.mi) {
                self.insert_sketch(&sketch, seq_id);
            }
        });
    }

    fn compute_sketch(&self, input: SketchInput, sketch: &mut Vec<u64>) {
        sketch.fill(self.mi);
        let mut empty_cell = self.s;
        let is_valid_dna = |c: &u8| matches!(*c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');

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
                if let Ok(mut file_reader) = parse_fastx_file(path) {
                    while let Some(Ok(record)) = file_reader.next() {
                        update_sketch(&record.seq());
                    }
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
    
    fn insert_sketch(&self, sketch: &[u64], genome_id: Gid) {
        for (i, &fp) in sketch.iter().enumerate() {
            if fp < self.fingerprint_range {
                let flat_index = (fp * self.sketch_size) + i as u64;
                let new_index = (flat_index >> self.b) as usize;
                let offset = (flat_index as u32) & self.b_mask;

                if genome_id < (1 << (32 - self.b)) {
                    let encoded_gid = (genome_id << self.b) | offset;
                    
                    if !self.writers.is_empty() {
                        // Use a hash of the index to pick a bucket file, to distribute writes.
                        let bucket_index = new_index % self.writers.len();
                        let mut writer_guard = self.writers[bucket_index].lock().unwrap();
                        
                        // Write the final index and the encoded value to the bucket file.
                        writer_guard.write_u64::<LittleEndian>(new_index as u64).unwrap();
                        writer_guard.write_u32::<LittleEndian>(encoded_gid).unwrap();
                    }
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
        if hashed == self.mi { return self.mi; }
    
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
    pub fn from_file(filestr: &str) -> std::io::Result<Self> {
        println!("Loading index from disk: {}...", filestr);
        let mut file = File::open(filestr)?;
        let mut compressed_buffer = Vec::new();
        file.read_to_end(&mut compressed_buffer)?;
        
        let uncompressed_buffer = zstd::decode_all(&compressed_buffer[..])?;
        let mut reader = std::io::Cursor::new(uncompressed_buffer);

        let ls = reader.read_u32::<LittleEndian>()?;
        let k = reader.read_u32::<LittleEndian>()?;
        let w = reader.read_u32::<LittleEndian>()?;
        let e = reader.read_u32::<LittleEndian>()?;
        let fingerprint_range = reader.read_u64::<LittleEndian>()?;
        let genome_numbers = reader.read_u64::<LittleEndian>()?;
        let sketch_mode = match reader.read_u8()? {
            1 => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let b = reader.read_u32::<LittleEndian>()?;
        
        let s = 1u64 << ls;
        let total_size = ((fingerprint_range * s) >> b) as usize;
        let mut data = Vec::with_capacity(total_size);
        
        for _ in 0..total_size {
            let original_count = reader.read_u32::<LittleEndian>()?;
            let compressed_len = reader.read_u32::<LittleEndian>()?;
            
            let mut compressed_gids = vec![0u8; compressed_len as usize];
            reader.read_exact(&mut compressed_gids)?;
            
            data.push((original_count, compressed_gids));
        }
        
        println!("Finished loading index.");
        Ok(Self {
            ls, w, k, e, s,
            sketch_size: s,
            fingerprint_range,
            mi: u64::MAX,
            b,
            data,
            genome_numbers,
            sketch_mode,
        })
    }

   pub fn dump_index_disk(&self, filestr: &str, zstd_level: i32) -> std::io::Result<()> {
    let file = File::create(filestr)?;
    let mut writer = zstd::stream::write::Encoder::new(file, zstd_level)?.auto_finish();

    // Write metadata.
    writer.write_u32::<LittleEndian>(self.ls)?;
    writer.write_u32::<LittleEndian>(self.k)?;
    writer.write_u32::<LittleEndian>(self.w)?;
    writer.write_u32::<LittleEndian>(self.e)?;
    writer.write_u64::<LittleEndian>(self.fingerprint_range)?;
    writer.write_u64::<LittleEndian>(self.genome_numbers)?;
    writer.write_u8(self.sketch_mode as u8)?;
    writer.write_u32::<LittleEndian>(self.b)?;

    // Write the compressed index data.
    for (original_count, compressed_gids) in &self.data {
        // Write the original number of GIDs in the list.
        writer.write_u32::<LittleEndian>(*original_count)?;
        // Write the length of the compressed byte vector.
        writer.write_u32::<LittleEndian>(compressed_gids.len() as u32)?;
        // Write the compressed data itself.
        writer.write_all(compressed_gids)?;
    }
    
    Ok(())
}

    pub fn print_stats(&self) {
        println!("\nINDEX STATISTICS");
        println!("=========================================================");

        // Use the stored original counts for statistics, avoiding decompression.
        let sizes: Vec<u64> = self.data.par_iter()
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

        let variance = sizes.iter().map(|&value| {
            let diff = value as f64 - mean;
            diff * diff
        }).sum::<f64>() / num_lists;
        let std_dev = variance.sqrt();
        
        let sum_sq: f64 = sizes.iter().map(|&s| (s as f64).powi(2)).sum();
        let l2_norm = sum_sq.sqrt();

        let sum_cubed: f64 = sizes.iter().map(|&s| (s as f64).powi(3)).sum();
        let l3_norm = sum_cubed.cbrt();

        println!("Number of non-empty lists: {}", sizes.len());
        println!("Total elements in index: {}", sum);
        println!("Min list size: {}", min_size);
        println!("Max list size: {}", max_size);
        println!("Mean list size: {:.4}", mean);
        println!("Standard Deviation of list sizes: {:.4}", std_dev);
        println!("L2 Norm of list sizes: {:.4}", l2_norm);
        println!("L3 Norm of list sizes: {:.4}", l3_norm);
    }
    
    /// Performs an all-versus-all comparison between this index (reference) and another (query).
    pub fn all_vs_all_comparison(&self, query_index: &Index, threshold: f64, is_matrix: bool, zstd_level: i32, output_file: &str) {
        println!("\n--- Starting All-Versus-All Comparison ---");
        let start_time = std::time::Instant::now();
    
        let num_queries = query_index.genome_numbers as usize;
        let num_refs = self.genome_numbers as usize;
    
        let total_indices = self.data.len();
        if total_indices != query_index.data.len() {
            panic!("Incompatible indices for comparison: different total sizes.");
        }
        
        // Parallel map-reduce for score calculation to avoid lock contention
        let final_scores_map: HashMap<(u32, u32), u32> = (0..total_indices)
            .into_par_iter()
            .fold(
                || HashMap::new(),
                |mut local_scores, flat_index| {
                    let (ref_count, ref_compressed) = &self.data[flat_index];
                    if *ref_count > 0 {
                        let (query_count, query_compressed) = &query_index.data[flat_index];
                        if *query_count > 0 {
                            // Decompress and decode reference GIDs
                            let mut ref_decoded_deltas = vec![0u32; *ref_count as usize];
                            decode::<Scalar>(ref_compressed, *ref_count as usize, &mut ref_decoded_deltas);
                            for j in 1..ref_decoded_deltas.len() {
                                ref_decoded_deltas[j] = ref_decoded_deltas[j].wrapping_add(ref_decoded_deltas[j-1]);
                            }
    
                            // Decompress and decode query GIDs
                            let mut query_decoded_deltas = vec![0u32; *query_count as usize];
                            decode::<Scalar>(query_compressed, *query_count as usize, &mut query_decoded_deltas);
                            for j in 1..query_decoded_deltas.len() {
                                query_decoded_deltas[j] = query_decoded_deltas[j].wrapping_add(query_decoded_deltas[j-1]);
                            }
                            
                            let b_mask = (1u32 << self.b) - 1;
                            let b_shift = self.b;
                            
                            let mut ref_gids_by_offset: HashMap<u32, Vec<u32>> = HashMap::new();
                            for &encoded_gid in &ref_decoded_deltas {
                                let offset = encoded_gid & b_mask;
                                let gid = encoded_gid >> b_shift;
                                ref_gids_by_offset.entry(offset).or_default().push(gid);
                            }
    
                            let mut query_gids_by_offset: HashMap<u32, Vec<u32>> = HashMap::new();
                            for &encoded_gid in &query_decoded_deltas {
                                let offset = encoded_gid & b_mask;
                                let gid = encoded_gid >> b_shift;
                                query_gids_by_offset.entry(offset).or_default().push(gid);
                            }
                            
                            // Find common offsets and increment scores for all pairs
                            for (offset, ref_gids) in ref_gids_by_offset.iter() {
                                if let Some(query_gids) = query_gids_by_offset.get(offset) {
                                    for &query_gid in query_gids {
                                        for &ref_gid in ref_gids {
                                            *local_scores.entry((query_gid, ref_gid)).or_insert(0) += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    local_scores
                },
            )
            .reduce(
                || HashMap::new(),
                |mut a, b| {
                    for (k, v) in b {
                        *a.entry(k).or_insert(0) += v;
                    }
                    a
                },
            );
    
        println!("Score calculation finished in {} seconds. Formatting and writing results...", start_time.elapsed().as_secs());
    
        // Convert HashMap to a format that is easier to process for output (grouped by query GID)
        let mut scores_by_query: Vec<Vec<(Gid, u32)>> = vec![Vec::new(); num_queries];
        for ((query_gid, ref_gid), count) in final_scores_map {
            if (query_gid as usize) < num_queries {
                scores_by_query[query_gid as usize].push((ref_gid, count));
            }
        }
        
        let mut writer = self.create_writer(output_file, zstd_level);
    
        for query_gid in 0..num_queries {
            let mut line = String::new();
            // Using generic "query_GID" as we don't have file names here
            write!(&mut line, "query_{}\t", query_gid).unwrap();
            
            let res = &scores_by_query[query_gid];
    
            if is_matrix {
                // Reconstruct the full row for the matrix format
                let mut full_row = vec![0.0; num_refs];
                for &(ref_gid, count) in res {
                    if (ref_gid as usize) < num_refs {
                        full_row[ref_gid as usize] = count as f64 / self.s as f64;
                    }
                }
                let similarities: Vec<String> = full_row.iter().map(|&similarity| {
                    if similarity >= threshold {
                        format!("{:.4}", similarity)
                    } else {
                        "0.0000".to_string()
                    }
                }).collect();
                line.push_str(&similarities.join(","));
            } else { // Sparse format
                let mut hits: Vec<(Gid, f64)> = res.iter()
                    .map(|&(ref_gid, count)| {
                        let sim = count as f64 / self.s as f64;
                        (ref_gid, sim)
                    })
                    .filter(|&(_, sim)| sim >= threshold)
                    .collect();
                
                hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                
                line.push_str(&hits.iter().map(|(gid, sim)| format!("{}:{:.4}", gid, sim)).collect::<Vec<String>>().join(","));
            }
            writeln!(writer, "{}", line).unwrap();
        }
        
        println!("All-vs-all comparison complete. Results written to {}.", output_file);
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