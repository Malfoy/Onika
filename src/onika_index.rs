use crate::utils::*;
use ahash::AHashMap;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use flate2::read::GzDecoder;
use io::BufWriter;
use needletail::parse_fastx_reader;
use nthash::NtHashIterator;
use num_traits::{Float, ToPrimitive};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicU32, AtomicU64, AtomicUsize, Ordering},
    mpsc, Arc, Mutex,
};
use std::thread;

use f128::f128;
use std::fmt::Write as FmtWrite;
use stream_vbyte::{decode::decode, scalar::Scalar};
use tempfile::NamedTempFile;
use zstd::stream::write::Encoder as ZstdEncoder;

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

/// Represents the source of the index data, either fully in memory or on disk.
pub enum IndexData {
    InMemory(Vec<Vec<(u32, Vec<u8>)>>),
    OnDisk {
        handle: Option<Arc<NamedTempFile>>,
        path: String,
        toc: Vec<(u64, u64)>, // offset, length
    },
}

/// The final, immutable, lock-free index. It contains `s` inverted indexes, one for each sketch position.
pub struct Index {
    ls: u32,
    w: u32,
    k: u32,
    e: u32,
    s: u64,
    sketch_size: u64,
    fingerprint_range: u64,
    mi: u64,
    data: IndexData,
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
    mi: u32,
    // A single, flat vector representing a 2D matrix (pos x doc_id).
    intermediate_data: Vec<AtomicU32>,
    // Store num_docs to know the matrix width.
    num_docs: usize, 
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
    pub fn new(ils: u32, ik: u32, iw: u32, ie: u32, mode: SketchMode) -> Self {
        let s = 1u64 << ils;
        let fingerprint_range = 1u64 << iw;

        Self {
            ls: ils,
            w: iw,
            k: ik,
            e: ie,
            s,
            sketch_size: s,
            fingerprint_range,
            mi: u32::MAX,
            intermediate_data: Vec::new(),
            num_docs: 0,
            genome_numbers: AtomicU64::new(0),
            empty_sketches_count: AtomicU64::new(0),
            sketch_mode: mode,
        }
    }

    /// Finalizes the index construction by building inverted indexes from the fingerprint lists.
    pub fn into_final_index(self, output_path: &str, compress: bool, zstd_level: i32) -> io::Result<Index> {
        println!("Building inverted indexes from fingerprint lists (in parallel)...");
        
        let output_dir = Path::new(output_path).parent().unwrap_or_else(|| Path::new("."));
        let temp_file = tempfile::Builder::new().prefix(".tmp-index-").tempfile_in(output_dir)?;
        let path_str = temp_file.path().to_str().unwrap().to_string();

        let position_buffers: Vec<Vec<u8>> = self.intermediate_data
            .par_chunks(self.num_docs)
            .map(|fingerprint_chunk| {
                let mut inverted_index: Vec<Vec<u32>> = vec![Vec::new(); self.fingerprint_range as usize];
                for (doc_id, fp_atomic) in fingerprint_chunk.iter().enumerate() {
                    let fp = fp_atomic.load(Ordering::Relaxed);
                    if fp != self.mi {
                        inverted_index[fp as usize].push(doc_id as u32);
                    }
                }
                
                let mut uncompressed_buffer = Vec::new();
                for mut gids in inverted_index {
                    let original_count = gids.len() as u32;
                    let gids_bytes = if compress && !gids.is_empty() {
                        gids.sort_unstable();
                        let mut prev = gids[0];
                        for i in 1..gids.len() {
                            let current = gids[i];
                            gids[i] = current.wrapping_sub(prev);
                            prev = current;
                        }
                        let mut encoded_data = vec![0; 5 * gids.len()];
                        let encoded_len = stream_vbyte::encode::encode::<Scalar>(&gids, &mut encoded_data);
                        encoded_data.truncate(encoded_len);
                        encoded_data
                    } else {
                        let mut bytes = Vec::with_capacity(gids.len() * 4);
                        for gid in gids {
                            bytes.write_u32::<LittleEndian>(gid).unwrap();
                        }
                        bytes
                    };
                    uncompressed_buffer.write_u32::<LittleEndian>(original_count).unwrap();
                    uncompressed_buffer.write_u32::<LittleEndian>(gids_bytes.len() as u32).unwrap();
                    uncompressed_buffer.write_all(&gids_bytes).unwrap();
                }
                zstd::encode_all(&uncompressed_buffer[..], zstd_level).unwrap()
            })
            .collect();
        
        let mut toc = Vec::with_capacity(self.s as usize);
        
        {
            let mut file = BufWriter::new(temp_file.as_file());
            let header_size = 4 + 4 + 4 + 4 + 8 + 8 + 1 + 1;
            let toc_size = (self.s * 16) as u64;
            let mut current_offset = header_size + toc_size;
            for buffer in &position_buffers {
                let len = buffer.len() as u64;
                toc.push((current_offset, len));
                current_offset += len;
            }
            file.seek(SeekFrom::Start(0))?;
            file.write_u32::<LittleEndian>(self.ls)?;
            file.write_u32::<LittleEndian>(self.k)?;
            file.write_u32::<LittleEndian>(self.w)?;
            file.write_u32::<LittleEndian>(self.e)?;
            file.write_u64::<LittleEndian>(self.fingerprint_range)?;
            file.write_u64::<LittleEndian>(self.genome_numbers.load(Ordering::Relaxed))?;
            file.write_u8(self.sketch_mode as u8)?;
            file.write_u8(compress as u8)?;
            for (offset, len) in &toc {
                file.write_u64::<LittleEndian>(*offset)?;
                file.write_u64::<LittleEndian>(*len)?;
            }
            for buffer in position_buffers {
                file.write_all(&buffer)?;
            }
        }
        
        let empty_count = self.empty_sketches_count.load(Ordering::Relaxed);
        if empty_count > 0 {
            eprintln!("\nWarning: {} input documents resulted in empty sketches and were not indexed.", empty_count);
        }
        println!("Index construction complete.");

        let (file_handle, path) = temp_file.into_parts();
        
        Ok(Index {
            ls: self.ls,
            w: self.w,
            k: self.k,
            e: self.e,
            s: self.s,
            sketch_size: self.sketch_size,
            fingerprint_range: self.fingerprint_range,
            mi: self.mi as u64,
            data: IndexData::OnDisk { handle: Some(Arc::new(NamedTempFile::from_parts(file_handle, path))), path: path_str, toc },
            genome_numbers: self.genome_numbers.load(Ordering::Relaxed),
            sketch_mode: self.sketch_mode,
            is_compressed: compress,
        })
    }




pub fn index_fasta_file(&mut self, filestr: &str) {
        let reader = open_compressed_file(filestr).expect("Could not open FASTA/FASTQ file");
        let mut fastx_reader =
            needletail::parse_fastx_reader(reader).expect("Invalid FASTA/FASTQ format");

        // Use a simple loop to collect sequences, which is the correct pattern.
        let mut sequences = Vec::new();
        while let Some(record) = fastx_reader.next() {
            let seqrec = record.expect("Invalid record in FASTA/FASTQ file");
            sequences.push(seqrec.seq().to_vec());
        }
        
        let num_docs = sequences.len();
        self.num_docs = num_docs;
        self.genome_numbers.store(num_docs as u64, Ordering::SeqCst);
        
        let total_size = self.s as usize * num_docs;
        self.intermediate_data.resize_with(total_size, || AtomicU32::new(self.mi));
        
        // This can now be safely shared across threads.
        let data_arc = Arc::new(&self.intermediate_data);

        sequences.into_par_iter().enumerate().for_each(|(i, seq)| {
            let seq_id = i as u32;
            let mut sketch = vec![self.mi as u64; self.s as usize];
            self.compute_sketch(SketchInput::Sequence(&seq), &mut sketch);

            if sketch.iter().all(|&x| x == self.mi as u64) {
                self.empty_sketches_count.fetch_add(1, Ordering::Relaxed);
                return;
            }

            for (pos, &fp) in sketch.iter().enumerate() {
                if fp < self.fingerprint_range {
                    let index = pos * num_docs + (seq_id as usize);
                    data_arc[index].store(fp as u32, Ordering::Relaxed);
                }
            }
        });
    }




    pub fn index_file_of_files(&mut self, fof_path: &str) {
        let reader = open_compressed_file(fof_path).expect("Could not open file of files");
        let files_to_process: Vec<String> = reader.lines().filter_map(Result::ok).collect();
        let num_docs = files_to_process.len();
        
        self.num_docs = num_docs;
        self.genome_numbers.store(num_docs as u64, Ordering::SeqCst);

        let total_size = self.s as usize * num_docs;
        self.intermediate_data.resize_with(total_size, || AtomicU32::new(self.mi));
        let data_arc = Arc::new(&self.intermediate_data);

        files_to_process.into_par_iter().enumerate().for_each(|(i, filepath)| {
            let seq_id = i as u32;
            let mut sketch = vec![self.mi as u64; self.s as usize];
            self.compute_sketch(SketchInput::FilePath(&filepath), &mut sketch);

            if sketch.iter().all(|&x| x == self.mi as u64) {
                self.empty_sketches_count.fetch_add(1, Ordering::Relaxed);
                return;
            }

            for (pos, &fp) in sketch.iter().enumerate() {
                if fp < self.fingerprint_range {
                    let index = pos * num_docs + (seq_id as usize);
                    data_arc[index].store(fp as u32, Ordering::Relaxed);
                }
            }
        });
    }

    fn compute_sketch(&self, input: SketchInput, sketch: &mut Vec<u64>) {
        sketch.fill(self.mi as u64);
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
                            if sketch[bucket_id] == self.mi as u64 {
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
                if *val != self.mi as u64 {
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
                if sketch[i] != self.mi as u64 {
                    let hash_pos = (hash_family(sketch[i], step) % size as u64) as usize;
                    if hash_pos < sketch.len() && sketch[hash_pos] == self.mi as u64 {
                        changes.push((hash_pos, sketch[i]));
                    }
                }
            }
            if changes.is_empty() && empty_cell > 0 { break; }

            for (pos, val) in changes {
                if sketch[pos] == self.mi as u64 {
                    sketch[pos] = val;
                    empty_cell -= 1;
                    if empty_cell == 0 { return; }
                }
            }
            step += 1;
        }
    }

    fn get_perfect_fingerprint(&self, hashed: u64) -> u64 {
        if hashed == self.mi as u64 { return self.mi as u64; }
        let b = f128::from(hashed >> 32);
        let twopowern = f128::from(1u128 << 32);
        let ratio = f128::from(self.e) / f128::from(self.s);
        let top = (twopowern - b).powf(ratio);
        let bottom = twopowern.powf(ratio);
        let scaled_top = f128::from(self.fingerprint_range) * top;
        let result_f = f128::from(self.fingerprint_range) - (scaled_top / bottom);
        if result_f < f128::from(0.0) { 0 } else { result_f.floor().to_u64().unwrap_or(0) }
    }
}

impl Index {
    pub fn from_file(filestr: &str) -> io::Result<Self> {
        println!("Loading index from disk: {}...", filestr);
        let mut file = BufReader::new(File::open(filestr)?);

        let ls = file.read_u32::<LittleEndian>()?;
        let k = file.read_u32::<LittleEndian>()?;
        let w = file.read_u32::<LittleEndian>()?;
        let e = file.read_u32::<LittleEndian>()?;
        let fingerprint_range = file.read_u64::<LittleEndian>()?;
        let genome_numbers = file.read_u64::<LittleEndian>()?;
        let sketch_mode = match file.read_u8()? {
            1 => SketchMode::Perfect,
            _ => SketchMode::Default,
        };
        let is_compressed = file.read_u8()? != 0;
        let s = 1u64 << ls;

        let mut toc = Vec::with_capacity(s as usize);
        for _ in 0..s {
            toc.push((file.read_u64::<LittleEndian>()?, file.read_u64::<LittleEndian>()?));
        }

        println!("Finished loading index metadata.");
        Ok(Self {
            ls, w, k, e, s, sketch_size: s, fingerprint_range, mi: u64::MAX,
            data: IndexData::OnDisk { handle: None, path: filestr.to_string(), toc },
            genome_numbers, sketch_mode, is_compressed,
        })
    }

    pub fn dump_index_disk(self, filestr: &str) -> io::Result<()> {
        match self.data {
            IndexData::OnDisk{ handle, path, ..} => {
                if let Some(temp_handle) = handle {
                    match Arc::try_unwrap(temp_handle) {
                        Ok(temp_file) => { temp_file.persist(filestr)?; },
                        Err(_) => return Err(io::Error::new(io::ErrorKind::Other, "Failed to unwrap Arc for temp file")),
                    }
                } else {
                    std::fs::copy(&path, filestr)?;
                }
            },
            IndexData::InMemory(_) => {
                panic!("dump_index_disk called on an in-memory index, which is not supported by the new design.");
            }
        }
        Ok(())
    }

    pub fn print_stats(&self) {
        if let IndexData::InMemory(data) = &self.data {
            println!("\nINDEX STATISTICS (In-Memory)");
            println!("=========================================================");

            let sizes: Vec<u64> = data
                .par_iter()
                .flat_map(|pos_index| pos_index.par_iter().map(|(count, _)| *count as u64))
                .filter(|&size| size > 0)
                .collect();

            if sizes.is_empty() { println!("Index is empty."); return; }
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
        } else {
            println!("\nINDEX STATISTICS (On-Disk)");
            println!("=========================================================");
            println!("Statistics are not computed for on-disk indexes to avoid slow file reads.");
            println!("Sketch Size (s): {}", self.sketch_size);
            println!("Fingerprint Range (2^w): {}", self.fingerprint_range);
            println!("Number of Genomes: {}", self.genome_numbers);
        }
    }
    
    pub fn all_vs_all_comparison(&self, query_index: &Index, threshold: f64, is_matrix: bool, zstd_level: i32, output_file: &str, num_io_threads: usize, io_buffer_size: usize, num_compute_threads: usize) {
        println!("\n--- Starting All-Versus-All Comparison ---");
        let start_time = std::time::Instant::now();

        let num_queries = query_index.genome_numbers as usize;
        let num_refs = self.genome_numbers as usize;

        let min_required_score = if threshold > 0.0 { (threshold * self.s as f64).ceil() as u32 } else { 0 };
        let perform_pruning = min_required_score > 1;
        if perform_pruning { println!("Similarity threshold of {} requires at least {} shared fingerprints. Applying pruning.", threshold, min_required_score); }

        let num_shards = 1024;
        let scores_shards: Vec<Mutex<AHashMap<(u32, u32), u32>>> =
            (0..num_shards).map(|_| Mutex::new(AHashMap::default())).collect();

        let pool = ThreadPoolBuilder::new()
            .num_threads(num_compute_threads)
            .build()
            .unwrap();

        pool.install(|| {
            if let (IndexData::InMemory(ref_data), IndexData::InMemory(query_data)) = (&self.data, &query_index.data) {
                (0..self.sketch_size as usize).into_par_iter().for_each(|pos| {
                    let remaining_hits_possible = self.sketch_size - pos as u64 - 1;
                    let ref_pos_index = &ref_data[pos];
                    let query_pos_index = &query_data[pos];

                    for fp in 0..self.fingerprint_range as usize {
                        let (ref_count, ref_bytes) = &ref_pos_index[fp];
                        if *ref_count > 0 {
                            let (query_count, query_bytes) = &query_pos_index[fp];
                            if *query_count > 0 {
                                let ref_gids = decode_gid_list(*ref_count, ref_bytes, self.is_compressed);
                                let query_gids = decode_gid_list(*query_count, query_bytes, query_index.is_compressed);
                                for &query_gid in &query_gids {
                                    for &ref_gid in &ref_gids {
                                        let key = (query_gid, ref_gid);
                                        let shard_index = key.0 as usize % num_shards;
                                        let mut shard = scores_shards[shard_index].lock().unwrap();
                                        if !perform_pruning {
                                            *shard.entry(key).or_insert(0) += 1;
                                        } else if let Some(score) = shard.get_mut(&key) {
                                            *score += 1;
                                            if *score as u64 + remaining_hits_possible < min_required_score as u64 { shard.remove(&key); }
                                        } else if 1 + remaining_hits_possible >= min_required_score as u64 {
                                            shard.insert(key, 1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            } else {
                type PosData = (usize, Vec<(u32, Vec<u8>)>, Vec<(u32, Vec<u8>)>);
                let (tx, rx) = mpsc::sync_channel::<PosData>(io_buffer_size);
                
                let position_counter = Arc::new(AtomicUsize::new(0));
                let mut producer_handles = Vec::new();

                for _ in 0..num_io_threads {
                    let tx_clone = tx.clone();
                    let ref_path = self.path_for_reading().to_string();
                    let query_path = query_index.path_for_reading().to_string();
                    let ref_toc = self.toc_for_reading().to_vec();
                    let query_toc = query_index.toc_for_reading().to_vec();
                    let s = self.s;
                    let fp_range = self.fingerprint_range;
                    let counter_clone = Arc::clone(&position_counter);

                    producer_handles.push(thread::spawn(move || {
                        let mut ref_file = File::open(ref_path).expect("Cannot open reference index file in producer.");
                        let mut query_file = File::open(query_path).expect("Cannot open query index file in producer.");
                        loop {
                            let pos = counter_clone.fetch_add(1, Ordering::SeqCst);
                            if pos >= s as usize { break; }

                            let ref_pos_data = read_pos_from_file(&mut ref_file, &ref_toc, pos, fp_range);
                            let query_pos_data = read_pos_from_file(&mut query_file, &query_toc, pos, fp_range);
                            if tx_clone.send((pos, ref_pos_data, query_pos_data)).is_err() { break; }
                        }
                    }));
                }
                drop(tx);

                rx.into_iter().par_bridge().for_each(|(pos, ref_pos_index, query_pos_index)| {
                    let remaining_hits_possible = self.sketch_size - pos as u64 - 1;
                    (0..self.fingerprint_range as usize).into_par_iter().for_each(|fp| {
                        let (ref_count, ref_bytes) = &ref_pos_index[fp];
                        if *ref_count > 0 {
                            let (query_count, query_bytes) = &query_pos_index[fp];
                            if *query_count > 0 {
                                let ref_gids = decode_gid_list(*ref_count, ref_bytes, self.is_compressed);
                                let query_gids = decode_gid_list(*query_count, query_bytes, query_index.is_compressed);
                                for &query_gid in &query_gids {
                                    for &ref_gid in &ref_gids {
                                        let key = (query_gid, ref_gid);
                                        let shard_index = key.0 as usize % num_shards;
                                        let mut shard = scores_shards[shard_index].lock().unwrap();
                                        if !perform_pruning {
                                            *shard.entry(key).or_insert(0) += 1;
                                        } else if let Some(score) = shard.get_mut(&key) {
                                            *score += 1;
                                            if *score as u64 + remaining_hits_possible < min_required_score as u64 { shard.remove(&key); }
                                        } else if 1 + remaining_hits_possible >= min_required_score as u64 {
                                            shard.insert(key, 1);
                                        }
                                    }
                                }
                            }
                        }
                    });
                });

                for handle in producer_handles {
                    handle.join().unwrap();
                }
            }
        });
        
        println!("Score calculation finished in {} seconds. Aggregating and formatting results...", start_time.elapsed().as_secs());

        let mut scores_by_query: Vec<Vec<(Gid, u32)>> = vec![Vec::new(); num_queries];
        for shard_mutex in scores_shards {
            let shard_map = shard_mutex.into_inner().unwrap();
            for ((query_gid, ref_gid), count) in shard_map {
                if count >= min_required_score && (query_gid as usize) < num_queries {
                    scores_by_query[query_gid as usize].push((ref_gid, count));
                }
            }
        }
        
        println!("Parallel formatting of results...");
        let output_lines: Vec<String> = (0..num_queries).into_par_iter().map(|query_gid| {
            let mut line = String::with_capacity(256);
            write!(&mut line, "query_{}\t", query_gid).unwrap();
            let mut res = scores_by_query[query_gid].clone();
            if is_matrix {
                let mut full_row = vec![0.0; num_refs];
                for &(ref_gid, count) in res.iter() {
                    if (ref_gid as usize) < num_refs { full_row[ref_gid as usize] = count as f64 / self.s as f64; }
                }
                line.push_str(&full_row.iter().map(|&sim| if sim >= threshold { format!("{:.4}", sim) } else { "0.0000".to_string() }).collect::<Vec<String>>().join(","));
            } else {
                res.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                line.push_str(&res.iter().map(|(gid, count)| format!("{}:{:.4}", gid, *count as f64 / self.s as f64)).collect::<Vec<String>>().join(","));
            }
            line
        }).collect();
        
        println!("Writing results to disk...");
        let mut writer = self.create_writer(output_file, zstd_level);
        for line in output_lines { writeln!(writer, "{}", line).unwrap(); }

        println!("All-vs-all comparison complete. Total time: {} seconds.", start_time.elapsed().as_secs());
    }

    fn create_writer<'a>(&self, output_file: &'a str, zstd_level: i32) -> Box<dyn Write + 'a> {
        let file = File::create(output_file).expect("Cannot create output file");
        if zstd_level > 0 { Box::new(ZstdEncoder::new(file, zstd_level).unwrap().auto_finish()) } else { Box::new(BufWriter::new(file)) }
    }

    fn path_for_reading(&self) -> &str {
        match &self.data {
            IndexData::OnDisk { path, .. } => path,
            IndexData::InMemory(_) => panic!("Cannot get file path for an in-memory index."),
        }
    }
    
    fn toc_for_reading(&self) -> &[(u64, u64)] {
        match &self.data {
            IndexData::OnDisk { toc, .. } => toc,
            IndexData::InMemory(_) => panic!("Cannot get TOC for an in-memory index."),
        }
    }
}

fn read_pos_from_file(file: &mut File, toc: &[(u64, u64)], pos: usize, fp_range: u64) -> Vec<(u32, Vec<u8>)> {
    let (offset, len) = toc[pos];
    if len == 0 { 
        return vec![(0, Vec::new()); fp_range as usize];
    }
    
    let mut compressed_buffer = vec![0; len as usize];
    file.seek(SeekFrom::Start(offset)).unwrap();
    file.read_exact(&mut compressed_buffer).unwrap();
    
    let buffer = zstd::decode_all(&compressed_buffer[..]).unwrap();
    
    let mut cursor = io::Cursor::new(&buffer);
    let mut pos_index = Vec::with_capacity(fp_range as usize);

    while (cursor.position() as usize) < buffer.len() {
        let original_count = cursor.read_u32::<LittleEndian>().unwrap();
        let gids_len = cursor.read_u32::<LittleEndian>().unwrap();
        let mut gids_bytes = vec![0u8; gids_len as usize];
        cursor.read_exact(&mut gids_bytes).unwrap();
        pos_index.push((original_count, gids_bytes));
    }
    pos_index
}

/// Helper to decode a GID list from bytes depending on whether it's compressed.
fn decode_gid_list(count: u32, bytes: &[u8], is_compressed: bool) -> Vec<u32> {
    if count == 0 { return Vec::new(); }
    if is_compressed {
        let mut gids = vec![0u32; count as usize];
        decode::<Scalar>(bytes, count as usize, &mut gids);
        for j in 1..gids.len() { gids[j] = gids[j].wrapping_add(gids[j - 1]); }
        gids
    } else {
        let mut gids = Vec::with_capacity(count as usize);
        let mut cursor = io::Cursor::new(bytes);
        for _ in 0..count { gids.push(cursor.read_u32::<LittleEndian>().unwrap()); }
        gids
    }
}