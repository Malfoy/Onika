use crate::utils::*;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use flate2::read::GzDecoder;
use io::BufWriter;
use nthash::NtHashIterator;
use num_traits::{Float, ToPrimitive};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rlimit;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicU32, AtomicU64, Ordering},
    Arc, Mutex,
};
use tempfile::Builder as TempBuilder;
use std::fs;
use f128::f128;
use std::fmt::Write as FmtWrite;
use tempfile::NamedTempFile;
use zstd::stream::{write::Encoder as ZstdEncoder, read::Decoder as ZstdDecoder};

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

/// Represents the source of the index data.
pub enum IndexData {
    OnDisk {
        handle: Option<Arc<NamedTempFile>>,
        path: String,
        toc: Vec<(u64, u64)>, // offset, length
    },
}

/// The final, immutable, lock-free index.
pub struct Index {
    _ls: u32,
    _w: u32,
    _k: u32,
    _e: u32,
    s: u64,
    sketch_size: u64,
    fingerprint_range: u64,
    _mi: u64,
    data: IndexData,
    genome_numbers: u64,
    _sketch_mode: SketchMode,
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
    intermediate_data: Vec<AtomicU32>,
    num_docs: usize,
    genome_numbers: AtomicU64,
    empty_sketches_count: AtomicU64,
    sketch_mode: SketchMode,
}

pub(crate) fn open_compressed_file(path: &str) -> io::Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let magic_bytes = reader.fill_buf()?;

    if magic_bytes.starts_with(&[0x1f, 0x8b]) {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else if magic_bytes.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
        let file = File::open(path)?;
        let decoder = zstd::stream::read::Decoder::new(file)?;
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(reader))
    }
}

impl IndexBuilder {
    pub fn new(ils: u32, ik: u32, iw: u32, ie: u32, mode: SketchMode) -> Self {
        let s = 1u64 << ils;
        let fingerprint_range = 1u64 << iw;

        Self {
            ls: ils, w: iw, k: ik, e: ie, s, sketch_size: s, fingerprint_range,
            mi: u32::MAX, intermediate_data: Vec::new(), num_docs: 0,
            genome_numbers: AtomicU64::new(0), empty_sketches_count: AtomicU64::new(0),
            sketch_mode: mode,
        }
    }

    pub fn into_final_index(self, output_path: &str, compress: bool, zstd_level: i32) -> io::Result<Index> {
        println!("Building inverted indexes from fingerprint lists (in parallel)...");
        
        let output_dir = Path::new(output_path).parent().unwrap_or_else(|| Path::new("."));
        let temp_file = tempfile::Builder::new().prefix(".tmp-index-").tempfile_in(output_dir)?;
        let path_str = temp_file.path().to_str().unwrap().to_string();

        let position_buffers: Vec<Vec<u8>> = self.intermediate_data
            .par_chunks(self.num_docs)
            .map(|fingerprint_chunk| {
                let mut inverted_index: Vec<Vec<Gid>> = vec![Vec::new(); self.fingerprint_range as usize];
                for (doc_id, fp_atomic) in fingerprint_chunk.iter().enumerate() {
                    let fp = fp_atomic.load(Ordering::Relaxed);
                    if fp != self.mi {
                        inverted_index[fp as usize].push(doc_id as Gid);
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
                        let encoded_len = stream_vbyte::encode::encode::<stream_vbyte::scalar::Scalar>(&gids, &mut encoded_data);
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
            _ls: self.ls, _w: self.w, _k: self.k, _e: self.e, s: self.s,
            sketch_size: self.sketch_size, fingerprint_range: self.fingerprint_range,
            _mi: self.mi as u64,
            data: IndexData::OnDisk { handle: Some(Arc::new(NamedTempFile::from_parts(file_handle, path))), path: path_str, toc },
            genome_numbers: self.genome_numbers.load(Ordering::Relaxed),
            _sketch_mode: self.sketch_mode, is_compressed: compress,
        })
    }

    pub fn index_fasta_file(&mut self, filestr: &str) {
        let reader = open_compressed_file(filestr).expect("Could not open FASTA/FASTQ file");
        let mut fastx_reader = needletail::parse_fastx_reader(reader).expect("Invalid FASTA/FASTQ format");
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
        let data_arc = Arc::new(&self.intermediate_data);

        sequences.into_par_iter().enumerate().for_each(|(i, seq)| {
            let seq_id = i as Gid;
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
        const CHUNK_SIZE: usize = 32;

        files_to_process
            .par_chunks(CHUNK_SIZE)
            .enumerate()
            .for_each(|(chunk_index, file_chunk)| {
                let base_index = chunk_index * CHUNK_SIZE;
                for (i_in_chunk, filepath) in file_chunk.iter().enumerate() {
                    let i = base_index + i_in_chunk;
                    let seq_id = i as Gid;
                    let mut sketch = vec![self.mi as u64; self.s as usize];
                    self.compute_sketch(SketchInput::FilePath(filepath), &mut sketch);

                    if sketch.iter().all(|&x| x == self.mi as u64) {
                        self.empty_sketches_count.fetch_add(1, Ordering::Relaxed);
                        continue;
                    }

                    for (pos, &fp) in sketch.iter().enumerate() {
                        if fp < self.fingerprint_range {
                            let index = pos * num_docs + (seq_id as usize);
                            data_arc[index].store(fp as u32, Ordering::Relaxed);
                        }
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
                        let fp = canon_hash;
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
                if let Ok(reader) = open_compressed_file(path) {
                    if let Ok(mut file_reader) = needletail::parse_fastx_reader(reader) {
                        while let Some(Ok(record)) = file_reader.next() {
                            update_sketch(&record.seq());
                        }
                    }
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
        let sketch_mode = match file.read_u8()? { 1 => SketchMode::Perfect, _ => SketchMode::Default };
        let is_compressed = file.read_u8()? != 0;
        let s = 1u64 << ls;

        let mut toc = Vec::with_capacity(s as usize);
        for _ in 0..s {
            toc.push((file.read_u64::<LittleEndian>()?, file.read_u64::<LittleEndian>()?));
        }

        println!("Finished loading index metadata.");
        Ok(Self {
            _ls: ls, _w: w, _k: k, _e: e, s, sketch_size: s, fingerprint_range, _mi: u64::MAX,
            data: IndexData::OnDisk { handle: None, path: filestr.to_string(), toc },
            genome_numbers, _sketch_mode: sketch_mode, is_compressed,
        })
    }

    pub fn dump_index_disk(self, filestr: &str) -> io::Result<()> {
        let IndexData::OnDisk { handle, path, .. } = self.data;
        if let Some(temp_handle) = handle {
            if let Ok(temp_file) = Arc::try_unwrap(temp_handle) {
                temp_file.persist(filestr)?;
            } else {
                return Err(io::Error::new(io::ErrorKind::Other, "Failed to unwrap Arc for temp file"));
            }
        } else {
            std::fs::copy(&path, filestr)?;
        }
        Ok(())
    }

    pub fn print_stats(&self) {
        println!("\nINDEX STATISTICS (On-Disk)");
        println!("=========================================================");
        println!("Statistics are not computed for on-disk indexes to avoid slow file reads.");
        println!("Sketch Size (s): {}", self.sketch_size);
        println!("Fingerprint Range (2^w): {}", self.fingerprint_range);
        println!("Number of Genomes: {}", self.genome_numbers);
    }

    pub fn all_vs_all_comparison(
        &self,
        query_index: &Index,
        threshold: f64,
        is_matrix: bool,
        zstd_level: i32,
        output_file: &str,
        num_threads: usize,
        shard_zstd_level: i32,
        temp_dir_path: Option<&str>,
    ) {
        if let Ok((soft_limit, hard_limit)) = rlimit::getrlimit(rlimit::Resource::NOFILE) {
            if soft_limit < hard_limit {
                println!("Attempting to increase open file limit from {} to {}.", soft_limit, hard_limit);
                if let Err(e) = rlimit::setrlimit(rlimit::Resource::NOFILE, hard_limit, hard_limit) {
                    eprintln!("Warning: Failed to increase open file limit: {}. This may cause errors.", e);
                }
            }
        }

        println!("\n--- Starting All-Versus-All Comparison (Low Memory Mode, Single Pass) ---");
        let start_time = std::time::Instant::now();

        let num_queries = query_index.genome_numbers as usize;
        let num_refs = self.genome_numbers as usize;

        let min_required_score = if threshold > 0.0 { (threshold * self.s as f64).ceil() as u32 } else { 0 };
        if min_required_score > 1 {
            println!("Similarity threshold of {} requires at least {} shared fingerprints.", threshold, min_required_score);
        }

        const NUM_SHARDS: usize = 2048;
        const SHARD_BUFFER_SIZE: usize = 256 * 1024;
        const LOCAL_BUFFER_FLUSH_THRESHOLD: usize = 1024;

        let mut writer = self.create_writer(output_file, zstd_level);
        
        let temp_dir = match temp_dir_path {
            Some(path) => {
                fs::create_dir_all(path).expect("Failed to create specified temporary directory path");
                TempBuilder::new().prefix("penta_comp_").tempdir_in(path)
            }
            None => TempBuilder::new().prefix("penta_comp_").tempdir(),
        }
        .expect("Failed to create temporary directory");

        println!("Using temporary directory for intermediate data: {:?}", temp_dir.path());
        if shard_zstd_level > 0 {
            println!("Using zstd level {} for temporary shard files.", shard_zstd_level);
        }

        enum ShardWriter {
            Uncompressed(BufWriter<File>),
            Compressed(ZstdEncoder<'static, BufWriter<File>>),
        }

        impl io::Write for ShardWriter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                match self {
                    ShardWriter::Uncompressed(writer) => writer.write(buf),
                    ShardWriter::Compressed(encoder) => encoder.write(buf),
                }
            }
            fn flush(&mut self) -> io::Result<()> {
                match self {
                    ShardWriter::Uncompressed(writer) => writer.flush(),
                    ShardWriter::Compressed(encoder) => encoder.flush(),
                }
            }
        }

        let mut shard_paths = Vec::with_capacity(NUM_SHARDS);
        let shard_writers: Vec<Mutex<ShardWriter>> = (0..NUM_SHARDS).map(|i| {
            let path = if shard_zstd_level > 0 {
                temp_dir.path().join(format!("shard_{}.bin.zst", i))
            } else {
                temp_dir.path().join(format!("shard_{}.bin", i))
            };
            let file = File::create(&path).expect("Failed to create shard file.");
            shard_paths.push(path);
            let writer = BufWriter::with_capacity(SHARD_BUFFER_SIZE, file);
            let shard_writer = if shard_zstd_level > 0 {
                let encoder = ZstdEncoder::new(writer, shard_zstd_level).expect("Failed to create zstd encoder");
                ShardWriter::Compressed(encoder)
            } else {
                ShardWriter::Uncompressed(writer)
            };
            Mutex::new(shard_writer)
        }).collect();

        let pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();

        let remaining_buffers: Vec<Vec<Vec<u64>>> = pool.install(|| {
            let num_pos = self.sketch_size as usize;
            let chunk_size = (num_pos + num_threads - 1) / num_threads;
            let pos_iterator: Vec<usize> = (0..num_pos).collect();
            let pos_chunks: Vec<&[usize]> = pos_iterator.chunks(chunk_size).collect();

            pos_chunks.into_par_iter().map(|pos_chunk| {
                let mut ref_file = File::open(&self.path_for_reading()).expect("Cannot open reference index file in thread.");
                let mut query_file = File::open(&query_index.path_for_reading()).expect("Cannot open query index file in thread.");
                let mut local_shard_buffers: Vec<Vec<u64>> = vec![Vec::new(); NUM_SHARDS];

                for &pos in pos_chunk {
                    let ref_pos_index = read_pos_from_file(&mut ref_file, &self.toc_for_reading(), pos, self.fingerprint_range);
                    let query_pos_index = read_pos_from_file(&mut query_file, &query_index.toc_for_reading(), pos, query_index.fingerprint_range);

                    for fp in 0..self.fingerprint_range as usize {
                        let (ref_count, ref_bytes) = &ref_pos_index[fp];
                        if *ref_count > 0 {
                            let (query_count, query_bytes) = &query_pos_index[fp];
                            if *query_count > 0 {
                                let ref_gids = decode_gid_list(*ref_count, ref_bytes, self.is_compressed);
                                let query_gids = decode_gid_list(*query_count, query_bytes, query_index.is_compressed);
                                for &query_gid in &query_gids {
                                    let shard_index = revhash64(query_gid as u64) as usize % NUM_SHARDS;
                                    for &ref_gid in &ref_gids {
                                        let key = (query_gid as u64) << 32 | (ref_gid as u64);
                                        local_shard_buffers[shard_index].push(key);
                                    }
                                }
                            }
                        }
                    }

                    for (shard_index, buffer) in local_shard_buffers.iter_mut().enumerate() {
                        if buffer.len() >= LOCAL_BUFFER_FLUSH_THRESHOLD {
                            if let Ok(mut writer_lock) = shard_writers[shard_index].try_lock() {
                                for &key in buffer.iter() {
                                    writer_lock.write_u64::<LittleEndian>(key).unwrap();
                                }
                                buffer.clear();
                            }
                        }
                    }
                }
                local_shard_buffers
            }).collect()
        });
        
        println!("Flushing remaining thread buffers in parallel...");
        (0..NUM_SHARDS).into_par_iter().for_each(|shard_index| {
            let mut writer_lock = shard_writers[shard_index].lock().unwrap();
            for thread_buffers in &remaining_buffers {
                let buffer = &thread_buffers[shard_index];
                if !buffer.is_empty() {
                    for &key in buffer {
                        writer_lock.write_u64::<LittleEndian>(key).unwrap();
                    }
                }
            }
        });

        for writer_mutex in shard_writers {
            let shard_writer = writer_mutex.into_inner().expect("Mutex lock failed");
            match shard_writer {
                ShardWriter::Uncompressed(mut writer) => writer.flush().expect("Failed to flush writer"),
                ShardWriter::Compressed(encoder) => { encoder.finish().expect("Failed to finish zstd stream"); },
            }
        }

        let aggregation_start_time = std::time::Instant::now();

        let aggregated_pairs: Vec<Vec<(u64, u32)>> = shard_paths.into_par_iter().filter_map(|path| {
            let file = File::open(&path).ok()?;
            if file.metadata().ok()?.len() == 0 { fs::remove_file(&path).ok(); return None; }

            let mut reader: Box<dyn Read> = if shard_zstd_level > 0 {
                Box::new(ZstdDecoder::new(file).ok()?)
            } else {
                Box::new(BufReader::new(file))
            };

            let mut all_bytes = Vec::new();
            reader.read_to_end(&mut all_bytes).ok()?;
            
            let mut keys = Vec::new();
            let mut cursor = io::Cursor::new(all_bytes);
            while let Ok(key) = cursor.read_u64::<LittleEndian>() { keys.push(key); }

            keys.par_sort_unstable();

            let mut results_for_shard = Vec::new();
            if keys.is_empty() { fs::remove_file(&path).ok(); return Some(results_for_shard); }

            let mut it = keys.into_iter();
            let mut current_key = it.next().unwrap();
            let mut count = 1;
            for key in it {
                if key == current_key {
                    count += 1;
                } else {
                    if count >= min_required_score { results_for_shard.push((current_key, count)); }
                    current_key = key;
                    count = 1;
                }
            }
            if count >= min_required_score { results_for_shard.push((current_key, count)); }
            
            fs::remove_file(&path).ok();
            Some(results_for_shard)
        }).collect();
        
        let mut all_results: Vec<(u64, u32)> = aggregated_pairs.into_iter().flatten().collect();
        
        all_results.par_sort_unstable_by(|a, b| {
            let a_query = a.0 >> 32;
            let b_query = b.0 >> 32;
            match a_query.cmp(&b_query) {
                std::cmp::Ordering::Equal => b.1.cmp(&a.1),
                other => other,
            }
        });

        println!("Aggregation finished in {} seconds.", aggregation_start_time.elapsed().as_secs());
        
        let output_lines: Vec<String> = (0..num_queries).into_par_iter().map(|query_idx| {
            let query_gid = query_idx as u64;
            let mut line = String::with_capacity(256);
            write!(&mut line, "query_{}\t", query_gid).unwrap();

            let start_idx = all_results.partition_point(|(key, _)| (*key >> 32) < query_gid);
            let end_idx = all_results.partition_point(|(key, _)| (*key >> 32) <= query_gid);
            let results_slice = &all_results[start_idx..end_idx];

            if !results_slice.is_empty() {
                if is_matrix {
                    let mut full_row = vec![0.0; num_refs];
                    for &(key, count) in results_slice {
                        let ref_gid = key as Gid;
                        if (ref_gid as usize) < num_refs {
                             full_row[ref_gid as usize] = count as f64 / self.s as f64;
                        }
                    }
                    line.push_str(&full_row.iter().map(|&sim| if sim >= threshold { format!("{:.4}", sim) } else { "0.0000".to_string() }).collect::<Vec<String>>().join(","));
                } else {
                    line.push_str(&results_slice.iter()
                        .map(|(key, count)| format!("{}:{:.4}", *key as Gid, *count as f64 / self.s as f64))
                        .collect::<Vec<String>>()
                        .join(","));
                }
            } else if is_matrix {
                let full_row = vec!["0.0000"; num_refs];
                line.push_str(&full_row.join(","));
            }
            line
        }).collect();

        for line in output_lines { writeln!(writer, "{}", line).unwrap(); }
        
        println!("\nAll-vs-all comparison complete. Total time: {} seconds.", start_time.elapsed().as_secs());
    }
   
    fn create_writer<'a>(&self, output_file: &'a str, zstd_level: i32) -> Box<dyn Write + 'a> {
        let file = File::create(output_file).expect("Cannot create output file");
        if zstd_level > 0 { Box::new(ZstdEncoder::new(file, zstd_level).unwrap().auto_finish()) } else { Box::new(BufWriter::new(file)) }
    }

    fn path_for_reading(&self) -> &str {
        let IndexData::OnDisk { path, .. } = &self.data;
        path
    }
    
    fn toc_for_reading(&self) -> &[(u64, u64)] {
        let IndexData::OnDisk { toc, .. } = &self.data;
        toc
    }
}

fn read_pos_from_file(file: &mut File, toc: &[(u64, u64)], pos: usize, fp_range: u64) -> Vec<(u32, Vec<u8>)> {
    let (offset, len) = toc[pos];
    if len == 0 { return vec![(0, Vec::new()); fp_range as usize]; }
    
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

fn decode_gid_list(count: u32, bytes: &[u8], is_compressed: bool) -> Vec<Gid> {
    if count == 0 { return Vec::new(); }
    if is_compressed {
        let mut gids = vec![0 as Gid; count as usize];
        stream_vbyte::decode::decode::<stream_vbyte::scalar::Scalar>(bytes, count as usize, &mut gids);
        for j in 1..gids.len() { gids[j] = gids[j].wrapping_add(gids[j - 1]); }
        gids
    } else {
        let mut gids = Vec::with_capacity(count as usize);
        let mut cursor = io::Cursor::new(bytes);
        for _ in 0..count { gids.push(cursor.read_u32::<LittleEndian>().unwrap()); }
        gids
    }
}