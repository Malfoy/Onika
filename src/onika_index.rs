use crate::utils::*;
use ahash::AHashMap;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use f128::f128;
use flate2::read::GzDecoder;
use io::BufWriter;
use nthash::NtHashIterator;
use num_traits::{Float, ToPrimitive};
use parking_lot::Mutex;
use probability::distribution::{Binomial, Discrete};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rlimit;
use std::collections::hash_map::Entry;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::mem;
use std::path::Path;
use std::sync::{
    atomic::{AtomicU32, AtomicU64, AtomicUsize, Ordering},
    Arc,
};
use tempfile::NamedTempFile;
use zstd::stream::write::Encoder as ZstdEncoder;

type Gid = u32;

const LIST_TO_MAP_THRESHOLD: usize = 128;

enum QueryData {
    List(Vec<(u32, u32)>),
    Map(AHashMap<u32, u32>),
}

struct QueryAccumulator {
    data: QueryData,
}

impl QueryAccumulator {
    fn new() -> Self {
        Self {
            data: QueryData::List(Vec::new()),
        }
    }

    fn promote_to_map(&mut self) {
        if let QueryData::List(list) = mem::replace(&mut self.data, QueryData::Map(AHashMap::new()))
        {
            if let QueryData::Map(ref mut map) = self.data {
                map.reserve(list.len());
                for (ref_gid, stats) in list {
                    map.insert(ref_gid, stats);
                }
            }
        }
    }

    fn apply_hit(
        &mut self,
        ref_gid: u32,
        enable_prob_heuristic: bool,
        min_required_score: u32,
        processed_positions: u64,
        future_after_current: u64,
        target_similarity: f64,
        prob_cutoff: f64,
    ) {
        match &mut self.data {
            QueryData::List(list) => {
                let mut promote = false;
                match list.binary_search_by_key(&ref_gid, |(gid, _)| *gid) {
                    Ok(idx) => {
                        if min_required_score == 0 {
                            list[idx].1 += 1;
                        } else {
                            list[idx].1 += 1;
                            let count = list[idx].1;
                            if !Index::probability_allows_threshold(
                                enable_prob_heuristic,
                                target_similarity,
                                prob_cutoff,
                                min_required_score,
                                count,
                                processed_positions,
                                future_after_current,
                            ) {
                                list.remove(idx);
                            }
                        }
                    }
                    Err(insert_idx) => {
                        if min_required_score == 0
                            || Index::probability_allows_threshold(
                                enable_prob_heuristic,
                                target_similarity,
                                prob_cutoff,
                                min_required_score,
                                1,
                                processed_positions,
                                future_after_current,
                            )
                        {
                            list.insert(insert_idx, (ref_gid, 1));
                            if list.len() > LIST_TO_MAP_THRESHOLD {
                                promote = true;
                            }
                        }
                    }
                }
                if promote {
                    self.promote_to_map();
                }
            }
            QueryData::Map(map) => match map.entry(ref_gid) {
                Entry::Occupied(mut occ) => {
                    if min_required_score == 0 {
                        *occ.get_mut() += 1;
                        return;
                    }
                    let val = occ.get_mut();
                    *val += 1;
                    let count = *val;
                    if !Index::probability_allows_threshold(
                        enable_prob_heuristic,
                        target_similarity,
                        prob_cutoff,
                        min_required_score,
                        count,
                        processed_positions,
                        future_after_current,
                    ) {
                        occ.remove_entry();
                    }
                }
                Entry::Vacant(vac) => {
                    if min_required_score == 0
                        || Index::probability_allows_threshold(
                            enable_prob_heuristic,
                            target_similarity,
                            prob_cutoff,
                            min_required_score,
                            1,
                            processed_positions,
                            future_after_current,
                        )
                    {
                        vac.insert(1);
                    }
                }
            },
        }
    }

    fn for_each<F: FnMut(u32, u32)>(&self, mut f: F) {
        match &self.data {
            QueryData::List(list) => {
                for &(gid, count) in list.iter() {
                    f(gid, count);
                }
            }
            QueryData::Map(map) => {
                for (&gid, &count) in map.iter() {
                    f(gid, count);
                }
            }
        }
    }

    fn collect_filtered(&self, min_score: u32) -> Vec<(u32, u32)> {
        match &self.data {
            QueryData::List(list) => list
                .iter()
                .filter(|(_, count)| *count >= min_score)
                .map(|&(gid, count)| (gid, count))
                .collect(),
            QueryData::Map(map) => map
                .iter()
                .filter(|(_, &count)| count >= min_score)
                .map(|(&gid, &count)| (gid, count))
                .collect(),
        }
    }
}

struct PendingRecord {
    count: u32,
    processed_positions: u64,
    future_after_current: u64,
}

impl Index {
    fn set_partition_limit(&mut self, limit: u64) {
        let IndexData::OnDisk { toc, .. } = &mut self.data;
        if limit == 0 || limit >= toc.len() as u64 {
            self.sketch_size = self.s;
        } else {
            toc.truncate(limit as usize);
            self.sketch_size = limit;
            self.s = limit;
        }
        self.s = self.sketch_size;
    }
}

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
    fn clone_for_reorder(&self) -> Self {
        let intermediate_data = self
            .intermediate_data
            .iter()
            .map(|val| AtomicU32::new(val.load(Ordering::Relaxed)))
            .collect();
        Self {
            ls: self.ls,
            w: self.w,
            k: self.k,
            e: self.e,
            s: self.s,
            sketch_size: self.sketch_size,
            fingerprint_range: self.fingerprint_range,
            mi: self.mi,
            intermediate_data,
            num_docs: self.num_docs,
            genome_numbers: AtomicU64::new(self.genome_numbers.load(Ordering::Relaxed)),
            empty_sketches_count: AtomicU64::new(self.empty_sketches_count.load(Ordering::Relaxed)),
            sketch_mode: self.sketch_mode,
        }
    }

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

    pub fn into_final_index(
        self,
        output_path: &str,
        compress: bool,
        zstd_level: i32,
    ) -> io::Result<Index> {
        let output_dir = Path::new(output_path)
            .parent()
            .unwrap_or_else(|| Path::new("."));
        let temp_file = tempfile::Builder::new()
            .prefix(".tmp-index-")
            .tempfile_in(output_dir)?;
        let path_str = temp_file.path().to_str().unwrap().to_string();

        let position_buffers: Vec<Vec<u8>> = self
            .intermediate_data
            .par_chunks(self.num_docs)
            .map(|fingerprint_chunk| {
                let mut inverted_index: Vec<Vec<Gid>> =
                    vec![Vec::new(); self.fingerprint_range as usize];
                for (doc_id, fp_atomic) in fingerprint_chunk.iter().enumerate() {
                    let fp = fp_atomic.load(Ordering::Relaxed);
                    if fp != self.mi {
                        inverted_index[fp as usize].push(doc_id as Gid);
                    }
                }

                let mut uncompressed_buffer = Vec::new();
                for mut gids in inverted_index {
                    let original_count = gids.len() as u32;
                    let gids_bytes = if !gids.is_empty() {
                        gids.sort_unstable();
                        let mut prev = gids[0];
                        for i in 1..gids.len() {
                            let current = gids[i];
                            gids[i] = current.wrapping_sub(prev);
                            prev = current;
                        }
                        if compress {
                            let mut encoded_data = vec![0; 5 * gids.len()];
                            let encoded_len = stream_vbyte::encode::encode::<
                                stream_vbyte::scalar::Scalar,
                            >(
                                &gids, &mut encoded_data
                            );
                            encoded_data.truncate(encoded_len);
                            encoded_data
                        } else {
                            let mut bytes = Vec::with_capacity(gids.len() * 4);
                            for delta in &gids {
                                bytes.write_u32::<LittleEndian>(*delta).unwrap();
                            }
                            bytes
                        }
                    } else {
                        Vec::new()
                    };
                    uncompressed_buffer
                        .write_u32::<LittleEndian>(original_count)
                        .unwrap();
                    uncompressed_buffer
                        .write_u32::<LittleEndian>(gids_bytes.len() as u32)
                        .unwrap();
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
            eprintln!(
                "\nWarning: {} input documents resulted in empty sketches and were not indexed.",
                empty_count
            );
        }
        println!("Index construction complete.");

        let (file_handle, path) = temp_file.into_parts();

        Ok(Index {
            _ls: self.ls,
            _w: self.w,
            _k: self.k,
            _e: self.e,
            s: self.s,
            sketch_size: self.sketch_size,
            fingerprint_range: self.fingerprint_range,
            _mi: self.mi as u64,
            data: IndexData::OnDisk {
                handle: Some(Arc::new(NamedTempFile::from_parts(file_handle, path))),
                path: path_str,
                toc,
            },
            genome_numbers: self.genome_numbers.load(Ordering::Relaxed),
            _sketch_mode: self.sketch_mode,
            is_compressed: compress,
        })
    }

    pub fn index_fasta_file(&mut self, filestr: &str) {
        let reader = open_compressed_file(filestr).expect("Could not open FASTA/FASTQ file");
        let mut fastx_reader =
            needletail::parse_fastx_reader(reader).expect("Invalid FASTA/FASTQ format");
        let mut sequences = Vec::new();
        while let Some(record) = fastx_reader.next() {
            let seqrec = record.expect("Invalid record in FASTA/FASTQ file");
            sequences.push(seqrec.seq().to_vec());
        }

        let num_docs = sequences.len();
        self.num_docs = num_docs;
        self.genome_numbers.store(num_docs as u64, Ordering::SeqCst);
        let total_size = self.s as usize * num_docs;
        self.intermediate_data
            .resize_with(total_size, || AtomicU32::new(self.mi));
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
        self.intermediate_data
            .resize_with(total_size, || AtomicU32::new(self.mi));
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

    pub fn reorder_by_similarity(
        &mut self,
        min_similarity: f64,
        sample_partitions: usize,
        threads: usize,
    ) -> io::Result<Vec<usize>> {
        let num_docs = self.num_docs;
        if num_docs <= 1 {
            return Ok((0..num_docs).collect());
        }
        let total_positions = self.s as usize;
        if total_positions == 0 {
            return Ok((0..num_docs).collect());
        }

        let min_similarity = min_similarity.clamp(0.0, 1.0);
        let temp_builder = self.clone_for_reorder();
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .unwrap();

        let order = pool.install(move || -> io::Result<Vec<usize>> {
            let temp_dir = tempfile::tempdir()?;
            let temp_index_path = temp_dir.path().join("reorder_index.bin");
            let mut temp_index = temp_builder.into_final_index(
                temp_index_path.to_str().expect("Invalid temp index path"),
                false,
                0,
            )?;

            if sample_partitions > 0 && sample_partitions < total_positions {
                temp_index.set_partition_limit(sample_partitions as u64);
            }

            let temp_output_path = temp_dir.path().join("reorder_output.tsv");
            temp_index.all_vs_self_comparison(
                min_similarity,
                false,
                0,
                temp_output_path.to_str().expect("Invalid temp output path"),
                threads.max(1),
                false,
                1.0 / 1_000_000.0,
            );
            drop(temp_index);

            let reader = BufReader::new(File::open(&temp_output_path)?);
            let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); num_docs];
            for line_res in reader.lines() {
                let line = line_res?;
                if line.is_empty() {
                    continue;
                }
                let mut parts = line.split('\t');
                let query_part = parts.next().unwrap_or_default();
                let entries_part = parts.next().unwrap_or_default();
                let query_gid = if let Some(idx_str) = query_part.strip_prefix("query_") {
                    idx_str.parse::<usize>().unwrap_or(0)
                } else {
                    continue;
                };
                if entries_part.is_empty() {
                    continue;
                }
                for entry in entries_part.split(',') {
                    if entry.is_empty() {
                        continue;
                    }
                    let mut kv = entry.split(':');
                    let ref_gid = kv
                        .next()
                        .and_then(|v| v.parse::<usize>().ok())
                        .unwrap_or(usize::MAX);
                    let sim = kv.next().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
                    if ref_gid == usize::MAX || ref_gid == query_gid || sim < min_similarity {
                        continue;
                    }
                    adjacency[query_gid].push((ref_gid, sim));
                    adjacency[ref_gid].push((query_gid, sim));
                }
            }

            let mut total_weight: Vec<f64> = adjacency
                .iter()
                .map(|neighbors| neighbors.iter().map(|(_, sim)| *sim).sum())
                .collect();
            let mut visited = vec![false; num_docs];
            let mut order = Vec::with_capacity(num_docs);

            let mut current = total_weight
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            for _ in 0..num_docs {
                order.push(current);
                visited[current] = true;
                total_weight[current] = f64::MIN;

                if order.len() == num_docs {
                    break;
                }

                let mut best_choice: Option<(usize, f64)> = None;
                for &(neighbor, sim) in &adjacency[current] {
                    if visited[neighbor] {
                        continue;
                    }
                    match best_choice {
                        Some((best_idx, best_sim)) => {
                            if sim > best_sim || (sim == best_sim && neighbor < best_idx) {
                                best_choice = Some((neighbor, sim));
                            }
                        }
                        None => best_choice = Some((neighbor, sim)),
                    }
                }

                current = if let Some((neighbor, _)) = best_choice {
                    neighbor
                } else {
                    total_weight
                        .iter()
                        .enumerate()
                        .filter(|(idx, _)| !visited[*idx])
                        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                        .map(|(idx, _)| idx)
                        .unwrap_or_else(|| (0..num_docs).find(|idx| !visited[*idx]).unwrap_or(0))
                };
            }

            Ok(order)
        })?;

        let raw_data = mem::take(&mut self.intermediate_data);
        let dense: Vec<u32> = raw_data.into_iter().map(|atom| atom.into_inner()).collect();
        let mut reordered = Vec::with_capacity(dense.len());
        for pos in 0..total_positions {
            let base = pos * num_docs;
            for &doc_idx in &order {
                reordered.push(dense[base + doc_idx]);
            }
        }
        self.intermediate_data = reordered
            .into_iter()
            .map(|value| AtomicU32::new(value))
            .collect();

        Ok(order)
    }

    fn compute_sketch(&self, input: SketchInput, sketch: &mut Vec<u64>) {
        sketch.fill(self.mi as u64);
        let mut empty_cell = self.s;
        let is_valid_dna =
            |c: &u8| matches!(*c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
        let mut update_sketch = |sequence: &[u8]| {
            for subseq in sequence.split(|c| !is_valid_dna(c)) {
                if subseq.len() >= self.k as usize {
                    let h_iter = NtHashIterator::new(subseq, self.k as usize).unwrap();
                    for canon_hash in h_iter {
                        let permuted = unrevhash64(canon_hash);
                        let bucket_id = (permuted >> (64 - self.ls)) as usize;
                        if bucket_id < sketch.len() {
                            if sketch[bucket_id] == self.mi as u64 {
                                empty_cell -= 1;
                                sketch[bucket_id] = canon_hash;
                            } else if sketch[bucket_id] > canon_hash {
                                sketch[bucket_id] = canon_hash;
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
            if empty_cell > 0 {
                self.sketch_densification(sketch, empty_cell as u32);
            }
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
            if changes.is_empty() && empty_cell > 0 {
                break;
            }

            for (pos, val) in changes {
                if sketch[pos] == self.mi as u64 {
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
        if hashed == self.mi as u64 {
            return self.mi as u64;
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
            toc.push((
                file.read_u64::<LittleEndian>()?,
                file.read_u64::<LittleEndian>()?,
            ));
        }

        Ok(Self {
            _ls: ls,
            _w: w,
            _k: k,
            _e: e,
            s,
            sketch_size: s,
            fingerprint_range,
            _mi: u64::MAX,
            data: IndexData::OnDisk {
                handle: None,
                path: filestr.to_string(),
                toc,
            },
            genome_numbers,
            _sketch_mode: sketch_mode,
            is_compressed,
        })
    }

    pub fn dump_index_disk(self, filestr: &str) -> io::Result<()> {
        let IndexData::OnDisk { handle, path, .. } = self.data;
        if let Some(temp_handle) = handle {
            if let Ok(temp_file) = Arc::try_unwrap(temp_handle) {
                temp_file.persist(filestr)?;
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to unwrap Arc for temp file",
                ));
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
        use_probabilistic_pruning: bool,
        prob_threshold_probability: f64,
    ) {
        let _ = shard_zstd_level;
        let _ = temp_dir_path;

        Self::run_pairwise_comparison(
            self,
            query_index,
            threshold,
            is_matrix,
            zstd_level,
            output_file,
            num_threads,
            use_probabilistic_pruning,
            prob_threshold_probability,
            false,
            "All-vs-all",
        );
    }

    pub fn all_vs_self_comparison(
        &self,
        threshold: f64,
        is_matrix: bool,
        zstd_level: i32,
        output_file: &str,
        num_threads: usize,
        use_probabilistic_pruning: bool,
        prob_threshold_probability: f64,
    ) {
        Self::run_pairwise_comparison(
            self,
            self,
            threshold,
            is_matrix,
            zstd_level,
            output_file,
            num_threads,
            use_probabilistic_pruning,
            prob_threshold_probability,
            true,
            "All-vs-self",
        );
    }

    fn run_pairwise_comparison(
        ref_index: &Self,
        query_index: &Self,
        threshold: f64,
        is_matrix: bool,
        zstd_level: i32,
        output_file: &str,
        num_threads: usize,
        use_probabilistic_pruning: bool,
        prob_threshold_probability: f64,
        skip_identical_hits: bool,
        label: &str,
    ) {
        if let Ok((soft_limit, hard_limit)) = rlimit::getrlimit(rlimit::Resource::NOFILE) {
            if soft_limit < hard_limit {
                if let Err(e) = rlimit::setrlimit(rlimit::Resource::NOFILE, hard_limit, hard_limit)
                {
                    eprintln!(
                        "Warning: Failed to increase open file limit: {}. This may cause errors.",
                        e
                    );
                }
            }
        }

        let start_time = std::time::Instant::now();

        let num_queries = query_index.genome_numbers as usize;
        let num_refs = ref_index.genome_numbers as usize;

        let min_required_score = if threshold > 0.0 {
            (threshold * ref_index.s as f64).ceil() as u32
        } else {
            0
        };
        let probabilistic_active = use_probabilistic_pruning && min_required_score > 0;
        let prob_cutoff = prob_threshold_probability.clamp(0.0, 1.0);

        let mut writer = ref_index.create_writer(output_file, zstd_level);

        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let shared_counts: Vec<Mutex<QueryAccumulator>> = (0..num_queries)
            .map(|_| Mutex::new(QueryAccumulator::new()))
            .collect();

        pool.install(|| {
            let total_positions = ref_index.sketch_size as u64;
            let num_pos = total_positions as usize;
            let chunks_per_thread = 4usize;
            let chunk_size = ((num_pos + (num_threads * chunks_per_thread) - 1)
                / (num_threads * chunks_per_thread))
                .max(1);
            let pos_iterator: Vec<usize> = (0..num_pos).collect();
            let pos_chunks: Vec<&[usize]> = pos_iterator.chunks(chunk_size).collect();

            let remaining_positions = AtomicUsize::new(num_pos);
            let enable_prob_heuristic = probabilistic_active;
            let shared_counts = &shared_counts;
            let remaining_positions = &remaining_positions;

            pos_chunks.into_par_iter().for_each(|pos_chunk| {
                let total_positions = total_positions;
                let enable_prob_heuristic = enable_prob_heuristic;
                let remaining_positions = remaining_positions;
                let mut ref_file = File::open(&ref_index.path_for_reading())
                    .expect("Cannot open reference index file in thread.");
                let mut query_file = File::open(&query_index.path_for_reading())
                    .expect("Cannot open query index file in thread.");

                let mut ref_gid_buf: Vec<Gid> = Vec::new();
                let mut query_gid_buf: Vec<Gid> = Vec::new();
                let mut pending_updates: AHashMap<u32, AHashMap<u32, PendingRecord>> =
                    AHashMap::new();

                for &pos in pos_chunk {
                    let future_positions = remaining_positions.load(Ordering::Relaxed) as u64;
                    let future_after_current = future_positions.saturating_sub(1);
                    let processed_positions = total_positions.saturating_sub(future_after_current);
                    let ref_pos_index = read_pos_from_file(
                        &mut ref_file,
                        &ref_index.toc_for_reading(),
                        pos,
                        ref_index.fingerprint_range,
                    );
                    let query_pos_index = read_pos_from_file(
                        &mut query_file,
                        &query_index.toc_for_reading(),
                        pos,
                        query_index.fingerprint_range,
                    );

                    let fp_range = ref_index.fingerprint_range as usize;
                    for fp in 0..fp_range {
                        let (ref_count, ref_bytes) = ref_pos_index.entry(fp);
                        if ref_count == 0 {
                            continue;
                        }
                        decode_gid_list_into(
                            &mut ref_gid_buf,
                            ref_count,
                            ref_bytes,
                            ref_index.is_compressed,
                        );
                        if ref_gid_buf.is_empty() {
                            continue;
                        }

                        let (query_count, query_bytes) = query_pos_index.entry(fp);
                        if query_count == 0 {
                            continue;
                        }
                        decode_gid_list_into(
                            &mut query_gid_buf,
                            query_count,
                            query_bytes,
                            query_index.is_compressed,
                        );
                        if query_gid_buf.is_empty() {
                            continue;
                        }

                        for &query_gid in &query_gid_buf {
                            if let Some(mutex) = shared_counts.get(query_gid as usize) {
                                if let Some(mut global_acc) = mutex.try_lock() {
                                    for &ref_gid in &ref_gid_buf {
                                        global_acc.apply_hit(
                                            ref_gid,
                                            enable_prob_heuristic,
                                            min_required_score,
                                            processed_positions,
                                            future_after_current,
                                            threshold,
                                            prob_cutoff,
                                        );
                                    }
                                } else {
                                    let query_pending = pending_updates
                                        .entry(query_gid)
                                        .or_insert_with(AHashMap::new);
                                    for &ref_gid in &ref_gid_buf {
                                        let entry =
                                            query_pending.entry(ref_gid).or_insert(PendingRecord {
                                                count: 0,
                                                processed_positions,
                                                future_after_current,
                                            });
                                        entry.count = entry.count.saturating_add(1);
                                        entry.processed_positions = processed_positions;
                                        entry.future_after_current = future_after_current;
                                    }
                                }
                            }
                        }

                        Self::flush_pending(
                            shared_counts,
                            &mut pending_updates,
                            enable_prob_heuristic,
                            min_required_score,
                            threshold,
                            prob_cutoff,
                        );
                    }

                    remaining_positions.fetch_sub(1, Ordering::Relaxed);
                }

                Self::flush_pending(
                    shared_counts,
                    &mut pending_updates,
                    enable_prob_heuristic,
                    min_required_score,
                    threshold,
                    prob_cutoff,
                );
            })
        });

        let line_results: Vec<(String, usize)> = (0..num_queries)
            .into_par_iter()
            .map(|query_idx| {
                let query_gid = query_idx as u32;
                let mut line = String::with_capacity(if is_matrix { num_refs * 8 } else { 512 });
                write!(&mut line, "query_{}\t", query_gid).unwrap();
                let mut scores_count = 0usize;

                if is_matrix {
                    let mut full_row = vec![0.0; num_refs];
                    {
                        let query_acc = shared_counts[query_idx].lock();
                        query_acc.for_each(|ref_gid, count| {
                            if count >= min_required_score {
                                let ref_idx = ref_gid as usize;
                                if ref_idx < num_refs {
                                    full_row[ref_idx] = count as f64 / ref_index.s as f64;
                                }
                            }
                        });
                    }

                    let mut first = true;
                    for sim in full_row {
                        if first {
                            first = false;
                        } else {
                            line.push(',');
                        }
                        if sim >= threshold {
                            scores_count += 1;
                            write!(&mut line, "{:.4}", sim).unwrap();
                        } else {
                            line.push_str("0.0000");
                        }
                    }
                } else {
                    let mut filtered: Vec<(u32, u32)> = {
                        let query_acc = shared_counts[query_idx].lock();
                        query_acc.collect_filtered(min_required_score)
                    };
                    filtered.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

                    let mut first = true;
                    for (ref_gid, count) in filtered {
                        if skip_identical_hits && ref_gid == query_gid {
                            continue;
                        }
                        if first {
                            first = false;
                        } else {
                            line.push(',');
                        }
                        let sim = count as f64 / ref_index.s as f64;
                        write!(&mut line, "{}:{:.4}", ref_gid, sim).unwrap();
                        scores_count += 1;
                    }
                }

                (line, scores_count)
            })
            .collect();

        for (line, _) in &line_results {
            writeln!(writer, "{}", line).unwrap();
        }

        let total_scores_above_threshold: usize =
            line_results.iter().map(|(_, count)| *count).sum();

        println!(
            "Total scores above threshold written to output: {}",
            total_scores_above_threshold
        );

        println!(
            "\n{} comparison complete. Total time: {} seconds.",
            label,
            start_time.elapsed().as_secs(),
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

    fn path_for_reading(&self) -> &str {
        let IndexData::OnDisk { path, .. } = &self.data;
        path
    }

    fn toc_for_reading(&self) -> &[(u64, u64)] {
        let IndexData::OnDisk { toc, .. } = &self.data;
        toc
    }

    fn flush_pending(
        shared_counts: &[Mutex<QueryAccumulator>],
        pending: &mut AHashMap<u32, AHashMap<u32, PendingRecord>>,
        enable_prob_heuristic: bool,
        min_required_score: u32,
        similarity_threshold: f64,
        prob_cutoff: f64,
    ) {
        if pending.is_empty() {
            return;
        }
        let query_ids: Vec<u32> = pending.keys().copied().collect();
        for query_gid in query_ids {
            if let Some(mutex) = shared_counts.get(query_gid as usize) {
                if let Some(mut global_acc) = mutex.try_lock() {
                    if let Some(records) = pending.remove(&query_gid) {
                        for (ref_gid, record) in records {
                            for _ in 0..record.count {
                                global_acc.apply_hit(
                                    ref_gid,
                                    enable_prob_heuristic,
                                    min_required_score,
                                    record.processed_positions,
                                    record.future_after_current,
                                    similarity_threshold,
                                    prob_cutoff,
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[inline]
    fn probability_allows_threshold(
        enable: bool,
        similarity_threshold: f64,
        prob_cutoff: f64,
        min_required: u32,
        count: u32,
        processed: u64,
        remaining: u64,
    ) -> bool {
        if !enable || min_required == 0 || prob_cutoff <= 0.0 {
            return true;
        }
        let prob_cutoff = prob_cutoff.min(1.0);
        if count >= min_required {
            return true;
        }

        let delta = min_required.saturating_sub(count);
        if delta as u64 > remaining {
            return false;
        }

        if !similarity_threshold.is_finite() {
            return true;
        }

        if processed == 0 {
            return true;
        }

        if processed < count as u64 {
            return true;
        }

        let processed_f = processed as f64;
        if processed_f <= 0.0 {
            return true;
        }

        let observed_ratio = (count as f64) / processed_f;
        if observed_ratio >= similarity_threshold {
            return true;
        }

        let threshold_prob = similarity_threshold.clamp(1e-12, 1.0 - 1e-12);

        let trials = match usize::try_from(processed) {
            Ok(v) if v > 0 => v,
            _ => return true,
        };

        let binomial = if threshold_prob <= 0.5 {
            Binomial::new(trials, threshold_prob)
        } else {
            Binomial::with_failure(trials, 1.0 - threshold_prob)
        };

        let probability = binomial.mass(count as usize);
        if !probability.is_finite() {
            return true;
        }

        probability >= prob_cutoff
    }
}

struct PositionEntry {
    count: u32,
    start: usize,
    end: usize,
}

struct PositionIndex {
    entries: Vec<PositionEntry>,
    buffer: Vec<u8>,
}

impl PositionIndex {
    #[inline]
    fn entry(&self, idx: usize) -> (u32, &[u8]) {
        let entry = &self.entries[idx];
        (entry.count, &self.buffer[entry.start..entry.end])
    }
}

fn read_pos_from_file(
    file: &mut File,
    toc: &[(u64, u64)],
    pos: usize,
    fp_range: u64,
) -> PositionIndex {
    let (offset, len) = toc[pos];
    if len == 0 {
        return PositionIndex {
            entries: (0..fp_range as usize)
                .map(|_| PositionEntry {
                    count: 0,
                    start: 0,
                    end: 0,
                })
                .collect(),
            buffer: Vec::new(),
        };
    }

    let mut compressed_buffer = vec![0; len as usize];
    file.seek(SeekFrom::Start(offset)).unwrap();
    file.read_exact(&mut compressed_buffer).unwrap();

    let buffer = zstd::decode_all(&compressed_buffer[..]).unwrap();
    let mut cursor = io::Cursor::new(&buffer);
    let mut entries = Vec::with_capacity(fp_range as usize);

    for _ in 0..fp_range as usize {
        let original_count = cursor.read_u32::<LittleEndian>().unwrap();
        let gids_len = cursor.read_u32::<LittleEndian>().unwrap() as usize;
        let start = cursor.position() as usize;
        let end = start + gids_len;
        cursor.set_position(end as u64);
        entries.push(PositionEntry {
            count: original_count,
            start,
            end,
        });
    }

    PositionIndex { entries, buffer }
}

fn decode_gid_list_into(output: &mut Vec<Gid>, count: u32, bytes: &[u8], is_compressed: bool) {
    output.clear();
    if count == 0 {
        return;
    }

    if is_compressed {
        output.resize(count as usize, 0);
        stream_vbyte::decode::decode::<stream_vbyte::scalar::Scalar>(bytes, count as usize, output);
        for idx in 1..output.len() {
            let prev = output[idx - 1];
            output[idx] = output[idx].wrapping_add(prev);
        }
    } else {
        output.reserve(count as usize);
        let mut cursor = io::Cursor::new(bytes);
        let mut prev = 0u32;
        for idx in 0..count {
            let delta = cursor.read_u32::<LittleEndian>().unwrap();
            if idx == 0 {
                prev = delta;
                output.push(delta);
            } else {
                prev = prev.wrapping_add(delta);
                output.push(prev);
            }
        }
    }
}
