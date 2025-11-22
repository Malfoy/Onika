# Onika

Onika is an open-source Rust MinHash sketcher and similarity engine built around compressed inverted indexes. It delivers output-sensitive all-vs-all comparisons whose runtime scales with matching sketch positions, adds deterministic pruning for pairs that cannot reach a target Jaccard threshold, and provides a probabilistic early-stop variant with tunable false-rejection control. A similarity-aware sketch reordering pass further shrinks indexes and improves locality, especially in redundant collections.

**Please cite:** Onika preprint – https://www.biorxiv.org/content/10.1101/2025.11.21.689685v1

## Highlights
- Compressed inverted posting lists over MinHash fingerprints; output-sensitive all-vs-all comparison.
- Early pruning: exact threshold-based elimination plus probabilistic pruning with controllable tail probability.
- Similarity-aware sketch reordering to tighten locality and index size for redundant datasets.
- Ingests plain, `.gz`, or `.zst` FASTA/FASTQ; accepts single files or file-of-files (one path per line).
- Stream-vbyte delta encoding plus optional Zstandard compression for the index and comparison output.
- Parallelized with Rayon; jemalloc is used as the global allocator.
- Experimental: collision-aware fingerprint scaling exists behind a sketch-mode flag and may change.

## Requirements
- Rust toolchain (edition 2021; `cargo` from stable is sufficient).
- Typical build: `rustup default stable` and `cargo --version`.

## Build
```bash
git clone https://github.com/Malfoy/Onika.git
cd Onika
cargo build --release
# Binary is at target/release/Onika
target/release/Onika --help
```

## Input formats
- **Single FASTA/Q**: `--input-fasta some_reads.fq[.gz|.zst]` treats each record as one document.
- **File-of-files (FOF)**: `--input-fof samples.fof` where `samples.fof` lists one FASTA/Q path per line. The order defines the genome IDs used in outputs (`0`-indexed).
- Both modes accept gzipped or zstd-compressed input transparently.

## Quick start
1) Inspect the bundled sample FOF (two genomes in `data/`):
   ```bash
   cat data/fof.txt
   ```
2) Build a sketch index from the sample FOF:
   ```./target/release/Onika sketch \
     --input-fof data/fof.txt \
     --k_size 31 --s_size 10 --w_size 16 \
     --zstd-level 3 \
     --reorder-similarity \
     -o sample.index.bin
   ```
3) Compare the same FOF against the built index (self-comparison demo):
   ```./target/release/Onika compare \     
     --ref-sketch sample.index.bin \
     --query-sketch  sample.index.bin \
     --threshold 0.05 \
     --prob-threshold-probability 0.001 \
     -o sample_ref_vs_queries.tsv.zst
   ```
   By default the output is Zstandard-compressed. Use `--zstd-level 0` to write plain text.

## CLI reference
**Common**
- `--threads <INT>`: override the auto-detected thread count.
- `--zstd-level <LEVEL>`: compression level for the index (sketch) or comparison output; `0` disables compression.

**`sketch`**
- `--input-fof <FILE>` | `--input-fasta <FILE>`: input source (one required).
- `-o, --output <FILE>`: where to write the sketch index (`.bin` is typical).
- `-k, --k-size <INT>`: k-mer size (default 31).
- `-s, --s-size <INT>`: log2 of sketch partitions; sketch has `2^s` buckets (default 10 → 1024 partitions).
- `-w, --w-size <INT>`: fingerprint width in bits; values are reduced into a `2^w` range (default 16).
- `-e, --e-size <INT>`: expected genome length (defaults to 5,000,000; can be auto-estimated from early records when needed).
- `--sketch-mode [default|perfect]`: choose fingerprinting strategy; `perfect` is experimental.
- `--no-compress`: skip stream-vbyte compression inside the index (larger files, faster build).
- `--reorder-similarity`: greedily reorder sketches by self-similarity; tune with `--reorder-threshold` (min similarity considered) and `--reorder-sample` (number of partitions sampled, `0` = full sketch).

**`compare`**
- Reference inputs: `--ref-sketch <FILE>` | `--ref-fof <FILE>` | `--ref-fasta <FILE>`.
- Query inputs: `--query-sketch <FILE>` | `--query-fof <FILE>` | `--query-fasta <FILE>`.
- `-o, --output <FILE>`: destination for similarity results (Zstandard-compressed unless `--zstd-level 0`).
- `--threshold <FLOAT>`: similarity cutoff. Below this value hits are pruned in sparse mode; also controls how many bucket matches are required.
- `--matrix`: emit a full similarity matrix row per query (comma-separated floats). Without it, output is sparse.
- `--prob-threshold-probability <FLOAT>`: probabilistic pruning tail probability; set to `0` to disable the heuristic (default 0.001).
- `--no-compress`: write uncompressed shard data when building indices on-the-fly (debugging).
- `--sketch-mode`, `--k-size`, `--s-size`, `--w-size`, `--e-size`, `--zstd-level`: must match how the sketches were built when constructing indices on-the-fly.

## Output formats
- **Index (`.bin`)**: binary file containing the sketch table and metadata; compressed per-position blocks with stream-vbyte deltas and optional Zstandard. Load with `--ref-sketch` or `--query-sketch`.
- **Sparse comparison (default)**: each line is `query_<qid>\tref_gid:score[,ref_gid:score...]`, sorted by score descending. IDs are zero-based in the order the inputs were read. `score` is the fraction of matching buckets (`matches / sketch_size`).
- **Matrix mode**: each line is `query_<qid>\t` followed by comma-separated similarity scores for every reference (values below the threshold are written as `0.0000`).

## Tips
- Larger `s-size` adds resolution but increases memory and runtime; `w-size` controls fingerprint collision rate in the index.
- `--reorder-similarity` can improve spatial locality of sketches before compression; leave it off for the fastest build.
- Outputs and indexes are Zstandard-compressed by default; set `--zstd-level 0` if you need plain text for downstream tools.

## License
This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0). See `LICENSE` for details.
