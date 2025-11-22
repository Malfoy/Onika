# Onika

Onika is an efficient Rust MinHash sketcher and similarity engine built around compressed inverted indexes.


**Please cite:** Onika preprint – https://www.biorxiv.org/content/10.1101/2025.11.21.689685v1



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
   ```./target/release/Onika sketch --input-fof data/fof.txt --k_size 31 --s_size 10 --w_size 16  --reorder-similarity -o index.bin
   ```
3) Compare  the built index against itself:
   ``` ./target/release/Onika compare --ref-sketch index.bin --query-sketch  index.bin -o out.tsv.zst
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

- `--reorder-similarity`: greedily reorder sketches by self-similarity.

**`compare`**
- Reference inputs: `--ref-sketch <FILE>` | `--ref-fof <FILE>` | `--ref-fasta <FILE>`.
- Query inputs: `--query-sketch <FILE>` | `--query-fof <FILE>` | `--query-fasta <FILE>`.
- `-o, --output <FILE>`: destination for similarity results (Zstandard-compressed unless `--zstd-level 0`).
- `--threshold <FLOAT>`: similarity cutoff. Below this value hits are pruned in sparse mode; also controls how many bucket matches are required.
- `--matrix`: emit a full similarity matrix row per query (comma-separated floats). Without it, output is sparse.
- `--prob-threshold-probability <FLOAT>`: probabilistic pruning tail probability; set to `0` to disable the heuristic (default 0.001).


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
