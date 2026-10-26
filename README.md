# dart-nf

A Nextflow pipeline for [DART](https://github.com/genomewalker/dart) (Damage-Aware Read
Translation). For each sample it measures the library damage profile, predicts protein-coding
ORFs with decoy peptides, searches real and decoy proteins with MMseqs2, fits the decoy null
model, and runs EM reassignment with damage annotation.

## Requirements

- Nextflow 25.10 or later, Java 17 or later
- DART v2.0.1 or later
- MMseqs2 15.6f452
- Python 3.8 or later

## Test

```bash
nextflow run genomewalker/dart-nf -profile test \
  --dart /path/to/dart --mmseqs /path/to/mmseqs --outdir results
```

The test uses the DART tutorial data (50,000 synthetic reads and the 13,681 reference proteins
they come from) and finishes in about a minute on 8 cores. The last step compares the per-read
table with `test/expected/annotated.tsv.gz` and fails if any row differs. With
`-profile test,conda`, MMseqs2 is installed from Bioconda; `dart` is still taken from `--dart`.

A step-by-step walkthrough, from the test run to running a study on Slurm, is in
[docs/tutorial.md](docs/tutorial.md).

## Usage

Make a samplesheet with one row per sample:

```csv
sample,reads,profile_reads
S1,/data/S1.dedup.fq.gz,/data/S1.merged.fq.gz
```

`reads` is the FASTQ or FASTA that genes are predicted from. `profile_reads` is the FASTQ the
damage profile is measured on; leave it empty to use `reads`. For a deduplicated library, give
the merged reads before deduplication. Relative paths are resolved against the samplesheet's
directory.

```bash
nextflow run genomewalker/dart-nf -profile slurm --samplesheet samples.csv \
  --db /db/uniref90/uniref90 --dart /path/to/dart --mmseqs /path/to/mmseqs \
  --outdir results
```

`--db` takes a protein FASTA (plain or gzipped) or an MMseqs2 database prefix. Build a
database once with `mmseqs createdb ref.faa ref_db && mmseqs createindex ref_db tmp`.

The `slurm` profile submits every task to Slurm with 8 CPUs and 16 GB, and 16 CPUs and 32 GB
for the search. For a large database, set the queue, time and search memory in your own config
(`-c my.config`).

MMSEQS_SEARCH and DAMAGE_ANNOTATE run in node-local scratch. By default (`--scratch true`) that
is `$TMPDIR`, or `/tmp` if it is unset. Pass a path to use another directory, or `false` to run
in the work directory. `damage-annotate` writes its EM tape with `O_DIRECT`, so the scratch must
be a local disk, not tmpfs.

## Pipeline

| Process | Command | Outputs |
|---|---|---|
| PROFILE | `dart profile` | `<s>.profile.json`, `<s>.library.tdmg` |
| PREDICT | `dart predict --decoy-every 10 --fasta-aa --damage-index` | `<s>.predictions.gff`, `<s>.proteins.faa`, `<s>.proteins.decoy.faa`, `<s>.predictions.agd` |
| MMSEQS_SEARCH | `mmseqs easy-search`, proteins and decoys with the same options | `<s>.hits.tsv`, `<s>.decoy_hits.tsv` |
| HITS2EMI | `dart hits2emi --decoy-hits --decoy-fasta` | `<s>.hits.emi`, `<s>.hits.emi.null` |
| DAMAGE_ANNOTATE | `dart damage-annotate --gene-summary --protein-summary` | `<s>.annotated.tsv`, `<s>.gene_summary.tsv`, `<s>.protein_summary.tsv` |
| REPORT | `bin/dart_nf_report.py` | `<s>.report.md`, `<s>.versions.txt` |
| COMPARE_EXPECTED | `bin/compare_annotated.py` (with `--expected`) | `<s>.compare.txt` |

Each process also writes a log. All files are published to `<outdir>/<sample>/`.

The decoy search provides the null model: `predict` writes the reversed protein of every 10th
read, the decoys are searched with the same options as the real proteins, and `hits2emi` fits
one Gumbel per aligned-length bin to the decoy best hits. `damage-annotate` flags reads whose
per-bin target-decoy q-value exceeds `--null_fdr` and leaves them out of gene and protein
abundance.

Samples searched against the same database with the same options can share one null. Run the
pipeline on a few representative samples, pool their decoy hits with
`dart null -o study.null --db-label mydb a.decoy_hits.tsv b.decoy_hits.tsv`, and run the rest
with `--annotate_params '--em-null-model /abs/path/study.null --auto-calibrate-spurious'`.

## Outputs

- `<s>.annotated.tsv`: one row per read whose best alignment has a mismatch, with the assigned
  reference, ancient-read posterior, classification and EM flags.
- `<s>.gene_summary.tsv` and `<s>.protein_summary.tsv`: use `eff_reads` for abundance. The
  protein ancient call is `auth_binom_call`; trust it when `protein_auth_tier` is `detected`.
- `<s>.report.md`: read, ORF, decoy and hit counts, the null fit, the EM summary and the top
  proteins.

Column definitions: [Output Formats](https://github.com/genomewalker/dart/wiki/Output-Formats)
and [Interpreting Results](https://github.com/genomewalker/dart/wiki/Interpreting-Results).

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | | CSV with columns `sample,reads[,profile_reads]` |
| `--db` | | protein FASTA or MMseqs2 database prefix |
| `--outdir` | `results` | output directory |
| `--dart`, `--mmseqs` | `dart`, `mmseqs` | executables |
| `--vtml20` | `test/VTML20.out` | substitution matrix for the search |
| `--decoy_every` | 10 | one decoy peptide per N reads |
| `--domain` | `meta` | `predict --domain` |
| `--library_type` | `auto` | `profile --library-type` |
| `--search_params` | tutorial options | MMseqs2 options other than matrix, output format and threads |
| `--null_fdr` | 0.05 | per-bin target-decoy FDR for chance hits |
| `--null_p_max` | | use the fixed rule `null_p > P` instead of the FDR |
| `--annotate_params` | `--auto-calibrate-spurious` | extra `damage-annotate` options, e.g. `--keep-null-flagged` |
| `--scratch` | `true` | scratch for the search and annotation tasks |
| `--expected` | | per-read table to compare against (set by `-profile test`) |

See the DART [Quick Start Tutorial](https://github.com/genomewalker/dart/wiki/Quick-Start-Tutorial)
and [Command Reference](https://github.com/genomewalker/dart/wiki/Command-Reference).

## License

MIT. See [LICENSE](LICENSE).
