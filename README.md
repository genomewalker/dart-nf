# dart-nf

Status as of 2026-10-23: first version. Runs DART v2.0.0 (dart 1.1.0-116-g344d19f) end to end on
the bundled test data and reproduces the release's expected per-read table read for read.

The workflow runs, per sample: `dart profile`, `dart predict` with decoy peptides, two MMseqs2
searches (real and decoy proteins), `dart hits2emi` with the decoy null model, `dart damage-annotate`,
and a report step. Each step publishes its files to `<outdir>/<sample>/`.

## Quick start

Requirements: Nextflow 25.10 or later, Java 17 or later, MMseqs2 15.6f452, Python 3.8 or later,
and a `dart` binary built from the v2.0.0 release.

```bash
cd dart-nf
nextflow run . -profile test,local --outdir results
```

The `local` profile sets `--dart` and `--mmseqs` to the binaries on the Copenhagen cluster. Anywhere
else, pass them yourself:

```bash
nextflow run . -profile test --dart /path/to/dart --mmseqs /path/to/mmseqs --outdir results
```

With `-profile test,conda`, Nextflow installs MMseqs2 from bioconda for the search step; `dart` is
still taken from `--dart`.

The test run takes about one minute (55 s wall on one node with 8 cores). It ends with the `COMPARE_EXPECTED`
process, which compares the per-read table against `test/expected/annotated.tsv.gz` and fails the
run if any row differs. The expected table was produced at 8 threads; the test profile fixes
every process at 8 cpus for that reason.

The test data are the DART tutorial inputs: `test/reads.fq.gz` (50,000 synthetic reads from the
KapK-12-1-37 sediment sample, about 25% terminal C to T damage), `test/ref_proteins.faa.gz`
(the 13,681 reference proteins these reads come from) and `test/VTML20.out` (the substitution
matrix for short ancient peptides).

## What happens

| Process | Command | Published files |
|---|---|---|
| PROFILE | `dart profile -i <profile_reads> --damage-summary <s>.library.tdmg` | `<s>.profile.json`, `<s>.library.tdmg`, `<s>.profile.log` |
| PREDICT | `dart predict --library-profile <tdmg> --decoy-every 10 --fasta-aa --damage-index` | `<s>.predictions.gff`, `<s>.proteins.faa`, `<s>.proteins.decoy.faa`, `<s>.predictions.agd`, `<s>.predict.log` |
| MMSEQS_SEARCH | `mmseqs easy-search` twice, proteins then decoys, same options | `<s>.hits.tsv`, `<s>.decoy_hits.tsv`, `<s>.mmseqs.log` |
| HITS2EMI | `dart hits2emi --decoy-hits --decoy-fasta` | `<s>.hits.emi`, `<s>.hits.emi.null`, `<s>.hits2emi.log` |
| DAMAGE_ANNOTATE | `dart damage-annotate --gene-summary --protein-summary --auto-calibrate-spurious --em-tape-dir` | `<s>.annotated.tsv`, `<s>.gene_summary.tsv`, `<s>.protein_summary.tsv`, `<s>.damage_annotate.log` |
| REPORT | `bin/dart_nf_report.py` | `<s>.report.md`, `<s>.versions.txt` |
| COMPARE_EXPECTED | `bin/compare_annotated.py` (only with `--expected`) | `<s>.compare.txt` |

The profile step runs on every sample. The tutorial script in the dart repository instead passes
`--allow-self-profile` to `predict`. On the test reads the two give byte-identical GFF, protein,
decoy and AGD files, so the workflow uses the explicit profile and keeps one code path for test
and production. On a deduplicated library the profile must come from the merged reads before
deduplication; give that file in the `profile_reads` column of the samplesheet.

The decoy search is what makes the null model. `predict --decoy-every 10` writes the reversed
protein of every 10th read with a `DECOY_` prefix. The decoys are searched with the same MMseqs2
options as the real proteins, and `hits2emi` fits one Gumbel per aligned-length bin to the decoy
best hits and writes `<s>.hits.emi.null`. `damage-annotate` finds that file next to the EMI by
itself.

MMSEQS_SEARCH and DAMAGE_ANNOTATE run on node-local disk through Nextflow's `scratch` directive:
inputs are staged in and outputs copied back. `--scratch` controls it. The default `true` uses the
node's `$TMPDIR`, or `/tmp` when that is unset, so the quick start runs on any machine. The `local`
and `slurm` profiles set it to `/local/scratch`, the node-local disk of the Copenhagen cluster. A path
names another directory and `false` runs both tasks in the work directory. The EM tape in
`damage-annotate` opens its file with `O_DIRECT`, which tmpfs refuses and some network file systems
refuse, so the directory should be a real local disk. All three forms were run on the test data on
a compute node and matched the expected table, including `--scratch false` with the tape in the
work directory on the Isilon NFS.

## Reading the outputs

Column schemas are on the DART wiki:
[Output Formats](https://github.com/genomewalker/dart/wiki/Output-Formats) and
[Interpreting Results](https://github.com/genomewalker/dart/wiki/Interpreting-Results). This section
covers the columns needed to use the tables.

### `<s>.annotated.tsv`, one row per read with a mismatch

The table has one row per read whose best alignment has at least one mismatch. Reads that align
perfectly carry no substitution evidence and are not listed; the EM still assigned them, and
the gene and protein summaries count them.

- `query_id` is the ORF id from `predict`: the read name followed by `_<strand>_<frame>_<ordinal>`.
- `target_id` is the reference the EM assigned the read to; `best_hit` is 1 when that is also the
  bit-score-best target.
- `posterior` is the probability that the read is ancient, combining terminal damage (`p_read`),
  damage-consistent substitutions and identity. `is_damaged` is 1 when `posterior >= 0.7`.
- `classification` is `ancient_confident`, `ancient_likely`, `undetermined` or `modern_confident`.
  `tier` says which evidence the read has: `none`, `sites`, `terminal` or `convergent` (both).
- `em_keep` is 1 when the read's probability for its assigned target is at least
  `em_keep_threshold`. Rows with `em_keep = 0` are kept in the table so that a consumer can weight
  them; do not filter on this column without reading `em_margin` too.
- `em_ambiguous` is 1 when the assignment is not trusted as damage evidence: `em_reject_reason` is
  non-zero, or `em_margin` (log probability gap between the first and second target) is below the
  minimum, or `em_ref_neff` (effective number of candidate references) is above the maximum.
  Ambiguous reads keep their row but do not count toward the damage rate of their protein.
- `null_p` is the chance that a decoy peptide reaches the read's bit score in the same aligned-length
  bin. `NA` means the bin has no fitted null. A read with `null_p > --null_p_max` (default 0.05)
  gets reject bit 1.
- `em_reject_reason` is a bit mask: 1 = `null_p` above the threshold, 2 = low EM margin,
  4 = high effective reference count. 0 means kept. Values add, so 6 means low margin and high
  reference count.

### `<s>.gene_summary.tsv`, one row per predicted gene with enough reads

`n_reads`, `n_ancient`, `n_modern`, `n_undetermined` count the assigned reads by classification;
`ancient_frac` with `ci_low` and `ci_high` is the ancient share with its interval;
`mean_posterior` and `mean_identity` average the per-read values. `breadth`, `depth_evenness`,
`start_diversity` and `positional_score` describe how the reads cover the gene.
`--auto-calibrate-spurious` sets the thresholds on `positional_score` and `terminal_middle_ratio`
from the supported genes of the sample; the values it chose are in `<s>.damage_annotate.log` under
`Auto-calibration`.

### `<s>.protein_summary.tsv`, one row per reference protein with assigned reads

The ancient call for a protein is `auth_binom_call`, a binomial test of damaged reads against the
modern rate; `auth_binom_p_ancient` and `auth_binom_pvalue` are its probability and p-value.
Trust the call only when `protein_auth_tier` is `detected`; `advisory` means the library has a
damage rate but the test did not reach the detection threshold, `unauthenticatable` means the
library has no usable damage rate (`library_damage_state` 3), and `no_library` means no profile
was available. A `_filtered` suffix means the protein failed a mapping filter listed in
`filter_fail`. `p_protein_damaged` and `n_damaged` are the ungated damage estimate and count;
`protein_damage_confidence` is the gated one and is 0 when the protein is not authenticatable.

### `<s>.damage_annotate.log`

The `null reject:` lines, one per aligned-length bin, give the reads with a best hit in that bin,
how many were flagged by `null_p`, and the bit score below which a read in that bin is flagged.
`NA` means the bin holds too small a share of the decoys to flag anything at the chosen
`--null_p_max`. On a large protein database only the shortest bins (10-12 aa) usually flag.
The `EM reassignment` block gives the alignment and read counts, the iterations and the keep
filter; `Auto-calibration` gives the thresholds used for the gene summary.

### `<s>.report.md`

Counts of reads, ORFs, decoys and hits; the `dart damage-info` dump of the null sidecar; the
`null reject:` table; the EM summary; the posterior histogram and the classification, tier,
`em_keep`, `em_ambiguous` and `em_reject_reason` counts; the protein authentication tiers and the
top proteins; and the gene count. The versions of dart, MMseqs2 and Nextflow are at the top.

### What the test run shows

- 50,000 reads give 147,871 ORFs and 14,814 decoy peptides.
- MMseqs2 returns 45,714 alignments for 21,904 ORFs from 21,864 reads. The decoys get 3 alignments
  from one peptide: the reference set is only the 13,681 proteins the reads come from, too small for
  chance hits.
- Because one decoy query is not enough to fit a Gumbel, every null bin is `unfit`. `hits2emi` and
  `damage-annotate` both say so, and every `null reject:` line reads `flagged=0` with `NA`. The run
  is then the same as `--no-null`. A real protein database fills the bins.
- The per-read table has 8,500 rows. 4,364 have `posterior >= 0.5` and 1,252 have `is_damaged = 1`;
  `classification` is `ancient_confident` for 3,137, `undetermined` for 3,367 and `modern_confident`
  for 1,996. `em_reject_reason` is 0 for 7,057 reads, 2 for 638 and 6 for 805; `em_keep` is 0 for 106.
- The EM converges in 8 iterations over 45,714 alignments and keeps 21,458 of 21,904 ORFs.
- 346 genes reach the gene-summary thresholds; 218 of them have `ancient_frac >= 0.5`.
- The library profile from these synthetic reads has `library_damage_state` 3, so every protein is
  `unauthenticatable` and `auth_binom_call` is `undetermined_low_coverage` for all 7,075 proteins.
  The report ranks proteins by the ungated `p_protein_damaged` and says so.
- `<s>.compare.txt` ends with `RESULT: IDENTICAL`.

### The example run: a database large enough for chance hits

`-profile example` searches the same reads against 1,002,005 AlphaFold DB proteins (a uniform random
subset of the 241M AFDB members) and writes a decoy for every 2nd read. On a compute node it took
under a minute, most of it the 33 s search.

- 73,968 decoys give 515 alignments from 434 decoy queries (0.59% of the decoys submitted).
- The null fit is `1 own, 5 pooled, 0 unfit`: the 7-9 aa bin has 344 decoy queries and its own
  Gumbel (mu 29.0, beta 1.76 bits); the other bins share the pooled fit.
- The real reads get 32,517 alignments from 4,185 reads; 3,098 reads have a mismatch.
- No read is flagged. `null_p` is the chance that a random query gets a best hit in the read's
  length bin scoring at least as high, and on this database even the fullest bin (7-9 aa) holds only
  344 of 73,968 decoys, 0.47%. A bin can only flag when its share of decoys exceeds `--null-p-max`,
  so at 0.05 every `null reject:` line reads `NA`. That is the intended answer: against a 1M-protein
  database a hit of this length is rarely a chance hit. On the 241M-protein database the 10-12 aa bin
  holds 15.6% of the decoys and the same threshold flags its weaker hits.

## Running on real data

Make a samplesheet with one row per sample:

```csv
sample,reads,profile_reads
KapK-12-1-37,/data/KapK-12-1-37.dedup.fq.gz,/data/KapK-12-1-37.merged.fq.gz
```

`reads` is the FASTQ (or FASTA) the genes are predicted from. `profile_reads` is the FASTQ the
damage profile is measured on; leave it empty to profile `reads`. For a deduplicated library it
must be the merged reads before deduplication. Relative paths are resolved against the samplesheet's
directory.

```bash
nextflow run . -profile slurm --samplesheet samples.csv \
  --db /db/uniref90/uniref90 --dart /path/to/dart --mmseqs /path/to/mmseqs \
  --outdir results
```

`--db` takes a protein FASTA (plain or gzipped) or an MMseqs2 database prefix; with a prefix,
every file starting with it is staged (`createindex` output included). Build the database once
with `mmseqs createdb ref.faa ref_db && mmseqs createindex ref_db tmp`.

The `slurm` profile submits every task to Slurm with the resources in `nextflow.config`
(8 cpus and 16 GB by default, 16 cpus and 32 GB for the search). Adjust `process.queue`,
`time` and the memory of `MMSEQS_SEARCH` to the database. Nextflow's HyperQueue executor works the
same way; set `process.executor = 'hq'` in a profile of your own. The `slurm` profile sets `--scratch /local/scratch`;
point it at the local disk of your compute nodes, since the MMseqs2 temporary directory and the EM
tape both go there.

For a study with many samples, one null model per database and search configuration is enough.
Run this workflow on a few representative samples, pool their `<s>.decoy_hits.tsv` with
`dart null -o study.null --db-label mydb a/a.decoy_hits.tsv b/b.decoy_hits.tsv`, and pass
`--annotate_params '--em-null-model study.null --auto-calibrate-spurious'` for the rest; see the
dart README. `dart damage-info study.null` prints the fitted bins and the provenance.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--samplesheet` | none | CSV: `sample,reads[,profile_reads]` |
| `--db` | none | protein FASTA or MMseqs2 database prefix |
| `--vtml20` | `test/VTML20.out` | substitution matrix for the search |
| `--outdir` | `results` | where `<sample>/` directories are published |
| `--dart`, `--mmseqs` | `dart`, `mmseqs` | binaries; the `local` profile sets the cluster paths |
| `--decoy_every` | 10 | every N-th read gives a decoy peptide; must be at least 1 |
| `--domain` | `meta` | `predict --domain` |
| `--library_type` | `auto` | `profile --library-type` |
| `--search_params` | tutorial options | MMseqs2 options other than matrix, output format and threads |
| `--null_p_max` | 0.05 | `damage-annotate --null-p-max` |
| `--annotate_params` | `--auto-calibrate-spurious` | extra `damage-annotate` options |
| `--scratch` | `true` | where the search and annotation tasks run: `true` = node `$TMPDIR` or `/tmp`, a path, or `false` = work dir; `local` and `slurm` profiles set `/local/scratch` |
| `--expected` | none | per-read table to compare against; set by `-profile test` |

Tutorial and command reference on the DART wiki:
[Quick Start Tutorial](https://github.com/genomewalker/dart/wiki/Quick-Start-Tutorial),
[Command Reference](https://github.com/genomewalker/dart/wiki/Command-Reference).
