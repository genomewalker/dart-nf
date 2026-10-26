# dart-nf tutorial

This tutorial runs the bundled test, explains each step and the files it writes, shows how to
read the results, and then runs real samples on a Slurm cluster.

## 1. Install

You need:

- **Nextflow** 25.10 or later (with Java 17 or later)
- **DART** v2.0.1 or later
- **MMseqs2** 15.6f452
- **Python** 3.8 or later, for the report step

Build DART:

```bash
git clone --branch v2.0.1 https://github.com/genomewalker/dart.git
cd dart
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
./build/dart --version
```

or download `dart-linux-x86_64` from the [DART release page](https://github.com/genomewalker/dart/releases/latest)
and `chmod +x` it.

You don't need to clone dart-nf: Nextflow fetches it from GitHub when you run
`nextflow run genomewalker/dart-nf`. Use `-r v2.0.1` to pin the version.

## 2. Run the test

```bash
nextflow run genomewalker/dart-nf -r v2.0.1 -profile test \
  --dart /path/to/dart --mmseqs /path/to/mmseqs --outdir results
```

If MMseqs2 is not installed, add `conda` to the profile (`-profile test,conda`) and Nextflow
installs it from Bioconda for the search step.

The test data are 50,000 synthetic ancient reads (about 25% terminal C→T damage) and the 13,681
reference proteins they come from. The run takes about a minute on 8 cores and should end with:

```
[SUCCESS] completed=7 failed=0 cached=0
```

The last step, COMPARE_EXPECTED, compares your per-read table with the one shipped in
`test/expected/` and fails the run if any row differs. `results/tutorial/tutorial.compare.txt`
should end with `RESULT: IDENTICAL`. The test fixes every step at 8 CPUs because that is how the
expected table was made.

All files land in `results/tutorial/`.

The test only checks that the pipeline runs and gives the expected output. Two things stay
idle in it. The chance-hit filter flags nothing, because the test database is too small to
produce chance matches. No protein is called ancient, because the profile of the synthetic
library does not yield a usable damage rate. Both are explained below, with what to expect on
real data.

## 3. What each step does

The pipeline runs these steps for every sample in the samplesheet:

```
reads ─► PROFILE ─► PREDICT ─► MMSEQS_SEARCH ─► HITS2EMI ─► DAMAGE_ANNOTATE ─► REPORT
```

### PROFILE: measure the library's damage

`dart profile` looks at the ends of the reads and measures how much deamination damage the
library carries: how often a C reads as T near the 5′ end (and G as A, or C as T, near the 3′
end), and how fast that decays along the read. It also decides whether the library looks
double- or single-stranded.

| File | What it is |
|---|---|
| `<s>.profile.json` | The full damage profile |
| `<s>.library.tdmg` | A compact binary summary that `predict` reads |
| `<s>.profile.log` | Log |

Damage must be measured on the reads **before deduplication**; the estimate is unreliable on
deduplicated reads. That is why the samplesheet has a separate `profile_reads` column.

### PREDICT: translate reads into proteins

`dart predict` finds the protein-coding frame in each read. Damage can turn a normal codon
into a stop codon (CAA → TAA, for example), so DART weighs each stop by how likely it is to be
damage before rejecting a frame. A read can give more than one candidate protein (ORF).

It also writes **decoys**: the reversed protein of every 10th read. Reversed proteins have the
length and composition of real ones but no true matches, so whatever they hit is a chance
match.

| File | What it is |
|---|---|
| `<s>.predictions.gff` | One row per predicted ORF |
| `<s>.proteins.faa` | The predicted proteins; ids are `<read>_<strand>_<frame>_<ordinal>` |
| `<s>.proteins.decoy.faa` | Decoy peptides, ids prefixed `DECOY_` |
| `<s>.predictions.agd` | Damage index: per-read damage evidence, read by the later steps |
| `<s>.predict.log` | Log |

In the test: 147,871 ORFs from 50,000 reads, and 14,814 decoys.

### MMSEQS_SEARCH: search the reference database

MMseqs2 searches the proteins against `--db`, then searches the decoys with exactly the same
options. The default options suit short ancient peptides (VTML20 matrix, minimum 12 aa
alignments, high identity).

| File | What it is |
|---|---|
| `<s>.hits.tsv` | Hits of the real proteins (16 columns) |
| `<s>.decoy_hits.tsv` | Hits of the decoys |
| `<s>.mmseqs.log` | Log |

In the test: 45,714 alignments for 21,904 ORFs. The decoys get only 3 alignments, because the
test database holds only the proteins the reads came from; there is almost nothing to match by
chance.

### HITS2EMI: index the hits and fit the chance-hit model

`dart hits2emi` turns the hits into a compact columnar index (EMI) and fits the chance-hit
model from the decoy hits: for each range of alignment lengths, how high a score a decoy
reaches by chance.

| File | What it is |
|---|---|
| `<s>.hits.emi` | Hit index |
| `<s>.hits.emi.null` | Chance-hit model (the "null") |
| `<s>.hits2emi.log` | Log |

In the test the log warns that no length bin could be fitted: one decoy with hits is not enough.
That is expected for this tiny database. The run continues without chance-hit filtering. With a
real database the decoys get thousands of hits and the bins are fitted.

### DAMAGE_ANNOTATE: assign reads and score damage

`dart damage-annotate` does three things:

1. **Assigns reads.** A short read often matches several related proteins. An EM algorithm
   shares each read among its candidates according to how well it matches and how much
   support each candidate has from other reads.
2. **Flags chance hits.** Using the null, each read's best hit gets a q-value (`null_q`): the
   expected share of chance matches among reads scoring at least as high in the same length
   range. Reads above 5% (`--null_fdr`) are flagged and left out of abundance.
3. **Scores damage.** Each read gets a `posterior`, the probability that it is ancient, from
   terminal damage, damage-type substitutions in the alignment, and identity. Each protein then
   gets a test of whether it has more damaged reads than modern DNA would.

| File | What it is |
|---|---|
| `<s>.annotated.tsv` | One row per read whose best alignment has a mismatch |
| `<s>.gene_summary.tsv` | One row per gene with enough reads |
| `<s>.protein_summary.tsv` | One row per reference protein with assigned reads |
| `<s>.damage_annotate.log` | Log, including the chance-hit and EM summaries |

### REPORT

`<s>.report.md` gathers the counts, the null model, the EM summary, the damage distributions and
the top proteins into one page. `<s>.versions.txt` records the dart, libtaph, MMseqs2,
Nextflow and dart-nf versions used.

## 4. Reading the results

Start with `report.md`, then go to the tables.

### report.md

- **Counts.** Reads in, ORFs, decoys, hits, and rows in each table. A low share of reads with a
  hit usually means the database lacks the organisms in the sample, or the reads are very short.
- **Null model.** The fitted bins and, per length range, how many reads were flagged as chance
  hits and the lowest score kept. Expect the short ranges (7–12 aa) to lose the most reads;
  long alignments are rarely chance hits. If every bin says `unfit`, no filtering happened.
- **EM assignment.** `Unique mappers` hit one reference only; `Multi-mappers` were shared by the
  EM. `Reassigned vs bitscore-best` is the share of reads whose final reference is not their
  top-scoring hit.
- **Per-read damage.** The posterior histogram and the read classes. An ancient sample has many
  reads above 0.5; a modern one piles up near 0.
- **Proteins.** The authentication tiers and calls, and the top proteins.
- **Genes.** How many genes pass the read thresholds, and how many look mostly ancient.

In the test you should see 413 genes, 262 of them with `ancient_frac >= 0.5`, and every protein
`unauthenticatable`. The synthetic reads give a library profile whose damage state is
undetermined (`library_damage_state` 3), and without a usable library damage rate no protein can
be called ancient. The report says so and ranks proteins by their ungated damage score instead.
On a real ancient library with clear damage the state should be `DETECTED` (0).

### annotated.tsv (per read)

The columns you will use most:

| Column | Meaning |
|---|---|
| `query_id`, `target_id` | The ORF and the reference the EM assigned it to |
| `posterior` | Probability that the read is ancient |
| `is_damaged` | 1 when `posterior ≥ 0.7` |
| `classification` | `ancient_confident`, `ancient_likely`, `undetermined` or `modern_confident` |
| `tier` | Evidence present: `terminal` (read ends), `sites` (substitutions), `convergent` (both) or `none` |
| `gamma` | EM weight of this assignment |
| `null_q` | Chance-hit q-value; above 0.05 the read is flagged |
| `em_reject_reason` | 0 = clean; 1 = chance hit; 2 = low EM margin; 4 = many competing references (values add) |

Reads that align perfectly carry no substitution evidence and have no row, but the EM still
assigned them and the summaries count them.

### gene_summary.tsv and protein_summary.tsv

- **Abundance:** use `eff_reads`, the EM-weighted read count. `n_reads` counts every candidate
  assignment, so a read shared by two proteins counts on both.
- **Ancient share of a gene:** `ancient_frac` with its interval `ci_low`–`ci_high`.
- **Is the protein ancient?** Use `auth_binom_call` together with `protein_auth_tier`. Trust an
  `ancient` call when the tier is `detected`. `advisory` means the library's damage was not
  firmly detected. `undetermined_low_coverage` means too few reads to decide, not "modern".
- **Coverage:** `breadth`, `depth_evenness` and `positional_score` show whether the reads spread
  along the protein or pile up in one spot (a sign of a shared domain or a spurious match).

Full column definitions are on the DART wiki:
[Output Formats](https://github.com/genomewalker/dart/wiki/Output-Formats) and
[Interpreting Results](https://github.com/genomewalker/dart/wiki/Interpreting-Results).

## 5. Running your own samples on Slurm

### Samplesheet

```csv
sample,reads,profile_reads
S1,/data/S1.dedup.fq.gz,/data/S1.merged.fq.gz
S2,/data/S2.dedup.fq.gz,/data/S2.merged.fq.gz
```

- `reads`: the reads to translate, usually deduplicated.
- `profile_reads`: the merged reads before deduplication, for the damage profile. Leave it
  empty to profile `reads`.
- Relative paths are resolved against the samplesheet's folder.

### Database

`--db` takes a protein FASTA (plain or gzipped) or an MMseqs2 database. For anything larger than
a few hundred thousand proteins, build an indexed database once and reuse it:

```bash
mmseqs createdb uniref90.faa.gz uniref90
mmseqs createindex uniref90 tmp
```

then pass `--db /db/uniref90/uniref90`. Every file starting with that prefix is staged.

### Resources

The `slurm` profile submits each step as a job: 8 CPUs, 16 GB and 4 h by default, 16 CPUs and
32 GB for the search. A large database needs more for the search. Put your cluster's settings in
a small config file:

```groovy
// my.config
process {
    queue = 'compute'
    withName: 'MMSEQS_SEARCH' {
        cpus = 32
        memory = 200.GB
        time = 24.h
    }
}
```

### Run

```bash
nextflow run genomewalker/dart-nf -r v2.0.1 -profile slurm -c my.config \
  --samplesheet samples.csv --db /db/uniref90/uniref90 \
  --dart /path/to/dart --mmseqs /path/to/mmseqs \
  --scratch /local/scratch --outdir results -resume
```

- `--scratch` is where the search and annotation steps run. Point it at the compute nodes'
  local disk. The default (`true`) uses `$TMPDIR`, or `/tmp` if that is unset.
- `-resume` reuses finished steps if you rerun after a failure or change.
- Run Nextflow itself inside `screen`/`tmux` or as a small long-running job, since it has to stay
  alive while the jobs run.

## 6. One null per study

The chance-hit model depends on the database and the search options, not on the sample. For a
study with many samples searched the same way, you can fit one null and share it:

1. Run the pipeline on a few representative samples.
2. Pool their decoy hits:

   ```bash
   dart null -o study.null --db-label uniref90 \
     results/S1/S1.decoy_hits.tsv results/S2/S2.decoy_hits.tsv results/S3/S3.decoy_hits.tsv
   dart damage-info study.null    # check that the bins are fitted
   ```

3. Run the remaining samples with the pooled model:

   ```bash
   --annotate_params '--em-null-model /abs/path/study.null --auto-calibrate-spurious'
   ```

   Use an absolute path. The pipeline still searches decoys for every sample; raise
   `--decoy_every` (for example to 100) to make that search cheap.

Refit the null when the database, the search options, the library type or the read lengths
change.

## 7. Keeping short hits: `--keep-null-flagged`

By default, reads flagged as chance hits are left out of gene and protein abundance. This keeps
spurious matches out of your tables, but it also removes some real short hits: on a test sample
of known content, about half of the true 10–12 aa hits were dropped.

If very short peptides matter for your question, keep them:

```bash
--annotate_params '--auto-calibrate-spurious --keep-null-flagged'
```

The flagged reads are always listed in `annotated.tsv` with their `null_q` either way, so you can
also apply your own cutoff afterwards. `damage_annotate.log` reports how many reads the default
dropped (`null drop:`).

## 8. Common failures

| Symptom | Cause and fix |
|---|---|
| `set --samplesheet` or `set --db` | A required parameter is missing. |
| `dart: command not found` or `mmseqs: command not found` | The binaries are not on the job's `PATH`. Pass full paths with `--dart` and `--mmseqs`. |
| `Unknown option` in a step's log | The dart binary is older than v2.0.1. Use v2.0.1 or later. |
| DAMAGE_ANNOTATE fails while opening the EM tape | The scratch directory is on tmpfs or a file system without `O_DIRECT`. Set `--scratch` to a local disk. |
| MMSEQS_SEARCH killed (exit 137 or 140) | Out of memory or time. Raise `memory` and `time` for `MMSEQS_SEARCH` in your config. |
| `no null bin fitted` warning in `hits2emi.log` | Too few decoy hits to fit the chance-hit model, so no filtering happens. Normal for the test. On real data, use a larger database or pool decoys across samples with `dart null` (section 6). |
| `library_not_authenticatable` warning | The library has no usable damage rate, so no protein can be called ancient. Check that `profile_reads` points to the reads before deduplication, and look at `profile.json`. |
| COMPARE_EXPECTED reports differences | Your dart version differs from the one that made `test/expected/`, or the test ran with a different CPU count. |
| A sample is missing from the results | Check its row in the samplesheet: the `sample` name and the paths. |

To see why a step failed, look at its log in `results/<sample>/`, or at `.command.err` in the
work directory that Nextflow prints for the failed task.
