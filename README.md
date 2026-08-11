# nf-macs3

`nf-macs3` calls ChIP-seq peaks with MACS3.

This module is normally the sixth step in the ChIP-seq pipeline:

```text
nf-fastqc -> nf-fastp -> nf-bwa -> nf-picard -> nf-chipfilter -> nf-macs3
```

`nf-chipfilter` produces MAPQ-filtered `*.nomulti.bam` files. `nf-macs3` uses those BAM files to call peaks for downstream IDR, consensus peak analysis, annotation, and visualization.

## Can It Connect To The Previous Modules?

Yes, by interface design.

The upstream `nf-chipfilter` output looks like:

```text
sample.sorted.dedup.nomulti.bam
sample.sorted.dedup.nomulti.bam.bai
```

`nf-macs3` resolves BAMs by `sample_id` prefix and `.nomulti.bam` suffix, so `sample_id=sample` correctly matches:

```text
sample*.nomulti.bam
```

This means it is compatible with the current outputs from:

```text
nf-fastqc -> nf-fastp -> nf-bwa -> nf-picard -> nf-chipfilter
```

Actual end-to-end execution still needs to be tested in an environment where `nextflow` is installed.

## What It Does

1. Finds treatment BAMs from the `nf-chipfilter` output folder.
2. Optionally finds matched control/input BAMs.
3. Runs MACS3 peak calling.
4. Runs one or more peak-calling branches:
   - `idr_q0.1`
   - `consensus_q0.05`
   - `strict_q0.01`
5. Optionally filters peaks against a blacklist BED file.
6. Writes outputs into branch-specific folders.

## Before You Run

You need:

- Nextflow installed
- MAPQ-filtered BAM files from `nf-chipfilter`
- A `samples_master` CSV or a MACS3 samplesheet
- For local runs: Docker available
- For HPC runs: Slurm and Singularity available
- For HPC runs: MACS3 and Bedtools Singularity image paths in `configs/slurm.config` must exist

Default HPC notification email:

```text
molendo.hpc@gmail.com
```

## Treatment And Input/Control BAMs

MACS3 peak calling usually compares a ChIP treatment BAM against an input/control BAM:

```bash
macs3 callpeak -t treatment.bam -c input_control.bam
```

In this module:

- `treatment_bam` is the ChIP/IP sample, for example H3K27ac, H3K4me3, TF ChIP, CUT&Tag-like target sample.
- `control_bam` is the matching input/control sample used by MACS3 with `-c`.
- If no `control_bam` is provided or resolved, MACS3 runs treatment-only with `-t` and no `-c`.

There are two common ways to provide the input/control BAM.

### Scenario A: Input/Control Is Part Of The Same Pipeline Run

Use this when the input/control sample has raw FASTQ files and should go through the same upstream workflow as the ChIP samples:

```text
nf-fastqc -> nf-fastp -> nf-bwa -> nf-picard -> nf-chipfilter
```

In this case, the input/control sample is listed in `samples_master` with:

```text
is_control=true
```

The ChIP sample points to it with:

```text
control_id=<input sample_id>
```

Example:

```csv
sample_id,library_type,is_control,control_id,enabled
WT_input,input,true,,true
WT_H3K27ac,chip,false,WT_input,true
```

After upstream processing, both samples must have MAPQ-filtered BAMs in `chipfilter_output`:

```text
WT_input*.nomulti.bam
WT_H3K27ac*.nomulti.bam
```

Then `nf-macs3` automatically uses:

```text
-t WT_H3K27ac*.nomulti.bam
-c WT_input*.nomulti.bam
```

### Scenario B: Input/Control BAM Already Exists

Use this when the input/control BAM is already available and should not be rerun through this pipeline.

This is the shared-input style used in cases like the Marjolein liver run: the ChIP treatment BAM comes from the current pipeline, while the input/control BAM may already exist from a previous project or shared reference analysis.

In this case, use `--macs3_samplesheet` and write the existing control BAM path directly:

```csv
sample_id,treatment_bam,control_bam
WT_H3K27ac,,/path/to/existing_liver_input.nomulti.bam
```

If `treatment_bam` is empty, `nf-macs3` finds it from `chipfilter_output` using `sample_id`:

```text
WT_H3K27ac*.nomulti.bam
```

The command effectively becomes:

```text
-t WT_H3K27ac*.nomulti.bam
-c /path/to/existing_liver_input.nomulti.bam
```

You can also provide both paths explicitly:

```csv
sample_id,treatment_bam,control_bam
WT_H3K27ac,/path/to/WT_H3K27ac.nomulti.bam,/path/to/existing_liver_input.nomulti.bam
```

Important: if you use an existing input/control BAM, make sure it is compatible with the treatment BAM:

- Same genome build
- Same chromosome naming style, for example `chr1` vs `1`
- Similar upstream processing policy when possible
- Sorted and readable BAM
- Preferably indexed with `.bai`, although MACS3 mainly needs the BAM

## Input Option 1: Samples Master Recommended

Use this when running the full ChIP-seq pipeline.

Required columns:

```csv
sample_id,is_control
```

Useful optional columns:

```csv
library_type,control_id,enabled
```

Example:

```csv
sample_id,library_type,is_control,control_id,enabled
WT_input,input,true,,true
WT_H3K27ac,chip,false,WT_input,true
```

Rules:

- Enabled non-control `chip` rows become treatment samples.
- `control_id` should match an enabled input/control `sample_id`.
- If `control_id` is empty and exactly one enabled control exists, that control is used automatically.
- If no control is resolved, MACS3 runs treatment-only without `-c`.
- Rows with `enabled=false` are skipped.
- Empty `enabled` values are treated as enabled.

Run with:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --project_folder /path/to/output_project \
  --macs3_output macs3_output
```

## Input Option 2: Explicit MACS3 Samplesheet

Use this when you want to manually specify treatment/control BAMs, especially when the input/control BAM already exists.

CSV header:

```csv
sample_id,treatment_bam,control_bam
```

Example:

```csv
sample_id,treatment_bam,control_bam
WT_H3K27ac,/path/to/WT_H3K27ac.nomulti.bam,/path/to/WT_input.nomulti.bam
```

Notes:

- `sample_id` is required.
- `treatment_bam` is optional. If empty, the module resolves it from `chipfilter_output` using `sample_id`.
- `control_bam` is optional. If empty, MACS3 runs treatment-only without `-c`.
- `control_bam` can point to a shared or previously generated input/control BAM.

Run with:

```bash
nextflow run main.nf -profile hpc \
  --macs3_samplesheet /path/to/macs3_samplesheet.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --project_folder /path/to/output_project \
  --macs3_output macs3_output
```

## Branches

By default, three branches are enabled:

```text
idr_q0.1
consensus_q0.05
strict_q0.01
```

Branch controls:

| Parameter | Meaning | Default |
| --- | --- | --- |
| `run_idr_branch` | Run IDR-oriented branch | `true` |
| `idr_qvalue` | q-value for IDR branch | `0.1` |
| `run_consensus_branch` | Run relaxed consensus branch | `true` |
| `consensus_qvalue` | q-value for relaxed consensus branch | `0.05` |
| `run_strict_branch` | Run strict branch | `true` |
| `strict_qvalue` | q-value for strict branch | `0.01` |

Typical downstream use:

- `idr_q0.1`: input for `nf-idr`
- `strict_q0.01`: input for strict consensus and DiffBind-style workflows
- `consensus_q0.05`: relaxed consensus branch

## Blacklist Filtering

Blacklist filtering is optional.

To enable it:

```bash
--peak_blacklist_bed /path/to/blacklist.bed
```

If provided, peaks overlapping the blacklist are removed with Bedtools and a report is written:

```text
sample_peaks.blacklist_applied.txt
```

If `peak_blacklist_bed` is empty or not provided, blacklist filtering is skipped.

## Output

Results are written to:

```text
${project_folder}/${macs3_output}/
```

Example:

```text
/path/to/output_project/macs3_output/
  idr_q0.1/
    WT_H3K27ac_peaks.narrowPeak
    WT_H3K27ac_peaks.xls
  consensus_q0.05/
    WT_H3K27ac_peaks.narrowPeak
    WT_H3K27ac_peaks.xls
  strict_q0.01/
    WT_H3K27ac_peaks.narrowPeak
    WT_H3K27ac_peaks.xls
```

Optional outputs:

- `*_summits.bed`
- `*_treat_pileup.bdg`
- `*_control_lambda.bdg`
- `*_peaks.blacklist_applied.txt`

BedGraph files can be large and are disabled by default:

```bash
--write_bedgraph false
```

## Recommended HPC Run

From inside the `nf-macs3` folder:

```bash
cd /path/to/nf-macs3

nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --project_folder /path/to/output_project \
  --macs3_output macs3_output \
  --peak_blacklist_bed /path/to/blacklist.bed
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc -resume \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --project_folder /path/to/output_project \
  --macs3_output macs3_output
```

Override the HPC notification email:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --project_folder /path/to/output_project \
  --mail_user molendo.hpc@gmail.com
```

## Local Test Run

Use local mode only for small test BAM files:

```bash
nextflow run main.nf -profile local \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/test_chipfilter_output \
  --project_folder ./test_output
```

## Key Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `samples_master` | Metadata CSV for automatic treatment/control pairing | `null` |
| `macs3_samplesheet` | Explicit treatment/control BAM sheet | `null` |
| `chipfilter_output` | Input folder containing `*.nomulti.bam` files | `null` |
| `macs3_output` | MACS3 output subfolder | `macs3_output` |
| `seq` | `paired` uses `BAMPE`; anything else uses `BAM` | `paired` |
| `genome_size` | MACS3 genome size | `mm` |
| `keep_dup` | MACS3 duplicate policy | `all` |
| `call_summits` | Add `--call-summits` | `true` |
| `write_bedgraph` | Write MACS3 bedGraph outputs | `false` |
| `peak_type` | Empty for narrow peaks; `--broad` for broad peaks | empty |
| `peak_blacklist_bed` | Optional blacklist BED | `null` |
| `peak_blacklist_fraction` | Bedtools overlap fraction | `0.5` |
| `cpus` | CPUs per MACS3 task | `4` |
| `memory` | Memory per MACS3 task | `16GB` |
| `time` | Runtime limit per MACS3 task | `8h` |
| `mail_user` | HPC notification email | `molendo.hpc@gmail.com` |

## Existing Results Are Skipped

For each sample and branch, the module checks whether expected peak outputs already exist:

```text
sample_peaks.narrowPeak
sample_peaks.xls
```

If blacklist filtering is enabled, this report must also exist:

```text
sample_peaks.blacklist_applied.txt
```

If all expected files exist, that sample/branch is skipped. If any file is missing, MACS3 runs again for that sample/branch.

## How To Check Results

Start with:

```text
${project_folder}/${macs3_output}/strict_q0.01/sample_peaks.narrowPeak
${project_folder}/${macs3_output}/strict_q0.01/sample_peaks.xls
```

Important checks:

- Peak count is not zero.
- Treatment samples have stronger peaks than control/input samples.
- Replicates show comparable peak numbers.
- If blacklist filtering was enabled, inspect `sample_peaks.blacklist_applied.txt`.

## Troubleshooting

If the run fails:

1. Check the main Nextflow log:

```bash
less .nextflow.log
```

2. Check the failed task error file:

```bash
less work/<hash>/.command.err
```

3. Common problems:

- `--chipfilter_output` does not point to the `nf-chipfilter` output folder.
- No matching `sample_id*.nomulti.bam` file exists.
- More than one BAM matches the same `sample_id`.
- `control_id` in `samples_master` does not match an enabled control row.
- `peak_blacklist_bed` is set but the BED file does not exist.
- `configs/slurm.config` points to a missing Singularity image.
- The HPC bind path in `extra_mounts` does not include BAM, blacklist, or output locations.
- Docker is not running for local mode.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
