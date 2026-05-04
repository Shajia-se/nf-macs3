# nf-macs3

`nf-macs3` is a Nextflow DSL2 module for ChIP-seq peak calling with MACS3.

## What This Module Does

For each ChIP sample:
1. Resolve treatment BAMs and optional control BAMs (explicit sheet or auto from `samples_master`).
2. Run MACS3 three times by default:
   - `idr_q<value>` for the IDR branch
   - `consensus_q<value>` for the relaxed consensus branch
   - `strict_q<value>` for the consensus / DiffBind branch
3. Filter called peaks against blacklist BED (`bedtools intersect -v`) for both profiles.
4. Write outputs into profile-specific subdirectories.
5. Skip only profile/sample outputs that already exist.

## Input Modes

1. `--macs3_samplesheet`
2. `--samples_master`

This module now supports both:

- standard ChIP-seq with control/input BAM
- exploratory / test runs without control/input BAM

### Explicit samplesheet

CSV header:
```text
sample_id,treatment_bam,control_bam
```

Notes:
- `sample_id` is required
- `treatment_bam` is optional
- `control_bam` is optional
- if `treatment_bam` is empty, the module resolves `${chipfilter_output}/${sample_id}*.nomulti.bam`
- if `control_bam` is empty, MACS3 runs without `-c`

This is useful when treatment BAMs are produced inside the same run, while control BAMs come from a shared external location.

### Auto from `samples_master`

Required columns:
```text
sample_id,is_control
```

Optional columns used:
```text
library_type,enabled,control_id
```

Auto-pairing rules:
- enabled non-control `chip` rows become treatment samples
- control BAM comes from `control_id` when available
- if `control_id` is empty and exactly one enabled control exists, it is used as default
- BAM paths are resolved from `${chipfilter_output}/${sample_id}*.nomulti.bam`
- missing or ambiguous treatment BAM matches fail early
- if no control is resolved, MACS3 runs without `-c`

## Output Layout

Under `${project_folder}/${macs3_output}`:

- `idr_q<value>/${sample}_peaks.narrowPeak`
- `idr_q<value>/${sample}_peaks.xls`
- `idr_q<value>/${sample}_peaks.blacklist_applied.txt`
- `consensus_q<value>/${sample}_peaks.narrowPeak`
- `consensus_q<value>/${sample}_peaks.xls`
- `consensus_q<value>/${sample}_peaks.blacklist_applied.txt`
- `strict_q<value>/${sample}_peaks.narrowPeak`
- `strict_q<value>/${sample}_peaks.xls`
- `strict_q<value>/${sample}_peaks.blacklist_applied.txt`
- plus optional side outputs in each profile directory:
  - `${sample}_summits.bed`
  - `${sample}_treat_pileup.bdg`
  - `${sample}_control_lambda.bdg`

## Key Parameters

- `macs3_samplesheet`: explicit treatment/control sheet
- `samples_master`: auto-pairing source
- `chipfilter_output`: location of `*.nomulti.bam`
- `seq`: `"paired"` (default) or `"single"`
- `genome_size`: effective genome size argument for MACS3 (`mm` default)
- `idr_qvalue`: q-value for IDR branch (default: `0.1`)
- `consensus_qvalue`: q-value for relaxed consensus branch (default: `0.05`)
- `strict_qvalue`: q-value for strict branch (default: `0.01`)
- `run_idr_branch`: enable IDR-oriented MACS3 branch (default: `true`)
- `run_consensus_branch`: enable relaxed consensus MACS3 branch (default: `true`)
- `run_strict_branch`: enable strict MACS3 branch (default: `true`)
- `keep_dup`: duplicate policy (default: `all`)
- `call_summits`: whether to output summits (default: `true`)
- `write_bedgraph`: whether to ask MACS3 to write `*_treat_pileup.bdg` and `*_control_lambda.bdg` (default: `false`)
- `peak_type`: `''` (narrow) or `'--broad'`
- `peak_blacklist_bed`: blacklist BED applied on called peaks (set to empty/null to disable)
- `peak_blacklist_fraction`: overlap fraction used by `bedtools intersect -f` (default: `0.5`)

## Run

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output
```

Explicit sheet:
```bash
nextflow run main.nf -profile hpc \
  --macs3_samplesheet /path/to/macs3_samplesheet.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output
```

No-input test mode via `samples_master`:
- omit the input row, or leave chip `control_id` empty
- MACS3 will run with treatment BAM only

## Downstream Branching

- `nf-idr` should use `idr_q0.1`
- `nf-diffbind` should use `strict_q0.01`
- consensus peak modules should also use `strict_q0.01`
- `strict_q0.01` is the default strict consensus branch for DiffBind and the existing consensus workflow
- `consensus_q0.05` is the relaxed consensus branch for a parallel downstream analysis path

## Interpretation Note

- Runs with control/input are preferred for standard peak calling.
- Runs without control are supported for testing or exploratory analysis, but they are expected to be less stringent.
- `*.bdg` files are optional large intermediate-style outputs and are disabled by default in this pipeline.
