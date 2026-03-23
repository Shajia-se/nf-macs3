# nf-macs3

`nf-macs3` is a Nextflow DSL2 module for ChIP-seq peak calling with MACS3.

## What This Module Does

For each ChIP sample:
1. Resolve treatment/control BAM pairs (explicit sheet or auto from `samples_master`).
2. Run MACS3 three times by default:
   - `idr_q0.1` for the IDR branch
   - `qc_q0.05` for QC / reporting peak-count summaries
   - `strict_q0.01` for the consensus / DiffBind branch
3. Filter called peaks against blacklist BED (`bedtools intersect -v`) for both profiles.
4. Write outputs into profile-specific subdirectories.
5. Skip only profile/sample outputs that already exist.

## Input Modes

1. `--macs3_samplesheet`
2. `--samples_master`

This module now assumes all ChIP samples have a control/input BAM.

### Explicit samplesheet

CSV header:
```text
sample_id,treatment_bam,control_bam
```

### Auto from `samples_master`

Required columns:
```text
sample_id,is_control,control_id
```

Optional columns used:
```text
library_type,enabled
```

Auto-pairing rules:
- enabled non-control `chip` rows become treatment samples
- control BAM comes from `control_id`
- if `control_id` is empty and exactly one enabled control exists, it is used as default
- BAM paths are resolved from `${chipfilter_output}/${sample_id}*.clean.bam`
- missing or ambiguous BAM/control matches fail early

## Output Layout

Under `${project_folder}/${macs3_output}`:

- `idr_q0.1/${sample}_peaks.narrowPeak`
- `idr_q0.1/${sample}_peaks.xls`
- `idr_q0.1/${sample}_peaks.blacklist_applied.txt`
- `qc_q0.05/${sample}_peaks.narrowPeak`
- `qc_q0.05/${sample}_peaks.xls`
- `qc_q0.05/${sample}_peaks.blacklist_applied.txt`
- `strict_q0.01/${sample}_peaks.narrowPeak`
- `strict_q0.01/${sample}_peaks.xls`
- `strict_q0.01/${sample}_peaks.blacklist_applied.txt`
- plus optional side outputs in each profile directory:
  - `${sample}_summits.bed`
  - `${sample}_treat_pileup.bdg`
  - `${sample}_control_lambda.bdg`

## Key Parameters

- `macs3_samplesheet`: explicit treatment/control sheet
- `samples_master`: auto-pairing source
- `chipfilter_output`: location of `*.clean.bam`
- `seq`: `"paired"` (default) or `"single"`
- `genome_size`: effective genome size argument for MACS3 (`mm` default)
- `idr_qvalue`: q-value for IDR branch (default: `0.1`)
- `qc_qvalue`: q-value for QC/reporting branch (default: `0.05`)
- `strict_qvalue`: q-value for strict branch (default: `0.01`)
- `keep_dup`: duplicate policy (default: `all`)
- `call_summits`: whether to output summits (default: `true`)
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

## Downstream Branching

- `nf-idr` should use `idr_q0.1`
- `nf-diffbind` should use `strict_q0.01`
- consensus peak modules should also use `strict_q0.01`
- `qc_q0.05` is an extra reporting branch; it is not wired into downstream modules by default
