# nf-macs3

`nf-macs3` is a Nextflow DSL2 module for ChIP-seq peak calling with MACS3.

## What This Module Does

For each treatment sample:
1. Resolve treatment/control BAM pairs (explicit sheet or auto from `samples_master`).
2. Run `macs3 callpeak` in paired-end (`BAMPE`) or single-end (`BAM`) mode.
3. Produce peak set and standard MACS3 side outputs.
4. Skip samples already completed in output directory.

## Input Modes (Priority Order)

1. `--macs3_samplesheet` (highest priority)
2. `--samples_master` auto-pairing
3. treatment-only fallback (requires `--allow_no_control true`)

### Mode 1: Explicit samplesheet

CSV header:
```text
sample_id,treatment_bam,control_bam
```

- `sample_id`: output prefix
- `treatment_bam`: treatment BAM path
- `control_bam`: control/input BAM path

### Mode 2: Auto from `samples_master`

Required columns:
```text
sample_id,is_control,control_id
```

Optional columns used:
```text
library_type,enabled
```

Auto-pairing rules:
- enabled non-control (`chip`) rows become treatment samples
- control BAM comes from `control_id`
- if `control_id` is empty and exactly one enabled control exists, it is used as default
- BAM paths are resolved from `${chipfilter_output}/${sample_id}*.clean.bam`
- ambiguous/missing BAM matches fail early with explicit error

### Mode 3: No-control fallback

- uses all `${chipfilter_output}/*.clean.bam` as treatment
- requires:
```bash
--allow_no_control true
```

## Output

Under `${project_folder}/${macs3_output}`:
- `${sample}_peaks.narrowPeak` (default)
- `${sample}_peaks.broadPeak` (when `--peak_type "--broad"`)
- `${sample}_peaks.xls`
- `${sample}_summits.bed` (if summits enabled)
- `${sample}_treat_pileup.bdg`
- `${sample}_control_lambda.bdg` (control mode)

## Key Parameters

- `macs3_samplesheet`: explicit treatment/control sheet (priority input)
- `samples_master`: auto-pairing source
- `chipfilter_output`: location of `*.clean.bam`
- `allow_no_control`: allow treatment-only mode (default: `false`)
- `seq`: `"paired"` (default) or `"single"`
- `genome_size`: effective genome size argument for MACS3 (`mm` default)
- `qvalue`: significance threshold (default: `0.01`)
- `keep_dup`: MACS3 duplicate policy (default: `all`)
- `call_summits`: whether to output summits (default: `true`)
- `peak_type`: `''` (narrow) or `'--broad'`

## Run

Explicit samplesheet:
```bash
nextflow run main.nf -profile hpc \
  --macs3_samplesheet /path/to/macs3_samplesheet.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output
```

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output
```

No-control mode:
```bash
nextflow run main.nf -profile hpc \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output \
  --allow_no_control true
```

Resume:
```bash
nextflow run main.nf -profile hpc -resume
```

## Notes

- Control mode is recommended for TF ChIP-seq to reduce false positives.
- Keep MACS3 parameters consistent across compared groups.
- Peak count changes are expected when using stricter upstream filtering (dedup/MAPQ/blacklist).
