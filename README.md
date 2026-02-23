# nf-macs3

`nf-macs3` is a Nextflow DSL2 module for MACS3 peak calling.

## Input Modes

Priority:
1. `--macs3_samplesheet` (explicit treatment/control pairs)
2. `--samples_master` (auto-generate treatment/control pairs)
3. treatment-only fallback (requires `--allow_no_control true`)

### 1) Explicit samplesheet (treatment + control)
CSV header:
```text
sample_id,treatment_bam,control_bam
```

### 2) Auto from `samples_master`
Required columns:
```text
sample_id,is_control,control_id
```
Optional columns used:
```text
library_type,enabled
```
Behavior:
- enabled chip samples are treated as treatment rows
- control BAM is resolved from `control_id`
- if exactly one enabled control exists and `control_id` is empty, that control is used as default
- BAM paths are resolved from `${chipfilter_output}/${sample_id}*.clean.bam`

### 3) Fallback: treatment-only
Reads `${chipfilter_output}/*.clean.bam`, requires:
```bash
--allow_no_control true
```

## Output

Under `${project_folder}/${macs3_output}`:
- `${sample}_peaks.narrowPeak` (default)
- `${sample}_peaks.broadPeak` (if `--peak_type "--broad"`)
- standard MACS3 side outputs (`*_summits.bed`, `*_peaks.xls`, bedGraph files)

## Key Parameters

- `macs3_samplesheet`: treatment/control pairing table
- `samples_master`: optional auto-pairing source
- `allow_no_control`: allow no-control mode (default: `false`)
- `qvalue`: default `0.01` (aligned to colleague script)
- `call_summits`: add `--call-summits` (default: `true`)
- `keep_dup`: MACS3 `--keep-dup` value (default: `all`)
- `peak_type`: `''` (narrow) or `--broad`

## Run

```bash
nextflow run main.nf -profile hpc --macs3_samplesheet /path/to/macs3_samplesheet.csv
```

Auto from `samples_master`:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output
```
