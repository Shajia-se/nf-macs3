# nf-macs3

`nf-macs3` is a Nextflow DSL2 module for MACS3 peak calling.

## Input Modes

### 1) Recommended: samplesheet (treatment + control)
CSV header:
```text
sample_id,treatment_bam,control_bam
```

### 2) Fallback: treatment-only
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
- `allow_no_control`: allow no-control mode (default: `false`)
- `qvalue`: default `0.01` (aligned to colleague script)
- `call_summits`: add `--call-summits` (default: `true`)
- `keep_dup`: MACS3 `--keep-dup` value (default: `all`)
- `peak_type`: `''` (narrow) or `--broad`

## Run

```bash
nextflow run main.nf -profile hpc --macs3_samplesheet /path/to/macs3_samplesheet.csv
```
