#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def macs3_output = params.macs3_output ?: "macs3_output"

process macs3 {
  tag "${bam.simpleName}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${macs3_output}", mode: 'copy'

  input:
    path bam

  script:
  """
  set -eux

  sample_name=\$(basename ${bam} .bam)

  if [[ "${params.seq}" == "paired" ]]; then
    FORMAT="BAMPE"
  else
    FORMAT="BAM"
  fi

  mkdir -p ${params.project_folder}/${macs3_output}

  macs3 callpeak \\
    -t ${bam} \\
    -f \$FORMAT \\
    -q ${params.qvalue} \\
    --keep-dup all \\
    -g ${params.genome_size} \\
    -n \$sample_name \\
    -B \\
    --outdir ${params.project_folder}/${macs3_output} \\
    ${params.peak_type ?: ''}
  """
}

workflow {

  Channel
    .fromPath("${params.chipfilter_output}/*.clean.bam")
    .filter { bam ->
      ! file("${params.project_folder}/${macs3_output}/${bam.simpleName}_peaks.narrowPeak").exists()
    }
    .set { bam_ch }

  macs3(bam_ch)
}
