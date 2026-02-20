#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def macs3_output = params.macs3_output ?: "macs3_output"
def allow_no_control = (params.allow_no_control == null) ? false : params.allow_no_control

process macs3_with_control {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${macs3_output}", mode: 'copy'

  input:
    tuple val(sample_id), path(treat_bam), path(control_bam)

  output:
    tuple val(sample_id), path("${sample_id}_peaks.*Peak")

  script:
  def keep_dup = params.keep_dup ?: 'all'
  def call_summits = (params.call_summits == null || params.call_summits) ? '--call-summits' : ''
  """
  set -eux
  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp

  if [[ "${params.seq}" == "paired" ]]; then
    FORMAT="BAMPE"
  else
    FORMAT="BAM"
  fi

  macs3 callpeak \\
    -t ${treat_bam} \\
    -c ${control_bam} \\
    -f \$FORMAT \\
    -q ${params.qvalue} \\
    --keep-dup ${keep_dup} \\
    -g ${params.genome_size} \\
    -n ${sample_id} \\
    -B \\
    ${call_summits} \\
    --outdir . \\
    ${params.peak_type ?: ''}
  """
}

process macs3_no_control {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${macs3_output}", mode: 'copy'

  input:
    tuple val(sample_id), path(treat_bam)

  output:
    tuple val(sample_id), path("${sample_id}_peaks.*Peak")

  script:
  def keep_dup = params.keep_dup ?: 'all'
  def call_summits = (params.call_summits == null || params.call_summits) ? '--call-summits' : ''
  """
  set -eux
  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp

  if [[ "${params.seq}" == "paired" ]]; then
    FORMAT="BAMPE"
  else
    FORMAT="BAM"
  fi

  macs3 callpeak \\
    -t ${treat_bam} \\
    -f \$FORMAT \\
    -q ${params.qvalue} \\
    --keep-dup ${keep_dup} \\
    -g ${params.genome_size} \\
    -n ${sample_id} \\
    -B \\
    ${call_summits} \\
    --outdir . \\
    ${params.peak_type ?: ''}
  """
}

workflow {
  def peak_ext = (params.peak_type ?: '').contains('--broad') ? 'broadPeak' : 'narrowPeak'
  def outdir = "${params.project_folder}/${macs3_output}"

  if (params.macs3_samplesheet) {
    def rows = Channel
      .fromPath(params.macs3_samplesheet, checkIfExists: true)
      .splitCsv(header: true)

    def with_control = rows
      .filter { row -> row.control_bam && row.control_bam.toString().trim() }
      .map { row ->
        def sid = row.sample_id?.toString()?.trim() ?: file(row.treatment_bam.toString()).simpleName
        tuple(sid, file(row.treatment_bam.toString()), file(row.control_bam.toString()))
      }
      .filter { sid, tbam, cbam ->
        !file("${outdir}/${sid}_peaks.${peak_ext}").exists()
      }

    macs3_with_control(with_control)

    if (allow_no_control) {
      def no_control = rows
        .filter { row -> !(row.control_bam && row.control_bam.toString().trim()) }
        .map { row ->
          def sid = row.sample_id?.toString()?.trim() ?: file(row.treatment_bam.toString()).simpleName
          tuple(sid, file(row.treatment_bam.toString()))
        }
        .filter { sid, tbam ->
          !file("${outdir}/${sid}_peaks.${peak_ext}").exists()
        }
      macs3_no_control(no_control)
    }
  } else {
    def treat_only = Channel
      .fromPath("${params.chipfilter_output}/*.clean.bam", checkIfExists: true)
      .ifEmpty { exit 1, "ERROR: No treatment BAM found under ${params.chipfilter_output}. Provide --macs3_samplesheet for treatment/control mode." }
      .map { bam -> tuple(bam.simpleName, bam) }
      .filter { sid, bam ->
        !file("${outdir}/${sid}_peaks.${peak_ext}").exists()
      }

    if (!allow_no_control) {
      exit 1, "ERROR: Control BAM is required. Please provide --macs3_samplesheet with columns: sample_id,treatment_bam,control_bam ; or set --allow_no_control true."
    }
    macs3_no_control(treat_only)
  }
}
