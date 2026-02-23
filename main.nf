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
    path "${sample_id}_peaks.xls"
    path "${sample_id}_summits.bed", optional: true
    path "${sample_id}_treat_pileup.bdg", optional: true
    path "${sample_id}_control_lambda.bdg", optional: true

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
    path "${sample_id}_peaks.xls"
    path "${sample_id}_summits.bed", optional: true
    path "${sample_id}_treat_pileup.bdg", optional: true
    path "${sample_id}_control_lambda.bdg", optional: true

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
        !(file("${outdir}/${sid}_peaks.${peak_ext}").exists() && file("${outdir}/${sid}_peaks.xls").exists())
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
          !(file("${outdir}/${sid}_peaks.${peak_ext}").exists() && file("${outdir}/${sid}_peaks.xls").exists())
        }
      macs3_no_control(no_control)
    }
  } else if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    def header = null
    def records = []
    master.eachLine { line, n ->
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (n == 1) {
        header = cols
      } else {
        def rec = [:]
        header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
        records << rec
      }
    }

    assert header : "samples_master header not found: ${params.samples_master}"
    assert header.contains('sample_id') : "samples_master missing required column: sample_id"
    assert header.contains('is_control') : "samples_master missing required column: is_control"
    assert header.contains('control_id') : "samples_master missing required column: control_id"

    def isEnabled = { rec ->
      def v = rec.enabled?.toString()?.trim()?.toLowerCase()
      (v == null || v == '' || v == 'true')
    }
    def isControl = { rec ->
      rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
    }
    def isChip = { rec ->
      def lt = rec.library_type?.toString()?.trim()?.toLowerCase()
      !isControl(rec) && (lt == null || lt == '' || lt == 'chip')
    }

    def enabledRecords = records.findAll { rec -> isEnabled(rec) }
    def controlRecords = enabledRecords.findAll { rec -> isControl(rec) }
    def controlIds = controlRecords.collect { it.sample_id?.toString()?.trim() }.findAll { it }
    def controlSet = controlIds as Set
    def defaultControl = (controlIds.size() == 1) ? controlIds[0] : null

    def bamDir = file(params.chipfilter_output)
    assert bamDir.exists() : "chipfilter_output directory not found: ${params.chipfilter_output}"

    def resolveCleanBam = { sid ->
      def hits = bamDir.listFiles()?.findAll { f ->
        f.isFile() && f.name.endsWith('.clean.bam') && f.name.startsWith("${sid}")
      } ?: []

      if (hits.isEmpty()) {
        throw new IllegalArgumentException("No clean BAM found for sample_id '${sid}' under: ${params.chipfilter_output}")
      }
      if (hits.size() > 1) {
        def names = hits.collect { it.name }.join(', ')
        throw new IllegalArgumentException("Multiple clean BAM files matched sample_id '${sid}': ${names}")
      }
      file(hits[0].absolutePath)
    }

    def autoWithControl = []
    def autoNoControl = []

    enabledRecords.findAll { rec -> isChip(rec) }.each { rec ->
      def sid = rec.sample_id?.toString()?.trim()
      if (!sid) return

      def treatBam = resolveCleanBam(sid)
      def cid = rec.control_id?.toString()?.trim()
      if (!cid && defaultControl) cid = defaultControl

      if (cid) {
        assert controlSet.contains(cid) : "control_id '${cid}' for sample '${sid}' is not an enabled control sample in samples_master"
        def controlBam = resolveCleanBam(cid)
        autoWithControl << tuple(sid, treatBam, controlBam)
      } else {
        if (!allow_no_control) {
          throw new IllegalArgumentException("No control_id found for sample '${sid}'. Add control_id in samples_master or set --allow_no_control true.")
        }
        autoNoControl << tuple(sid, treatBam)
      }
    }

    def with_control_ch = Channel.fromList(autoWithControl)
    if (!allow_no_control) {
      with_control_ch = with_control_ch.ifEmpty { exit 1, "ERROR: No treatment/control pairs generated from samples_master: ${params.samples_master}" }
    }
    with_control_ch = with_control_ch.filter { sid, tbam, cbam ->
        !(file("${outdir}/${sid}_peaks.${peak_ext}").exists() && file("${outdir}/${sid}_peaks.xls").exists())
      }

    macs3_with_control(with_control_ch)

    if (allow_no_control) {
      def no_control_ch = Channel
        .fromList(autoNoControl)
        .filter { sid, tbam ->
          !(file("${outdir}/${sid}_peaks.${peak_ext}").exists() && file("${outdir}/${sid}_peaks.xls").exists())
        }
      macs3_no_control(no_control_ch)
    }
  } else {
    def treat_only = Channel
      .fromPath("${params.chipfilter_output}/*.clean.bam", checkIfExists: true)
      .ifEmpty { exit 1, "ERROR: No treatment BAM found under ${params.chipfilter_output}. Provide --macs3_samplesheet for treatment/control mode." }
      .map { bam -> tuple(bam.simpleName, bam) }
      .filter { sid, bam ->
        !(file("${outdir}/${sid}_peaks.${peak_ext}").exists() && file("${outdir}/${sid}_peaks.xls").exists())
      }

    if (!allow_no_control) {
      exit 1, "ERROR: Control BAM is required. Please provide --macs3_samplesheet with columns: sample_id,treatment_bam,control_bam ; or set --allow_no_control true."
    }
    macs3_no_control(treat_only)
  }
}
