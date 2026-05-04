#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def macs3_output = params.macs3_output ?: "macs3_output"

def resolveBaseDir = { p ->
  def fp = file(p.toString())
  fp.isAbsolute() ? fp : file("${params.project_folder}/${p}")
}

process macs3_callpeak {
  tag "${profile_name}:${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir { "${params.project_folder}/${macs3_output}/${profile_name}" }, mode: 'copy'

  input:
    tuple val(profile_name), val(qval), val(sample_id), path(treat_bam), val(control_bam)

  output:
    tuple val(profile_name), val(sample_id), path("${sample_id}_peaks.*Peak"), emit: peaks
    path "${sample_id}_peaks.xls"
    path "${sample_id}_summits.bed", optional: true
    path "${sample_id}_treat_pileup.bdg", optional: true
    path "${sample_id}_control_lambda.bdg", optional: true

  script:
  def keep_dup = params.keep_dup ?: 'all'
  def call_summits = (params.call_summits == null || params.call_summits) ? '--call-summits' : ''
  def controlArg = control_bam ? "-c ${control_bam}" : ''
  def bedgraphArg = (params.write_bedgraph == null || params.write_bedgraph) ? '-B' : ''
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
    ${controlArg} \\
    -f \$FORMAT \\
    -q ${qval} \\
    --keep-dup ${keep_dup} \\
    -g ${params.genome_size} \\
    -n ${sample_id} \\
    ${bedgraphArg} \\
    ${call_summits} \\
    --outdir . \\
    ${params.peak_type ?: ''}
  """
}

process macs3_blacklist_peak {
  tag "${profile_name}:${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir { "${params.project_folder}/${macs3_output}/${profile_name}" }, mode: 'copy', overwrite: true

  input:
    tuple val(profile_name), val(sample_id), path(peak_file), path(blacklist_bed)

  output:
    tuple val(profile_name), val(sample_id), path("${sample_id}_peaks.*Peak", includeInputs: true)
    path "${sample_id}_peaks.blacklist_applied.txt"

  script:
  def frac = (params.peak_blacklist_fraction ?: 0.5).toString()
  """
  set -euo pipefail
  in_peak="${peak_file}"
  blacklist_bed="${blacklist_bed}"
  ext="\${in_peak##*.}"
  out_peak="${sample_id}_peaks.\${ext}"
  tmp_out="${sample_id}_peaks.filtered.\${ext}"

  bedtools intersect -a "\$in_peak" -b "\$blacklist_bed" -v -f ${frac} | awk '\$1 ~ /^chr/' > "\$tmp_out"
  mv "\$tmp_out" "\$out_peak"

  before_n=\$(wc -l < "\$in_peak" || echo 0)
  after_n=\$(wc -l < "\$out_peak" || echo 0)

  cat > ${sample_id}_peaks.blacklist_applied.txt << EOF
profile\t${profile_name}
sample_id\t${sample_id}
input_peak\t\$in_peak
output_peak\t\$out_peak
blacklist_bed\t\$blacklist_bed
blacklist_fraction\t${frac}
before_peaks\t\$before_n
after_peaks\t\$after_n
EOF
  """
}

workflow {
  def peak_ext = (params.peak_type ?: '').contains('--broad') ? 'broadPeak' : 'narrowPeak'
  def outBase = resolveBaseDir(macs3_output)

  def runIdrBranch = (params.run_idr_branch == null) ? true : params.run_idr_branch.toString().toLowerCase() == 'true'
  def runConsensusBranch = (params.run_consensus_branch == null) ? true : params.run_consensus_branch.toString().toLowerCase() == 'true'
  def runStrictBranch = (params.run_strict_branch == null) ? true : params.run_strict_branch.toString().toLowerCase() == 'true'

  def profiles = []
  if (runIdrBranch) {
    profiles << [profile_name: "idr_q${(params.idr_qvalue ?: 0.1).toString()}", qval: (params.idr_qvalue ?: 0.1).toString()]
  }
  if (runConsensusBranch) {
    profiles << [profile_name: "consensus_q${(params.consensus_qvalue ?: 0.05).toString()}", qval: (params.consensus_qvalue ?: 0.05).toString()]
  }
  if (runStrictBranch) {
    profiles << [profile_name: "strict_q${(params.strict_qvalue ?: 0.01).toString()}", qval: (params.strict_qvalue ?: 0.01).toString()]
  }
  assert !profiles.isEmpty() : "At least one MACS3 branch must be enabled (run_idr_branch/run_consensus_branch/run_strict_branch)"

  if (params.peak_blacklist_bed) {
    def bl = file(params.peak_blacklist_bed.toString())
    assert bl.exists() : "peak_blacklist_bed not found: ${params.peak_blacklist_bed}"
  }

  def rows
  if (params.macs3_samplesheet) {
    rows = Channel
      .fromPath(params.macs3_samplesheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        assert row.sample_id && row.treatment_bam : "macs3_samplesheet must contain: sample_id,treatment_bam"
        tuple(
          row.sample_id.toString().trim(),
          file(row.treatment_bam.toString().trim()),
          row.control_bam?.toString()?.trim() ?: ''
        )
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

    def bamDir = resolveBaseDir(params.chipfilter_output)
    assert bamDir.exists() : "chipfilter_output directory not found: ${params.chipfilter_output}"

    def resolveCleanBam = { sid ->
      def hits = bamDir.listFiles()?.findAll { f ->
        f.isFile() && f.name.endsWith('.nomulti.bam') && f.name.startsWith("${sid}")
      } ?: []

      if (hits.isEmpty()) {
        throw new IllegalArgumentException("No MAPQ-filtered BAM found for sample_id '${sid}' under: ${params.chipfilter_output}")
      }
      if (hits.size() > 1) {
        throw new IllegalArgumentException("Multiple MAPQ-filtered BAM files matched sample_id '${sid}': ${hits*.name.join(', ')}")
      }
      file(hits[0].toString())
    }

    def autoRows = enabledRecords.findAll { rec -> isChip(rec) }.collect { rec ->
      def sid = rec.sample_id?.toString()?.trim()
      if (!sid) return null

      def treatBam = resolveCleanBam(sid)
      def cid = rec.control_id?.toString()?.trim()
      if (!cid && defaultControl) cid = defaultControl
      def controlBam = ''
      if (cid) {
        assert controlSet.contains(cid) : "control_id '${cid}' for sample '${sid}' is not an enabled control sample in samples_master"
        controlBam = resolveCleanBam(cid).toString()
      }
      tuple(sid, treatBam, controlBam)
    }.findAll { it != null }

    rows = Channel
      .fromList(autoRows)
      .ifEmpty { exit 1, "ERROR: No treatment rows generated from samples_master: ${params.samples_master}" }
  } else {
    exit 1, "ERROR: Provide --macs3_samplesheet or --samples_master."
  }

  def jobs = rows.flatMap { sample_id, treat_bam, control_bam ->
    profiles.collect { prof ->
      tuple(prof.profile_name, prof.qval, sample_id, treat_bam, control_bam)
    }
  }
  .filter { profile_name, qval, sample_id, treat_bam, control_bam ->
    def profDir = file("${outBase}/${profile_name}")
    def peakOk = file("${profDir}/${sample_id}_peaks.${peak_ext}").exists()
    def xlsOk = file("${profDir}/${sample_id}_peaks.xls").exists()
    def blOk = !params.peak_blacklist_bed || file("${profDir}/${sample_id}_peaks.blacklist_applied.txt").exists()
    !(peakOk && xlsOk && blOk)
  }

  def called = macs3_callpeak(jobs)
  if (params.peak_blacklist_bed) {
    def blPath = file(params.peak_blacklist_bed.toString())
    def blCh = Channel.value(blPath)
    def blInputs = called.peaks.combine(blCh).map { profile_name, sample_id, peak_file, blacklist_bed ->
      tuple(profile_name, sample_id, peak_file, blacklist_bed)
    }
    macs3_blacklist_peak(blInputs)
  }
}
