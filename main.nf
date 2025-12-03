#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process macs3 {
  tag "${sample}.Rep_${replicate}"
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val sample
    val replicate
    val treatment_file
    val control_file
    val peak_type

  script:
  """
  set -eux

  mkdir -p ${params.project_folder}/macs3_output
  mkdir -p ${params.project_folder}/tmp

  sample_rep="${sample}.Rep_${replicate}"

  cd ${params.project_folder}/bowtie2_output

  if [[ ${params.seq} == "paired" ]] && [[ ${params.input} == 'no' ]] ; then

    echo "macs3 callpeak -t ${treatment_file} -f BAMPE -q 0.05 --keep-dup all --gsize ${params.GeTAG} --outdir ${params.project_folder}/macs3_output -n \${sample_rep} -B --tempdir=${params.project_folder}/tmp ${peak_type}"

    macs3 callpeak \
      -t ${treatment_file} \
      -f BAMPE \
      -q 0.05 \
      --keep-dup all \
      --gsize ${params.GeTAG} \
      --outdir ${params.project_folder}/macs3_output \
      -n \${sample_rep} \
      -B \
      --tempdir=${params.project_folder}/tmp \
      ${peak_type}

  elif [[ ${params.seq} == "single" ]] && [[ ${params.input} == 'no' ]] ; then

    echo "macs3 callpeak -t ${treatment_file} -f BAM -q 0.05 --keep-dup all --gsize ${params.GeTAG} --outdir ${params.project_folder}/macs3_output -n \${sample_rep} -B --tempdir=${params.project_folder}/tmp ${peak_type}"

    macs3 callpeak \
      -t ${treatment_file} \
      -f BAM \
      -q 0.05 \
      --keep-dup all \
      --gsize ${params.GeTAG} \
      --outdir ${params.project_folder}/macs3_output \
      -n \${sample_rep} \
      -B \
      --tempdir=${params.project_folder}/tmp \
      ${peak_type}

  elif [[ ${params.seq} == "paired" ]] && [[ ${params.input} == 'yes' ]] ; then

    echo "macs3 callpeak -t ${treatment_file} -c ${control_file} -f BAMPE -q 0.05 --keep-dup all --gsize ${params.GeTAG} --outdir ${params.project_folder}/macs3_output -n \${sample_rep} -B --tempdir=${params.project_folder}/tmp ${peak_type}"

    macs3 callpeak \
      -t ${treatment_file} \
      -c ${control_file} \
      -f BAMPE \
      -q 0.05 \
      --keep-dup all \
      --gsize ${params.GeTAG} \
      --outdir ${params.project_folder}/macs3_output \
      -n \${sample_rep} \
      -B \
      --tempdir=${params.project_folder}/tmp \
      ${peak_type}

  elif [[ ${params.seq} == "single" ]] && [[ ${params.input} == 'yes' ]] ; then

    echo "macs3 callpeak -t ${treatment_file} -c ${control_file} -f BAM -q 0.05 --keep-dup all --gsize ${params.GeTAG} --outdir ${params.project_folder}/macs3_output -n \${sample_rep} -B --tempdir=${params.project_folder}/tmp ${peak_type}"

    macs3 callpeak \
      -t ${treatment_file} \
      -c ${control_file} \
      -f BAM \
      -q 0.05 \
      --keep-dup all \
      --gsize ${params.GeTAG} \
      --outdir ${params.project_folder}/macs3_output \
      -n \${sample_rep} \
      -B \
      --tempdir=${params.project_folder}/tmp \
      ${peak_type}

  fi
  """
}


workflow {

  if ( params.peak_type == "none" ) {
    peak_type = ""
  } else {
    peak_type = params.peak_type
  }

  if ( params.input == "yes" ) {
    file_ending = ".chip.md.bam"
  } else {
    file_ending = ".md.bam"
  }

  rows = Channel
    .fromPath("${params.samples_csv}", checkIfExists: true)
    .splitCsv(sep: ',', skip: 1)

  rows = rows.filter { n ->
    ! file("${params.project_folder}/macs3_output/${n[1]}.Rep_${n[2]}_peaks.xls").exists()
  }

  sample        = rows.flatMap { n -> n[1] }
  replicate     = rows.flatMap { n -> n[2] }
  treatment_file= rows.flatMap { n -> n[3] }
  control_file  = rows.flatMap { n -> n[5] }

  macs3(sample, replicate, treatment_file, control_file, peak_type)
}
