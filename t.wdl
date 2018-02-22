## Copyright UCSC, 2018
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.


workflow TOPMed {
  String sample
  String sample_bam_pair_name1
  String sample_bam_pair_name2
  File human_ref_fasta
  File human_ref_index
  File human_ref_index2
  File human_ref_dict
  File initial_bam_sample
  File bwa_sa
  File bwa_pac
  File bwa_bwt
  File bwa_ann
  File bwa_amb
  File index_folder

  call o0_init_sort_bam_file {
  input: 
      sampleName=sample,
      inputBAM=initial_bam_sample
  }

  call o1_get_header {
  input: 
      sorted_bam=o0_init_sort_bam_file.outputBAM,
      sampleName=sample
  }

  call o1_get_fastq {
  input: 
      sorted_bam=o0_init_sort_bam_file.outputBAM,
      sampleName=sample
  }

  call o2_perform_alignment {
  input: 
      sorted_BAM1=o1_get_fastq.BAMfastq_1,
      sorted_BAM2=o1_get_fastq.BAMfastq_2,
      sampleName=sample_bam_pair_name1,
      read_group_header=o1_get_header.BAMheader,
      align_sa=bwa_sa,
      align_pac=bwa_pac,
      align_bwt=bwa_bwt,
      align_ann=bwa_ann,
      align_amb=bwa_amb,
      reference_fa=human_ref_fasta,
      reference_index=human_ref_index,
      reference_index2=human_ref_index2,
      reference_dict=human_ref_dict,
      index_path=index_folder
  }

  call o3_mark_duplicates {
  input: 
      sampleName=sample,
      duped_bam_input=o2_perform_alignment.aligned_bam
  }

  call o4_sort_bam {
  input:
      sampleName=sample,
      input_bamfile=o3_mark_duplicates.deduped_bam_output
  }

  call o5_recalibrate_quality_scores {
  input:
      ref_fasta=human_ref_fasta,
      ref_index=human_ref_index,
      ref_index2=human_ref_index2,
      ref_dict=human_ref_dict,
      deduped_marked_bam=o4_sort_bam.indexsorted_bam_output
  }

  call o6_bin_quality_scores {
  input:
      ref_fasta=human_ref_fasta,
      ref_index=human_ref_index,
      ref_index2=human_ref_index2,
      ref_dict=human_ref_dict,
      input_bam_to_bin=o5_recalibrate_quality_scores.recal_bam_output,
      binned_out_name=sample
  }

  call o7_convert_to_cram {
  input:
      ref_fasta=human_ref_fasta,
      ref_index=human_ref_index,
      ref_index2=human_ref_index2,
      ref_dict=human_ref_dict,
      input_end_bam=o6_bin_quality_scores.binned_output_bam,
      output_cram_name=sample
  }
}

task o0_init_sort_bam_file {
  File inputBAM
  String sampleName

  command {
    samtools sort ${inputBAM} -T tmp.srt -O bam -o ${sampleName}_sorted.bam
  }

  output {
    File outputBAM = "${sampleName}_sorted.bam"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c'
  }
}

task o1_get_header {
  File sorted_bam
  String sampleName

  command {
    samtools view -H ${sorted_bam} > ${sampleName}_sorted.bam.header
  }

  output {
    File BAMheader = "${sampleName}_sorted.bam.header"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c'
  }
}

task o1_get_fastq {
  File sorted_bam
  String sampleName

  command {
    samtools fastq -1 ${sampleName}_1.fq -2 ${sampleName}_2.fq -0 ${sampleName}_unpaired.fq ${sorted_bam}
  }

  output {
    File BAMfastq_1 = "${sampleName}_1.fq"
    File BAMfastq_2 = "${sampleName}_2.fq"
    File BAMfastq_u = "${sampleName}_unpaired.fq"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c'
  }
}

task o2_perform_alignment {
  File sorted_BAM1
  File sorted_BAM2
  File read_group_header
  String header = read_string(read_group_header)
  File reference_fa
  File reference_index
  File reference_index2
  File reference_dict
  String sampleName
  File align_sa
  File align_pac
  File align_bwt
  File align_ann
  File align_amb
  File seqtk_path
  File bwa_path
  File samblaster_path
  File samtools_path
  File index_path

  command {
    ${seqtk_path} mergepe ${sorted_BAM1} ${sorted_BAM2} | ${bwa_path} mem -p -R '@RG\tID:B\tSM:NA12878\tLB:Solexa-NA12878\tPL:illumina\tPU:H06HDADXX130110.1.ATCACGAT\n@HD\tVN:1.5\tSO:coordinate' ${reference_fa} - 2> outx.log.bwamem | ${samblaster_path} --addMateTags 2> outx.log.dedup | ${samtools_path} sort - ${sampleName}_aligned.out
  }

  output {
    File aligned_bam = "${sampleName}_aligned.out.bam"
  }
}

task o3_mark_duplicates {
  File duped_bam_input
  String sampleName
  
  command {
    java -jar /opt/picard-tools/picard.jar MarkDuplicates I=${duped_bam_input} O=${sampleName}.deduped.bam M=deduplicated_metrics.txt
  }

  output {
    File deduped_bam_output = "${sampleName}.deduped.bam"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/picardtools:2.9.0--4d726c4a1386d4252a0fc72d49b1d3f5b50b1e23'
  }
}

task o4_sort_bam {
  File input_bamfile
  String sampleName

  command {
    java -jar /opt/picard-tools/picard.jar SortSam I=${input_bamfile} O=${sampleName}.indexsorted.sam SORT_ORDER=coordinate
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/picardtools:2.9.0--4d726c4a1386d4252a0fc72d49b1d3f5b50b1e23'
  }

  output {
    File indexsorted_bam_output = "${sampleName}.indexsorted.sam"
  }
}

task o5_recalibrate_quality_scores {
  File ref_fasta
  File ref_index
  File ref_index2
  File ref_dict
  File deduped_marked_bam
  String recalibration_report_name
  File known_site_1
  File known_site_2
  File known_site_3
  File known_site_1_index
  File known_site_2_index
  File known_site_3_index

  command {
    java -jar /opt/gatk/gatk.jar -T BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${deduped_marked_bam} \
    -o ${recalibration_report_name} \
    -knownSites "${known_site_1}" \
    -knownSites "${known_site_2}" \
    -knownSites "${known_site_3}" \
    --downsample_to_fraction .1 -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 \
    --interval_padding 200 -rf BadCigar \
    --preserve_qscores_less_than 6 \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    --disable_bam_indexing \
    --useOriginalQualities
  }

  output {
    File recal_bam_output = "${recalibration_report_name}.recal.bam"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/gatk:3.7--e931e1ca12f6d930b755e5ac9c0c4ca266370b7b'
  }
}

task o6_bin_quality_scores {
  File ref_fasta
  File ref_index
  File ref_index2
  File ref_dict
  File input_bam_to_bin
  String binned_out_name
  command {
  java -jar /opt/gatk/gatk.jar -T PrintReads -R ${ref_fasta} \
                   -I ${input_bam_to_bin} \
                   -o ${binned_out_name}.binned.bam \
                   -BQSR /data/bqsr_report_name \
                   -SQQ 10 -SQQ 20 -SQQ 30 \
                   --globalQScorePrior -1.0 \
                   --preserve_qscores_less_than 6 \
                   --disable_indel_quals \
                   --useOriginalQualities \
                   -rf BadCigar
  }
  output {
    File binned_output_bam = "${binned_out_name}.binned.bam"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/gatk:3.7--e931e1ca12f6d930b755e5ac9c0c4ca266370b7b'
  }
}

task o7_convert_to_cram {
  File ref_fasta
  File ref_index
  File ref_index2
  File ref_dict
  File input_end_bam
  String output_cram_name

  command {
    samtools view -C -T ${ref_fasta} -o ${output_cram_name} ${input_end_bam}
  }

  output {
    File final_cram = "${output_cram_name}.final.cram"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c'
  }
}
