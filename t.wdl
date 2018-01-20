workflow TOPMed {
  String sample
  String sample_bam_pair_name1
  String sample_bam_pair_name2
  File samtools
  File human_ref_fasta

  call sort_cram_file {
  input: 
      samtools_path=samtools,
      sampleName=sample
  }

  call get_header {
  input: 
      samtools_path=samtools,
      sorted_cram=sort_cram_file.outputCRAM,
      sampleName=sample
  }

  call get_fastq {
  input: 
      samtools_path=samtools,
      sorted_cram=sort_cram_file.outputCRAM,
      sampleName=sample
  }

  call perform_alignment as align_1 {
  input: 
      sampleName=sample_bam_pair_name1,
      read_group_header=get_header.CRAMheader,
      reference_fa=human_ref_fasta,
      fq_file=get_fastq.CRAMfastq_1
  }

  call perform_alignment as align_2 {
  input: 
      sampleName=sample_bam_pair_name2,
      read_group_header=get_header.CRAMheader,
      reference_fa=human_ref_fasta,
      fq_file=get_fastq.CRAMfastq_2
  }

  call merge_bam_files {
  input: 
      samtools_path=samtools,
      sampleName=sample,
      BAM1=align_1.aligned_bam,
      BAM2=align_2.aligned_bam
  }

  call sort_bam_file {
  input: 
      inputBAM=merge_bam_files.merged_bam,
      samtools_path=samtools,
      sampleName=sample
  }

  call mark_duplicates {
  input: 
      sampleName=sample,
      duped_bam_input=sort_bam_file.outputBAM
  }

  call recalibrate_quality_scores {
  input:
      ref_fasta=human_ref_fasta,
      deduped_marked_bam=mark_duplicates.deduped_bam_output
  }
}

task sort_cram_file {
  File samtools_path
  File inputCRAM
  String sampleName

  command {
    ${samtools_path} sort ${inputCRAM} -T tmp.srt -O cram -o ${sampleName}_sorted.cram
  }

  output {
    File outputCRAM = "${sampleName}_sorted.cram"
  }
}

task get_header {
  File samtools_path
  File sorted_cram
  String sampleName

  command {
    ${samtools_path} view -H ${sorted_cram} > ${sampleName}_sorted.cram.header
  }

  output {
    File CRAMheader = "${sampleName}_sorted.cram.header"
  }
}

task get_fastq {
  File samtools_path
  File sorted_cram
  String sampleName

  command {
    ${samtools_path} fastq -1 ${sampleName}_1.fq -2 ${sampleName}_2.fq -0 ${sampleName}_unpaired.fq ${sorted_cram}
  }

  output {
    File CRAMfastq_1 = "${sampleName}_1.fq"
    File CRAMfastq_2 = "${sampleName}_2.fq"
    File CRAMfastq_u = "${sampleName}_unpaired.fq"
  }
}

task perform_alignment {
  File read_group_header
  File reference_fa
  File fq_file
  String sampleName

  command {
    /opt/bwa.kit/bwa mem -K 100000000 -R ${read_group_header} -Y ${reference_fa} ${fq_file} | /opt/samblaster-v.0.1.24/samblaster -a --addMateTags | /opt/bwa.kit/samtools view -Sb > ${sampleName}.out.bam
  }

  output {
    File aligned_bam = "${sampleName}.out.bam"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/bwakit:0.7.15--ed1aeaaebf4d88ba51042da83f02ef8373a822b9'
  }
}

task merge_bam_files {
  File samtools_path
  File BAM1
  File BAM2
  String sampleName

  command {
  ${samtools_path} merge -c -p ${sampleName}_merged.bam ${BAM1} ${BAM2}
  }

  output {
    File merged_bam = "${sampleName}_merged.bam"
  }
}

task sort_bam_file {
  File samtools_path
  File inputBAM
  String sampleName

  command {
    ${samtools_path} sort ${inputBAM} -T tmp.srt -O bam -o ${sampleName}_sorted.bam
  }

  output {
    File outputBAM = "${sampleName}_sorted.bam"
  }
}

task mark_duplicates {
  File duped_bam_input
  String sampleName
  
  command {
    java -jar /home/lifeisaboutfishtacos/Desktop/toil-workflows/tools/picard.jar MarkDuplicates I=${duped_bam_input} O=${sampleName}.deduped.bam M=deduplicated_metrics.txt
  }

  output {
    File deduped_bam_output = "${sampleName}.deduped.bam"
  }
}

task recalibrate_quality_scores {
  File GATK
  File ref_fasta
  File ref_index
  File ref_dict
  File deduped_marked_bam
  String recalibration_report_name
  File known_site_1
  File known_site_2
  File known_site_3

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
    File x = "${recalibration_report_name}"
  }

  runtime {
    docker: 'quay.io/ucsc_cgl/gatk:3.7--e931e1ca12f6d930b755e5ac9c0c4ca266370b7b'
  }
}