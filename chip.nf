/*
 * 'ChIP_nf' - A Nextflow pipeline for ChIP-seq data analysis
 *
 * This pipeline deals with ChIP-seq data 
 * input is reads in FASTQ format
 *
 * Feng Yan
 * feng.yan@monash.edu
 */


/*
 * Define the default parameters
 */

params.reads      = "$baseDir/data/sample_R{1,2}.fq.gz"
params.outdir     = "results"
params.aligner    = "bwa_mem"
params.species    = "mm10"
params.samplesheet= "$baseDir/data/samplesheet.csv"
params.trim       = "trim"
params.fasta      = ""
params.blacklist  = "mm10.blacklist.bed"

log.info """\
R R B S -  N F    v 1.0
================================
species  	: $params.species
reads    	: $params.reads
outdir   	: $params.outdir
samplesheet	: $params.samplesheet
trim            : $params.trim
fasta           : $params.fasta
blacklist       : $params.blacklist
"""

/*
 *  Parse the input parameters
 */

Channel
        .fromPath( params.fasta )
        .into{ fasta_ch; fasta_index_ch; ch_bam_filter }

samplesheet     = file(params.samplesheet)
species         = Channel.from(params.species)
blacklist       = file(params.blacklist)


/*
 * PART 0: Preparation
 */
process '0A_get_software_versions' {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    executor 'local'

    output:
    file '*.txt'

    script:
    """
    echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt
    """
}

/*
 * Create a channel for input read files
 */
Channel
        .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { raw_reads_fastqc_ch; raw_reads_align_ch }


/**********
 * PART 1: Preprocessing
 *
 * Process 1A: fastqc report for raw data
 */
process '1A_pre_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/pre_fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from raw_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/*
 * Process 1B: 2-step triming for NuGen RRBS data
 */
//process '1B_trim' {
//    tag "$name"
//    label 'big'
//    publishDir "${params.outdir}/trim"
//
//    input:
//    set val(name), file(reads) from raw_reads_trim_ch
//    file(trimpy) from numet
//    // val(clip3) from params.clip3prime
//
//    output:
//    set val(name), file('*.fq_trimmed.fq.gz') into clean_reads_bismark_ch, clean_reads_fastqc_ch
//    file('*report.txt') into ch_trimgalore_results_for_multiqc optional true
//    file('*.log') optional true 
//     
//    script:
//    if( params.library == "nugen" ) {
//    if( params.single_end ) {
//            """
//            trim_galore -a AGATCGGAAGAGC $reads --cores ${task.cpus}
//            echo $trimpy
//            python2 $trimpy -1 ${reads.simpleName}_trimmed.fq.gz &> ${reads.simpleName}_trimpy.log
//            """
//        } else {
//            """
//            trim_galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC \\
//            --paired $reads --cores ${task.cpus}
//            echo $trimpy
//            python2 $trimpy -1 ${reads[0].simpleName}_val_1.fq.gz \\
//            -2 ${reads[1].simpleName}_val_2.fq.gz &> ${name}_trimpy.log
//            """
//        }
//    } else if (params.library == "epic" && params.trim == "trim") { // adde trim for Epic becasue we noticed bias toward 3 prime end
//    if( params.single_end ) {
//            // clip3 = params.clip3prime == 0 ? "" : "--three_prime_clip_R1 $params.clip3prime"
//            """
//            ## leave to auto detection
//            trim_galore $reads --cores ${task.cpus}
//            mv ${reads.simpleName}_trimmed.fq.gz ${reads.simpleName}.fq_trimmed.fq.gz
//            """
//        } else {
//            // clip3 = params.clip3prime == 0 ? "" : "--three_prime_clip_R1 $params.clip3prime --three_prime_clip_R2 $params.clip3prime"
//            """
//            trim_galore --paired $reads --cores ${task.cpus}
//            mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}.fq_trimmed.fq.gz 
//            mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}.fq_trimmed.fq.gz
//            """
//        }
//    } else if (params.library == "epic" && params.trim == "skip") {
//    if( params.single_end ) {
//            """
//            mv ${reads} ${reads.simpleName}.fq_trimmed.fq.gz
//            """
//        } else {
//            """
//            mv ${reads[0]} ${reads[0].simpleName}.fq_trimmed.fq.gz
//            mv ${reads[1]} ${reads[1].simpleName}.fq_trimmed.fq.gz
//            """
//        }
//    }
//}

/**********
 * Process 1C: fastqc report for trimmed data
 */
//process '1C_post_fastqc' {
//    tag "$name"
//    label 'big'
//    publishDir "${params.outdir}/post_fastqc", mode: 'copy',
//        saveAs: { filename ->
//                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
//                }
//
//    input:
//    set val(name), file(reads) from clean_reads_fastqc_ch
//
//    output:
//    file '*_fastqc.{zip,html}' into ch_fastqc2_results_for_multiqc
//     
//    when:
//    params.trim == "trim"
// 
//    script:
//    """
//    fastqc --quiet --threads ${task.cpus} $reads
//    """
//}


/**********
 * PART 2: Alignment
 *
 * Process 2A: Create a genome index (.fai) 
 */

process '2A_index_genome' {
  tag "$fasta.baseName"
  label 'bismark'

  input:
      file fasta from fasta_ch
  output:
      file 'bwa_index' into ch_bwa_index

  """
  bwa index $fasta 
  mkdir bwa_index && mv ${fasta}* bwa_index
  """
}


//raw_reads_align_ch.view()
lastPath = params.fasta.lastIndexOf(File.separator)
bwa_base = params.fasta.substring(lastPath+1)

/**********
 * Process 2B: Align to the genome
 */

process '2B_mapping' {
  tag "$name"
  label 'bismark'
  publishDir "${params.outdir}/rawbamfiles", mode: 'symlink'

  input:
      //file index from ch_bwa_index
      set val(name), file(reads), file(index) from raw_reads_align_ch.combine(ch_bwa_index)
      
  output:
      set val(name), file('*.bam'), file('*.bai') into ch_bwa_bam
      //set val(name), file('*_dups.txt'), file('*_insert.txt'), file('*.flagstat'), file('*.idxstats') into ch_bamqc_for_multiqc, ch_bismark_align_log_for_Rsummary

  script:
  """
  bwa mem -t ${task.cpus} ${index}/${bwa_base} $reads | samtools view -@ ${task.cpus} -Sb - | samtools sort -@ ${task.cpus} - > ${name}_sorted.bam 
  samtools index -@ ${task.cpus} ${name}_sorted.bam
  """
}



/**********
 * Process 2C: Post-alignment processing BAM files
 */

process '2B_filter_bam' {
  tag "$name"
  label 'bismark'
  publishDir "${params.outdir}/bamfiles", mode: 'copy'

  input:
      set val(name), file(bam), file(bai), file(fasta) from ch_bwa_bam.combine(ch_bam_filter)

  output:
      set val(name), file('*final.bam'), file('*final.bam.bai') into ddup_bam_ch
      set val(name), file("${name}.flagstat"), file("${name}.idxstats"), file("${name}_dups.txt"), file("${name}_alignmetrics.txt"), file("${name}.final.flagstat") into ch_bamqc_for_multiqc
      file('*.pdf') optional true 
      set val(name), file("${name}_insert.txt") optional true into ch_insert_multiqc
      // Double-quoted strings support variable interpolations, while single-quoted strings do not.

  script:
  flag = params.single_end ? "" : "-f 2"
  """
  ## filter 
  ## samtools view -h -F 1804 ${flag} -q 30 ${bam} > ${name}.filtered.bam
   
  picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true \
                        INPUT=${bam} OUTPUT=${name}_sorted_mdups.bam \
                        METRICS_FILE=${name}_dups.txt

  picard CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=LENIENT \
                                        REFERENCE_SEQUENCE=${fasta} \
                                        INPUT=${bam} \
                                        OUTPUT=${name}_alignmetrics.txt

  picard CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT \
                                  INPUT=${bam} \
                                  OUTPUT=${name}_insert.txt \
                                  HISTOGRAM_FILE=${name}_insert_hist.pdf \
                                  M=0.5

  samtools flagstat ${bam} > ${name}.flagstat
  samtools idxstats ${bam} > ${name}.idxstats

  # PBC File output
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
  
  ## remove duplicates and final bam
  samtools view -@ ${task.cpus} -h -F 1804 ${flag} -q 30 -Sb ${name}_sorted_mdups.bam > ${name}.final.bam
  samtools flagstat ${name}.final.bam > ${name}.final.flagstat
  samtools index -@ ${task.cpus} ${name}.final.bam
  """
}

// bedgraph_bismark_ch.view()

/**********
 * PART 3: Summary
 *
 * Process 3A: MultiQC
 */
process '3A_multiqc' {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    // only use fastqc from trimmed reads
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('bamfiles/*') from ch_bamqc_for_multiqc.collect().ifEmpty([])
    file ('bamfiles/_insert.txt') from ch_insert_multiqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html"
    file "*_data"

    script:
    """
    multiqc -f .
    """
}

/**********
 * PART 4: Visualization
 *
 * Process 4A: Generate fasta index
 */
process '4A_faidx' {
    tag "$fasta.baseName"
    label 'big'

    input:
    file fasta from fasta_index_ch

    output:
    file "chrom.sizes" into chr_size_ch

    script:
    if( fasta.extension ==~ /fa|fasta/ ) {
            """
            samtools faidx ${fasta}
            cut -f1,2 ${fasta}.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       } else if( fasta.extension == 'gz' ) {
	    """
            zcat ${fasta} | bgzip -c > ${fasta.simpleName}.fa.bgz
            samtools faidx ${fasta.simpleName}.fa.bgz
            cut -f1,2 ${fasta.simpleName}.fa.bgz.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       }
}

// chr_size_ch.view()
// bedgraph_bismark_ch.view()
// bedgraph_bismark_ch.combine(chr_size_ch).view()

/**********
 * Process 4B: Generate bigwig files
 */
process '4B_BAMtoBigWig' {
    tag "$name"
    label 'bismark'
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    set val(name), file(bam), file(bai) from ddup_bam_ch

    output:
    file "*.bw"

    script:
    effectiveGenomeSize = params.species == 'mm10' ? '2652783500' : '2913022398'
    // for mm10 and hg38 currently
    blacklist = params.blacklist == '' ? '' : "--blackListFileName ${blacklist}"
    """
    bamCoverage -b ${bam} -o ${name}.bw -p ${task.cpus} --normalizeUsing RPGC --effectiveGenomeSize $effectiveGenomeSize
 
    """
}

//covgz_for_Rsummary
//  .join(ch_bismark_align_log_for_Rsummary).collect().view()

///**********
// * Process 4C: Generate summary statistics
// */
//process '4c_toRSummary' {
//    tag "summaryplot"
//    label 'bismark'
//    publishDir "${params.outdir}/summaryplot", mode: 'copy'
//
//    input:
//    file("*") from covgz_for_Rsummary.join(ch_bismark_align_log_for_Rsummary).collect()
//    file(samplesheet) from samplesheet
//    val(species) from species
//    file(summary) from summary
//    
//    output:
//    file "*.png"
//    // file "*.RData"
//
//    script:
//    """
//    module load R
//    Rscript --vanilla $summary $samplesheet $species
//    """
//}
