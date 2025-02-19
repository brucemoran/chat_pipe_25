nextflow.enable.dsl=2

/*
 * Example Nextflow DSL2 pipeline for:
 *   1) Reading single/paired FASTQ paths from a CSV
 *   2) Running fastp for QC
 *   3) Aligning to GRCh38 with bwa mem
 *   4) Mark duplicates and Base Quality Score Recalibration (BQSR)
 *   5) Variant calling with Mutect2
 *   6) Annotating final VCF with PCGR
 *
 * Make sure you have Nextflow installed and Docker (or another container engine)
 * configured so that Nextflow can pull the specified images. 
 */

/////////////////////////////////////////////////////
// Pipeline parameters
/////////////////////////////////////////////////////
params.samples_csv = (params.samples_csv ?: 'samples.csv')
params.outdir      = (params.outdir ?: 'results')

// URL and filename for the GRCh38 reference genome. In practice, you can
// supply your own pre-downloaded reference.
params.genome_url  = 'https://storage.googleapis.com/genomics-public-data/references/hg38/v0'
//params.genome_base = 'GRCh38.fa'
//make the real name, use a base and download everything with that base
params.genome_base = 'Homo_sapiens_assembly38'

/////////////////////////////////////////////////////
// Processes
/////////////////////////////////////////////////////

//workflow for reference
workflow reference {
    take:
    genome_url_ch
    genome_base_ch

    main:
    CURL_REFERENCE(genome_url_ch, genome_base_ch)
    INDEX_REFERENCE(CURL_REFERENCE.out)
    GATK_DICT_REFERENCE(CURL_REFERENCE.out)

    emit:
    ref_fa = INDEX_REFERENCE.out[0]
    ref_fai = INDEX_REFERENCE.out[1]
    ref_dict = GATK_DICT_REFERENCE.out
}

// Download the reference genome if it does not exist, then index it
process CURL_REFERENCE {
    tag 'download_grch38'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val genome_url
    val genome_base

    output:
    path "${genome_base}.fasta", emit: ref_fa

    script:
    """
    ## get fasta
    echo "Downloading: ${genome_base}.fasta"
    curl -L -o ${genome_base}.fasta ${genome_url}/${genome_base}.fasta
    """
}

process INDEX_REFERENCE {
    tag 'download_grch38'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref_fa

    output:
    path "${ref_fa}", emit: ref_fa
    path "${ref_fa}.fai", emit: ref_fai

    script:
    """
    ## get fasta
    echo "Indexing: ${ref_fa}"
    samtools faidx ${ref_fa}
    """
}

process GATK_DICT_REFERENCE {
    memory '20 GB'
    tag 'download_grch38'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref_fa

    output:
    path "${ref_fa}.dict", emit: ref_dict

    script:
    """
    ## get fasta
    echo "CreateSequenceDictionary: ${ref_fa}"
    gatk --java-options "-Xmx20G" CreateSequenceDictionary -R ${ref_fa} -O ${ref_fa}.dict
    """
}

// Run fastp QC on single- or paired-end FASTQs
process FASTP_QC {
    tag { sample_id }
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("cleaned_*.fastq.gz")

    script:
    """
    ## If fq2 is empty, run single-end; otherwise, run paired-end
    if [ -f ${fq2} ]; then
        fastp \\
            --in1 ${fq1} --in2 ${fq2} \\
            --out1 cleaned_R1.fastq.gz --out2 cleaned_R2.fastq.gz \\
            --thread 4
    else
        fastp \\
            --in1 ${fq1} \\
            --out1 cleaned_R1.fastq.gz \\
            --thread 4
    fi
    """
}

// Align reads with bwa mem and produce a sorted BAM
process ALIGN_BWA {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_reads)
    path ref_fa   
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.bam")
    path ref_fa
    path ref_fai

    script:
    // cleaned_reads might have 1 or 2 fastq files
    def r1 = cleaned_reads.find { it.name.contains('R1') }
    def r2 = cleaned_reads.find { it.name.contains('R2') }
    """
    # index fasta
    bwa index ${ref_fa}

    if [ -f "${r2}" ]; then
        bwa mem -t 4 ${ref_fa} ${r1} ${r2} > ${sample_id}.bam
    else
        bwa mem -t 4 ${ref_fa} ${r1} > ${sample_id}.bam
    fi
    """
}

process INDEX_BAM {
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")
    path ref_fa
    path ref_fai

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

// Mark duplicates & perform Base Quality Score Recalibration
process GATK_MARKDUP {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam")
    path ref_fa
    path ref_fai

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.dedup.bam \\
        -M ${sample_id}.dedup.metrics.txt \\
        --VALIDATION_STRINGENCY SILENT \\
        --REMOVE_SEQUENCING_DUPLICATES false
    """
}

// Mark duplicates & perform Base Quality Score Recalibration
process INDEX_MARKDUP {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam")
    path ref_fa
    path ref_fai

    script:
    """
    samtools index ${sample_id}.dedup.bam
    """
}

// Mark duplicates & perform Base Quality Score Recalibration
process GATK_BQSR {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.dedup.recal.bam")
    path ref_fa
    path ref_fai

    script:
    def outBamPrefix = "${sample_id}.dedup.recal"
    """
    ln -s ${fa_dict} \$(basename ${fa_dict} | sed 's/fasta.//')
    # Example known-sites from GATK resources (b37 used for demonstration).
    # Adjust for GRCh38 best-practice known sites in a real pipeline
    gatk BaseRecalibrator \\
        -R ${ref_fa} \\
        -I ${sample_id}.dedup.bam \\
        --known-sites https://storage.googleapis.com/gatk-best-practices/somatic-b37/1000G_phase1.indels.b37.vcf.gz \\
        --known-sites https://storage.googleapis.com/gatk-best-practices/somatic-b37/dbsnp_138.b37.vcf.gz \\
        -O ${sample_id}.recal_data.table

    gatk ApplyBQSR \\
        -R ${ref_fa} \\
        -I ${sample_id}.dedup.bam \\
        --bqsr-recal-file ${sample_id}.recal_data.table \\
        -O ${outBamPrefix}.bam
    """
}

// Mark duplicates & perform Base Quality Score Recalibration
process INDEX_BQSR {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.dedup.recal.bam")
    path ref_fa
    path ref_fai

    script:
    def outBamPrefix = "${sample_id}.dedup.recal"
    """
    samtools index ${outBamPrefix}.bam
    """
}

// Call somatic variants with GATK Mutect2
process GATK_MUTECT2_CALL {
    tag { sample_id }
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path ref_fa
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.mutect2.vcf.gz")

    script:
    """
    gatk Mutect2 \\
        -R ${ref_fa} \\
        -I ${bam} \\
        -tumor ${sample_id} \\
        -O ${sample_id}.mutect2.vcf.gz
    """
}

// Annotate the final VCF with PCGR
process PCGR_ANNOTATE {
    tag { sample_id }
    publishDir "${params.outdir}/pcgr", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    // All PCGR outputs typically go into a directory named <sample_id>
    // plus an annotated VCF. Adjust as needed.
    path "${sample_id}_out/*"

    script:
    def outPrefix = "${sample_id}.pcgr_anno"
    """
    BUNDLE="https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20240927.grch38.tgz"
    wget \${BUNDLE}
    gzip -dc *tgz | tar -xvf -
    pcgr --input_vcf ${vcf} \\
         --genome_assembly grch38 \\
         --pcgr_dir data \\
         --sample_id ${sample_id} \\
         --output_dir ${sample_id}_out \\
         --output_vcf ${outPrefix}.vcf \\
         --force_overwrite
    """
}

/////////////////////////////////////////////////////
// The main pipeline workflow
/////////////////////////////////////////////////////

workflow {
    // Download reference files (runs once)
    //take: 
    //download_reference = DOWNLOAD_REFERENCE(genome_url_ch, genome_base_ch)
    // Simple channels for the reference parameters
    genome_url_ch = Channel.value(params.genome_url)

    genome_base_ch = Channel.value(params.genome_base)

    reference(genome_url_ch, genome_base_ch)

    // Read the CSV and create a channel of [sample_id, fastq1, fastq2]
    sample_info_ch = Channel
                        .fromPath(params.samples_csv)
                        .splitCsv(header:true)
                        .map { row ->
                            // For single-end data, fastq_2 might be empty
                            def sample_id = row.sample_id
                            def fq1       = row.fastq_1
                            def fq2       = row.fastq_2
                            return [ sample_id, fq1, fq2 ]
                        }
    // Now process each sample in turn
    //sample_info_ch.view()
    FASTP_QC(sample_info_ch)
    ALIGN_BWA(FASTP_QC.out,
              reference.out.ref_fa, 
              reference.out.ref_fai,
              reference.out.ref_dict)
    INDEX_BAM(ALIGN_BWA.out)
    GATK_MARKDUP(INDEX_BAM.out)
    INDEX_MARKDUP(GATK_MARKDUP.out)
    GATK_BQSR(INDEX_MARKDUP.out, reference.out.ref_dict)
    INDEX_BQSR(GATK_BQSR.out)
    GATK_MUTECT2_CALL(INDEX_BQSR.out)
    PCGR_ANNOTATE(GATK_MUTECT2_CALL.out)
}
