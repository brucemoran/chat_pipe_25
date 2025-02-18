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
// Container images
def containers = [
    fastp   : 'quay.io/biocontainers/fastp:0.23.2--h0b8a92a_1',
    bwa     : 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7',
    samtools: 'quay.io/biocontainers/samtools:1.11--h6270b1f_0',
    gatk    : 'broadinstitute/gatk:4.2.6.1',
    pcgr    : 'sigven/pcgr:latest'
]

/////////////////////////////////////////////////////
// Processes
/////////////////////////////////////////////////////

//workflow for reference
workflow reference {
    take:
    genome_url_ch
    genome_base_ch

    main:
    DOWNLOAD_REFERENCE(genome_url_ch, genome_base_ch)
    
    emit:
    ref_fasta = DOWNLOAD_REFERENCE.out[0]
    ref_fai = DOWNLOAD_REFERENCE.out[1]
    ref_amb = DOWNLOAD_REFERENCE.out[2]
    ref_ann = DOWNLOAD_REFERENCE.out[3]
    ref_bwt = DOWNLOAD_REFERENCE.out[4]
    ref_pac = DOWNLOAD_REFERENCE.out[5]
    ref_sa = DOWNLOAD_REFERENCE.out[6]
    ref_dict = DOWNLOAD_REFERENCE.out[7]
}

// Download the reference genome if it does not exist, then index it
process DOWNLOAD_REFERENCE {
    tag 'download_grch38'
    container 'curlimages/curl:latest'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val genome_url
    val genome_base

    output:
    path "${genome_base}.fasta", emit: ref_fa
    path "${genome_base}.fasta.fai", emit: ref_fai
    path "${genome_base}.fasta.amb", emit: ref_amb
    path "${genome_base}.fasta.ann", emit: ref_ann
    path "${genome_base}.fasta.bwt", emit: ref_bwt
    path "${genome_base}.fasta.pac", emit: ref_pac
    path "${genome_base}.fasta.sa", emit: ref_sa
    path "${genome_base}.dict", emit: ref_dict

    script:
    """
    ## get fasta
    echo "Downloading: ${genome_base}.fasta"
    curl -L -o ${genome_base}.fasta ${genome_url}/${genome_base}.fasta

    ##get rest from API also
    ##amb, ann, bwt, pac, sa, fai
    for x in "amb" "ann" "bwt" "fai" "pac" "sa"; do
        echo "Downloading: ${genome_base}.fasta.\$x"
        curl -L -o ${genome_base}.fasta.\$x ${genome_url}/${genome_base}.fasta.\$x
    done
    echo "Downloading: ${genome_base}.dict"
    curl -L -o ${genome_base}.dict ${genome_url}/${genome_base}.dict
    """
}

// Run fastp QC on single- or paired-end FASTQs
process FASTP_QC {
    tag { sample_id }
    container containers.fastp
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("cleaned_*.fastq.gz")

    script:
    // If fq2 is empty, run single-end; otherwise, run paired-end
    def cmd
    if (fq2) {
        cmd = """
        fastp \\
            --in1 ${fq1} --in2 ${fq2} \\
            --out1 cleaned_R1.fastq.gz --out2 cleaned_R2.fastq.gz \\
            --thread 4
        """
    } else {
        cmd = """
        fastp \\
            --in1 ${fq1} \\
            --out1 cleaned_R1.fastq.gz \\
            --thread 4
        """
    }
    """
    $cmd
    """
}

// Align reads with bwa mem and produce a sorted BAM
process ALIGN_BWA {
    tag { sample_id }
    container containers.bwa
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_reads)
    path ref_fa   
    path ref_fai
    path ref_dict
    path ref_ann
    path ref_amb
    path ref_bwt
    path ref_pac
    path ref_sa

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")
    path ref_fa
    path ref_fai
    path ref_dict

    script:
    // cleaned_reads might have 1 or 2 fastq files
    def r1 = cleaned_reads.find { it.name.contains('R1') }
    def r2 = cleaned_reads.find { it.name.contains('R2') }
    """
    # Ensure bwa index is present
    bwa index ${ref_fasta} || true

    if [ -f "${r2}" ]; then
        bwa mem -t 4 ${ref_fasta} ${r1} ${r2} | samtools sort -o ${sample_id}.sorted.bam
    else
        bwa mem -t 4 ${ref_fasta} ${r1}       | samtools sort -o ${sample_id}.sorted.bam
    fi

    samtools index ${sample_id}.sorted.bam
    """
}

// Mark duplicates & perform Base Quality Score Recalibration
process MARKDUP_BQSR {
    tag { sample_id }
    container containers.gatk
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
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.dedup.bam \\
        -M ${sample_id}.dedup.metrics.txt \\
        --VALIDATION_STRINGENCY SILENT \\
        --REMOVE_SEQUENCING_DUPLICATES false

    samtools index ${sample_id}.dedup.bam

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

    samtools index ${outBamPrefix}.bam
    """
}

// Call somatic variants with GATK Mutect2
process MUTECT2_CALL {
    tag { sample_id }
    container containers.gatk
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
    container containers.pcgr
    publishDir "${params.outdir}/pcgr", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    // All PCGR outputs typically go into a directory named <sample_id>
    // plus an annotated VCF. Adjust as needed.
    path "${sample_id}.*"

    script:
    def outPrefix = "${sample_id}.annotated"
    """
    pcgr --input_vcf ${vcf} \\
         --genome_assembly grch38 \\
         --pcgr_dir /pcgr_data \\
         --sample_id ${sample_id} \\
         --output_dir . \\
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
    sample_info_ch.view()
    //FASTP_QC(sample_info_ch)
    // ALIGN_BWA(sample_info_ch,
    //           reference.out.ref_fa, 
    //           reference.out.ref_fai,
    //           reference.out.ref_dict,
    //           reference.out.ref_ann,
    //           reference.out.ref_amb,
    //           reference.out.ref_bwt,
    //           reference.out.ref_pac,
    //           reference.out.ref_sa)
    // MARKDUP_BQSR(ALIGN_BWA.out)
    // MUTECT2_CALL(MARKDUP_BQSR.out)
    // PCGR_ANNOTATE(MUTECT2_CALL.out)
}
