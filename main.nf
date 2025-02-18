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
params.genome_url  = 'https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta'
params.genome_name = 'GRCh38.fa'

// Container images
def containers = [
    fastp   : 'quay.io/biocontainers/fastp:0.23.2--h0b8a92a_1',
    bwa     : 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7',
    samtools: 'quay.io/biocontainers/samtools:1.11--h6270b1f_0',
    gatk    : 'broadinstitute/gatk:4.2.6.1',
    pcgr    : 'sigven/pcgr:latest'
]

/////////////////////////////////////////////////////
// Channel definitions
/////////////////////////////////////////////////////

// Read the CSV and create a channel of [sample_id, fastq1, fastq2]
Channel
    .fromPath(params.samples_csv)
    .splitCsv(header:true)
    .map { row ->
        // For single-end data, fastq_2 might be empty
        def sample_id = row.sample_id
        def fq1       = row.fastq_1
        def fq2       = row.fastq_2
        return [ sample_id, fq1, fq2 ]
    }
    .set { sample_info_ch }


/////////////////////////////////////////////////////
// Processes
/////////////////////////////////////////////////////

//workflow for reference
workflow reference {
    take:
    genome_url_ch
    genome_name_ch

    main:
    DOWNLOAD_REFERENCE(genome_url_ch, genome_name_ch)
    
    emit:
    DOWNLOAD_REFERENCE.out[0] //ref_fasta
    DOWNLOAD_REFERENCE.out[1] //ref_fai
    DOWNLOAD_REFERENCE.out[2] //ref_dict
}

// Download the reference genome if it does not exist, then index it
process DOWNLOAD_REFERENCE {
    tag 'download_grch38'
    container 'curlimages/curl:latest'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val genome_url
    val genome_name

    output:
    path genome_name        , emit: ref_fasta
    path "${genome_name}.fai", emit: ref_fai
    path "${genome_name.replaceAll('.fa$|.fasta$', '')}.dict", emit: ref_dict

    script:
    """
    if [ ! -f ${genome_name} ]; then
        echo "Downloading reference genome..."
        curl -L -o ${genome_name} ${genome_url}
    else
        echo "Reference genome file already exists."
    fi

    samtools faidx ${genome_name}
    gatk CreateSequenceDictionary -R ${genome_name} -O ${genome_name.replaceAll('.fa$|.fasta$', '')}.dict
    """
}

// Run fastp QC on single- or paired-end FASTQs
process FASTP_QC {
    tag { sample_id }
    container containers.fastp
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), val(fq1), val(fq2)

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
    path ref_fasta     from download_reference.out.ref_fasta
    path ref_fai       from download_reference.out.ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

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
    path ref_fasta from download_reference.out.ref_fasta
    path ref_fai   from download_reference.out.ref_fai
    path ref_dict  from download_reference.out.ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.dedup.recal.bam")

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
        -R ${ref_fasta} \\
        -I ${sample_id}.dedup.bam \\
        --known-sites https://storage.googleapis.com/gatk-best-practices/somatic-b37/1000G_phase1.indels.b37.vcf.gz \\
        --known-sites https://storage.googleapis.com/gatk-best-practices/somatic-b37/dbsnp_138.b37.vcf.gz \\
        -O ${sample_id}.recal_data.table

    gatk ApplyBQSR \\
        -R ${ref_fasta} \\
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
    path ref_fasta from download_reference.out.ref_fasta

    output:
    tuple val(sample_id), path("${sample_id}.mutect2.vcf.gz")

    script:
    """
    gatk Mutect2 \\
        -R ${ref_fasta} \\
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
    //download_reference = DOWNLOAD_REFERENCE(genome_url_ch, genome_name_ch)
    // Simple channels for the reference parameters
    genome_url_ch = Channel.value(params.genome_url)

    genome_name_ch = Channel.value(params.genome_name)

    reference(genome_url_ch, genome_name_ch)

    // Now process each sample in turn
    sample_info_ch \
        | FASTP_QC \
        | ALIGN_BWA(reference.out.ref_fasta, reference.out.ref_fai) \
        | MARKDUP_BQSR(reference.out.ref_fasta, 
                       reference.out.ref_fai, 
                       reference.out.ref_dict) \
        | MUTECT2_CALL(reference.out.ref_fasta) \
        | PCGR_ANNOTATE
}
