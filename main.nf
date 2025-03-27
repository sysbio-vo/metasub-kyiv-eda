#!/usr/bin/env nextflow

/*
 * Pipeline params
 */

// Parameters
params.reads_exp = "data/haib*_{1,2}.fastq.gz"
params.krakenDB = "~/krakenDBs/standard8DB"
params.threads = 20
params.kmer_len = 35 
params.read_len = 150
params.classif_lvl = "S"
params.bracken_thresh = 10

// Output directory
params.outdir = "${projectDir}/out"

// Input: Create a channel for paired-end reads
Channel
    .fromFilePairs(params.reads_exp, checkIfExists: true)
    .set { reads_ch }

/*
 * QC with fastp
 */
process FASTP_QC {

    conda "bioconda::fastp=0.24.0"
    publishDir "${params.outdir}/quality-control", mode: "symlink"

    input:
        tuple val(sample_id), path(reads_list)

    output:
        tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: trimmed_reads
        tuple val(sample_id), path("${sample_id}.html"), path("${sample_id}.json"), emit: qc_reports

    script:
    """
    fastp \
        --in1 ${reads_list[0]} --in2 ${reads_list[1]} \
        --out1 ${sample_id}_1.trimmed.fastq.gz \
        --out2 ${sample_id}_2.trimmed.fastq.gz \
        -h ${sample_id}.html -j ${sample_id}.json
    """
}


/*
 * Combine fastp reports with MultiQC
 */
process MULTOQC_JOINTREPORT {
    
    conda "bioconda::multiqc=1.27.1"
    publishDir  "${params.outdir}/quality-control/", mode: "symlink"

    input:
        path(qc_reports)
    output:
        path "multiqc_report.html"

    script:
    """
    multiqc ${params.outdir}/quality-control/
    """
}

/*
 * Run Kraken2 for taxonomic classification
 */
process KRAKEN2 {
    
    conda "bioconda::kraken2=2.1.2"
    publishDir  "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        tuple val(sample_id), path(read1), path(read2) 

    output:
        tuple val(sample_id), path("${sample_id}.report.kraken"), emit: reports
        tuple val(sample_id), path("${sample_id}.out.kraken"), emit: out
        tuple val(sample_id), path("${sample_id}_1.classified.fastq"), path("${sample_id}_2.classified.fastq"), emit: class_reads
        tuple val(sample_id), path("${sample_id}_1.unclassified.fastq"), path("${sample_id}_2.unclassified.fastq"), emit: unclass_reads

    script:
    """
    kraken2 --db ${params.krakenDB} \
        --gzip-compressed --paired --threads 20 \
        --classified-out ${sample_id}#.classified.fastq \
        --unclassified-out ${sample_id}#.unclassified.fastq \
        --report ${sample_id}.report.kraken \
        --output ${sample_id}.out.kraken \
        ${read1} ${read2}
    """
}

/*
 * Run Bracken to reestimate abundances
 */
process BRACKEN {
    
    conda "bioconda::kraken2=2.1.2 bioconda::bracken=2.8"
    publishDir  "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        tuple val(sample_id), path(kraken_report)

    output:
        tuple val(sample_id), path("${sample_id}.report.bracken"), emit: reports
        tuple val(sample_id), path("${sample_id}.out.bracken"), emit: out
        tuple val(sample_id), path("${sample_id}.out.bracken.filtered"), emit: out_filt

    script:
    """
    bracken -d ${params.krakenDB} -i ${kraken_report} \
        -o ${sample_id}.out.bracken  -w ${sample_id}.report.bracken \
        -r ${params.read_len} -l ${params.classif_lvl} \
        -t ${params.bracken_thresh}
    
    grep -v "uncultured" ${sample_id}.out.bracken  > ${sample_id}.out.bracken.filtered
    """
}

/*
 * Combine Kracken reports
 */
process COMBINE_KRAKEN {
    conda "bioconda::krakentools=1.2"
    publishDir "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        path(kraken_reports)
    
    output:
        path("combined_report.kraken") 

    script:
    """
    combine_kreports.py -r ${kraken_reports.join(' ')} \
        -o combined_report.kraken
    """
}

/*
 * Combine bracken output files
 */
process COMBINE_BRACKEN{
    conda "bioconda::bracken=2.8"
    publishDir  "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        path(bracken_out_filt)
    
    output:
        path("combined_out.bracken"), emit: out
    
    script:
    """
    combine_bracken_outputs.py --files ${bracken_out_filt.join(' ')} \
        --output combined_out.bracken
    """
}

/*
 * MultiQC Kraken report
 */
process KRAKEN_MULTIQC {
    conda "bioconda::multiqc=1.27.1"
    publishDir "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        path(kraken_reports)
    
    output:
        path("multiqc_kraken/") 

    script:
    """
    multiqc ${kraken_reports.join(' ')} \
        -o multiqc_kraken
    """
}

/*
 * MultiQC BRACKEN report
 */
process BRACKEN_MULTIQC {
    conda "bioconda::multiqc=1.27.1"
    publishDir "${params.outdir}/taxonomy-analysis/", mode: "symlink"

    input:
        path(bracken_reports)
    
    output:
        path("multiqc_bracken/") 

    script:
    """
    multiqc ${bracken_reports.join(' ')} \
        -o multiqc_bracken
    """
}

/*
 * Taxonomy table from combined Bracken out
 */
process TAXONKIT {
    conda "bioconda::taxonkit=0.19.0"
    publishDir "${params.outdir}/DA_analysis/", mode: "symlink"

    input:
        path(combined_bracken)
    
    output:
        path("taxtable.tsv")

    script:
    """
    cut -f 2 ${combined_bracken} | tail -n+2 | taxonkit reformat -I 1 -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" > taxtable.tsv
    """

}

workflow{
    
    // input channel for the paired-end samples

    fastp_ch = FASTP_QC(reads_ch)
    MULTOQC_JOINTREPORT(fastp_ch.qc_reports.collect( t -> t[1]))

    kraken_ch = KRAKEN2(fastp_ch.trimmed_reads)
    bracken_ch = BRACKEN(kraken_ch.reports)

    COMBINE_KRAKEN(kraken_ch.reports.collect( t -> t[1]))
    comb_bracken_ch = COMBINE_BRACKEN(bracken_ch.out_filt.collect( t -> t[1]))

    KRAKEN_MULTIQC(kraken_ch.reports.collect( t -> t[1]))
    BRACKEN_MULTIQC(bracken_ch.reports.collect( t -> t[1]))

    TAXONKIT(comb_bracken_ch.out)
}

