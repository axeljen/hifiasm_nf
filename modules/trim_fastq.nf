#!/usr/bin/env nextflow

/*
Module for trimming HiC fastq reads using fastp
*/

process fastp_trim_hic {
    cpus params.fastp_threads
    time params.fastp_time
    queue params.fastp_queue

    publishDir "results/${sample}/temps/trimmed_hic", mode: 'copy'

    input:
        tuple val(sample), path(reads1), path(reads2)

    output:
        tuple val(sample), path("${sample}_R1.trimmed.fastq.gz"), path("${sample}_R2.trimmed.fastq.gz"), emit: trimmed_reads
        path "${sample}_fastp.html", emit: fastp_report

    script:
    """
    echo "Trimming HiC reads for sample '${sample}' using fastp"
    fastp -i ${reads1} -I ${reads2} -o ${sample}_R1.trimmed.fastq.gz -O ${sample}_R2.trimmed.fastq.gz -w ${params.fastp_threads} -h ${sample}_fastp.html
    """
}