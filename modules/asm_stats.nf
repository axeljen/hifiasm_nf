#!/usr/bin/env nextflow

/*
Module for calculating basic assembly summary statiscs
*/

process asmStats {
    cpus '1'
    time '1h'
    queue params.asmStats_queue
//    conda params.main_conda

    publishDir "results/${sample}/assembly_stats", mode: 'copy'

    input:
       tuple val(sample), path(fasta)

    output:
        path "${fasta.baseName}.stats.txt", emit: fasta

    script:

    """
    mkdir -p stats
    asm_stats.py ${fasta} > ${fasta.baseName}.stats.txt
    """
}