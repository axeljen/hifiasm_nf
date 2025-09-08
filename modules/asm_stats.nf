#!/usr/bin/env nextflow

/*
Module for calculating basic assembly summary statiscs
*/

process asmStats {
    cpus '1'
    time '1h'
    queue params.asmStats_queue

    publishDir 'results/assemblies/', mode: 'copy'

    input:
       tuple val(sample), path(fasta)

    output:
        path "${sample}/stats/${fasta.baseName}.stats.txt", emit: fasta

    script:

    """
    mkdir -p ${sample}/stats
    asm_stats.py ${fasta} > ${sample}/stats/${fasta.baseName}.stats.txt
    """
}