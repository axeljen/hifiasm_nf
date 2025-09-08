#!/usr/bin/env nextflow

/*
Module for converting the output from hifiasm (gfa) to fasta format
*/

process gfa2fa {
    cpus '1'
    time '1h'
    queue params.gfa2fa_queue

    publishDir 'results/assemblies', mode: 'copy'

    input:
       tuple val(sample), path(gfa)

    output:
        tuple val(sample), path("${sample}/${gfa.baseName}.fa"), emit: fasta

    script:
    """
    mkdir -p ${sample}
    gfa2fa.py ${gfa} ${sample}/${gfa.baseName}.fa
    """
}