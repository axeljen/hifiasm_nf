#!/usr/bin/env nextflow

/*
Module for converting the output from hifiasm (gfa) to fasta format
*/

process gfa2fa {
    cpus '1'
    time '1h'
    queue params.gfa2fa_queue
//    conda params.main_conda
    
    publishDir "results/${sample}/assemblies", mode: 'copy'

    input:
       tuple val(sample), path(gfa)

    output:
        tuple val(sample), path("${gfa.baseName}.fa.gz"), emit: fasta
        path("${gfa.baseName}.fa.gz.fai"), emit: fasta_index

    script:
    """
    gfa2fa.py ${gfa} ${gfa.baseName}.fa
    # zip and index with samtools
    bgzip ${gfa.baseName}.fa
    samtools faidx ${gfa.baseName}.fa.gz
    """
}