#!/usr/bin/env nextflow

/*
Module for generating meryl k-mer count data base from input reads. Later used by merqury to assess assembly quality.
*/

process merylCounts {
    cpus params.merylCounts_threads
    time params.merylCounts_time
    queue params.merylCounts_queue
//    conda params.main_conda

    // publishDir "results/$sample/meryl_db", mode: 'copy'

    input:
       // using the reads channel here, we can ignore the mode though
       tuple val(sample), path(fasta), val(mode)

    output:
        tuple val(sample), path("${fasta.baseName}.meryl"), emit: meryl_db

    script:

    """
    meryl k=${params.meryl_kmer_size} count output ${fasta.baseName}.meryl ${fasta}
    """
}