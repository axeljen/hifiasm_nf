#!/usr/bin/env nextflow

/*
Module for assessing assembly completeness using BUSCO.
*/

process busco {
    cpus params.busco_threads
    time params.busco_time
    queue params.busco_queue
//    conda params.busco_conda

    publishDir "results/$sample/busco", mode: 'copy'

    input:
       // using the reads channel here, we can ignore the mode though
       tuple val(sample), path(fasta)

    output:
        path "${fasta.baseName}.busco", emit: busco_summary

    script:

    """
    busco -i ${fasta} -o ${fasta.baseName}.busco -l ${params.busco_lineage} -m genome -c ${params.busco_threads}
    """
}