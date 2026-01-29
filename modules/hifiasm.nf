#!/usr/bin/env nextflow

/*
Module for running hifiasm, either in hifi or nanopore mode
*/

process runHifiasm {
    cpus params.hifiasm_threads
    memory params.hifiasm_memory
    time params.hifiasm_time
    queue params.hifiasm_queue
    //conda params.main_conda

    publishDir "results/${sample}/assemblies", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(mode)

    output:
        tuple val(sample), path("${sample}*.p_ctg.gfa"), path("${sample}*.ec.fq.gz"), emit: files

    script:
    """
    # create the output directory
    # if nanopore, run hifiasm in nanopore mode
    if [ "${mode}" == "ont" ]; then
        echo "Running hifiasm in nanopore mode for sample '${sample}'"
        hifiasm -t ${params.hifiasm_threads} --ont -o ${sample}.asm --dual-scaf --write-ec ${reads}
        # compress error corrected reads
         pigz -p ${params.hifiasm_threads} ${sample}.asm.ec.fq
    elif [ "${mode}" == "hifi" ]; then
        echo "Running hifiasm in hifi mode for sample '${sample}'"
        hifiasm -t ${params.hifiasm_threads} -o ${sample}.asm --dual-scaf ${reads}
        # make empty ec file for hifi mode
        touch ${sample}.asm.ec.fq.gz
    else
        # if mode is not recognized, print an error message
        echo "Error: Unrecognized mode '${mode}' for sample '${sample}'. Must be specified as either hifi (for pacbio hifi reads) or ont (for oxford nanopore)"
        exit 1
    fi
        #hifiasm
    """
}