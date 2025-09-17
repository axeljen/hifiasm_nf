#!/usr/bin/env nextflow

/*
Module for running hifiasm, either in hifi or nanopore mode
*/

process runHifiasm {
    cpus params.hifiasm_threads
    time params.hifiasm_time
    queue params.hifiasm_queue
    //conda params.main_conda

    publishDir "results/${sample}/assemblies", mode: 'copy'

    input:
        tuple val(sample), path(reads), val(mode)

    output:
        tuple val(sample), path("${sample}*.p_ctg.gfa"), emit: files

    script:
    """
    # create the output directory
    # if nanopore, run hifiasm in nanopore mode
    if [ "${mode}" == "ont" ]; then
        echo "Running hifiasm in nanopore mode for sample '${sample}'"
        hifiasm -t ${params.hifiasm_threads} --ont -o ${sample}.asm ${reads}
    elif [ "${mode}" == "hifi" ]; then
        echo "Running hifiasm in hifi mode for sample '${sample}'"
        hifiasm -t ${params.hifiasm_threads} -o ${sample}.asm ${reads}
    else
        # if mode is not recognized, print an error message
        echo "Error: Unrecognized mode '${mode}' for sample '${sample}'. Must be specified as either hifi (for pacbio hifi reads) or ont (for oxford nanopore)"
        exit 1
    fi
        #hifiasm
    """
}