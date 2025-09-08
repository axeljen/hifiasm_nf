#!/usr/bin/env nextflow

/*
Module for running busco completeness check with compleasm
*/

process downloadOdb {
    cpus '1'
    time '1h'
    queue params.compleasm_queue

    publishDir "${params.compleasm_db_loc}", mode: 'copy'

    script:

    """
    compleasm download --odb ${params.compleasm_odb} --library_path ${params.compleasm_db_loc} ${params.compleasm_lineage}
    """
}

process runCompleasm {
    cpus params.compleasm_threads
    time params.compleasm_time
    queue params.compleasm_queue

    publishDir 'results/assemblies', mode: 'copy'

    input:
       tuple val(sample), path(fasta)

    output:
        path "${sample}/stats/${fasta.baseName}.compleasm_results.txt", emit: busco

    script:
    """
    mkdir -p compleasm_tmp
    mkdir -p ${sample}/stats
    compleasm run -a ${fasta} -o compleasm_tmp \
        -l ${params.compleasm_lineage} \
        --odb ${params.compleasm_odb} \
        --library_path ${params.compleasm_db_loc} \
        --threads ${params.compleasm_threads}

    mv compleasm_tmp/summary.txt ${sample}/stats/${fasta.baseName}.compleasm_results.txt

    """
}