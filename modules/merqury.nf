#!/usr/bin/env nextflow

/*
Module for assessing assembly completeness and quality using merqury.
*/

process merqury {
    cpus params.merqury_threads
    time params.merqury_time
    queue params.merqury_queue
//    conda params.main_conda

    publishDir "results/${sample}/merqury", mode: 'copy'

    input:
       tuple val(sample), path(fastas), path(meryl_db)

    output:
        path "${sample}_merq*"

    script:
    // Separate FASTAs into primary and haplotypes
    def primary = fastas.find { it.name.endsWith('.p_ctg.fa') && !it.name.contains('hap') }
    def haps    = fastas.findAll { it.name.endsWith('.hap1.p_ctg.fa') || it.name.endsWith('.hap2.p_ctg.fa') }

    // Build strings for command
    def primary_str = primary ? primary.toString() : ''
    def hap_str     = haps ? haps.join(' ') : ''

    """
    # run primary-only merqury (if primary exists)
    ${ primary_str ? "merqury.sh ${meryl_db} ${primary_str} ${sample}_merq_primary" : "" }

    # run diploid merqury (if haplotypes exist)
    ${ hap_str ? "merqury.sh ${meryl_db} ${hap_str} ${sample}_merq_diploid" : "" }
    """
}


// process merqury {
//     cpus params.merqury_threads
//     time params.merqury_time
//     queue params.merqury_queue
//     conda params.merqury_conda_env

//     publishDir "results/${sample}/merqury", mode: 'copy'

//     input:
//        tuple val(sample), path(fastas), path(meryl_db)

//     output:
//         path "${sample}_merq*"

//     script:
//     // `fastas` may be a list (diploid) or a single file (primary)
//     def fasta_str = fastas instanceof List ? fastas.join(' ') : fastas.toString()

//     """
//     merqury.sh ${meryl_db} ${fasta_str} ${sample}_merq
//     """
// }