#!/usr/bin/env nextflow

/*
Module for YaHS (Yet Another Hi-C Scaffolding) scaffolding
Takes HiC-mapped BAM file and assembly to produce scaffolded genome
*/

process yahs_scaffold {
    tag "${sample}_${assembly.baseName}"

    cpus params.yahs_threads
    time params.yahs_time
    queue params.yahs_queue
    memory params.yahs_memory

    publishDir "results/${sample}/scaffolded_assemblies", mode: 'copy'

    input:
        tuple val(sample), path(assembly), path(hic_bam)

    output:
        tuple val(sample), path("${assembly.baseName}_scaffolds_final.fa"), emit: scaffolded_assembly
        tuple val(sample), path("${assembly.baseName}_scaffolds_final.agp"), emit: agp_file
        path("${assembly.baseName}_yahs.out"), emit: yahs_log
        path("${assembly.baseName}_scaffolds_final.bin"), optional: true, emit: bin_file

    script:
    """
    # Index the BAM file if not already indexed
    samtools index ${hic_bam}

    # Index the assembly FASTA
    samtools faidx ${assembly}

    # Run YaHS scaffolding
    yahs \\
        -o ${assembly.baseName}_scaffolds \\
        ${assembly} \\
        ${hic_bam} \\
        > ${assembly.baseName}_yahs.out

    # Rename output files for clarity
    mv ${assembly.baseName}_scaffolds_scaffolds_final.fa ${assembly.baseName}_scaffolds_final.fa
    mv ${assembly.baseName}_scaffolds_scaffolds_final.agp ${assembly.baseName}_scaffolds_final.agp

    # Optional: move bin file if it exists
    if [ -f ${assembly.baseName}_scaffolds_scaffolds_final.bin ]; then
        mv ${assembly.baseName}_scaffolds_scaffolds_final.bin ${assembly.baseName}_scaffolds_final.bin
    fi
    """
}
