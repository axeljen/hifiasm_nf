#!/usr/bin/env nextflow

/*
Module for BWA indexing the genome assembly
Used for optional HiC mapping and scaffolding
*/

process bwa_index {
    tag "${sample}"

    cpus params.bwa_threads
    time params.bwa_time
    queue params.bwa_queue

    // Publish only the index files; assembly is already published by hifiasm
    publishDir "results/${sample}/assemblies", mode: 'copy', pattern: "*.{amb,ann,bwt,pac,sa}"

    input:
        tuple val(sample), path(genome_fasta)

    output:
        tuple val(sample), path(genome_fasta), path("${genome_fasta}.{amb,ann,bwt,pac,sa}"), emit: indexed_genome

    script:
    """
    echo "Indexing genome for sample '${sample}' using BWA"
    bwa index -a bwtsw ${genome_fasta}
    """
    }