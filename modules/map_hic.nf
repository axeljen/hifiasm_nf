#!/usr/bin/env nextflow

/*
Module for mapping trimmed HiC reads to indexed assembly using BWA
Outputs sorted BAM file, pairs file, and mapping stats for scaffolding
*/

process map_hic {
    tag "${sample}_${assembly.baseName}"

    cpus params.map_hic_threads
    time params.map_hic_time
    queue params.map_hic_queue
    memory params.map_hic_memory

    publishDir "results/${sample}/hic_mapping", mode: 'copy'

    input:
        tuple val(sample), path(assembly), path(index_files), path(reads1), path(reads2)

    output:
        tuple val(sample), path(assembly), path("${assembly.baseName}.bam"), emit: hic_bam
        tuple val(sample), path("${assembly.baseName}.pairs"), emit: hic_pairs
        path("${assembly.baseName}.stats"), emit: mapping_stats

    script:
    """
    # Create fasta index for genome file creation
    samtools faidx ${assembly}

    # Create genome file from the fasta index (chromosome names and sizes)
    cut -f 1,2 ${assembly}.fai > ${assembly.baseName}.genome

    # Map HiC reads with BWA and process with pairtools pipeline
    bwa mem -5SP -T0 -t ${params.map_hic_threads} ${assembly} ${reads1} ${reads2} | \\
        pairtools parse --min-mapq 40 --walks-policy 5unique --add-columns mapq \\
        --max-inter-align-gap 30 --nproc-in ${params.map_hic_threads} --nproc-out ${params.map_hic_threads} --chroms-path ${assembly.baseName}.genome | \\
        pairtools sort --tmpdir=\${TMPDIR:-.} --nproc ${params.map_hic_threads} | \\
        pairtools dedup --nproc-in ${params.map_hic_threads} \\
        --nproc-out ${params.map_hic_threads} --mark-dups --output-stats ${assembly.baseName}.stats | \\
        pairtools split --nproc-in ${params.map_hic_threads} \\
        --nproc-out ${params.map_hic_threads} --output-pairs ${assembly.baseName}.pairs --output-sam - | \\
        samtools view -bS -@ ${params.map_hic_threads} | \\
        samtools sort -@ ${params.map_hic_threads} -o ${assembly.baseName}.bam
    """
}