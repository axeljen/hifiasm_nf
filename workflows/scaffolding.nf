#!/usr/bin/env nextflow

// Scaffolding workflow - handles HiC-based scaffolding with YAHS

include { bwa_index } from '../modules/bwa_index.nf'
include { fastp_trim_hic } from '../modules/trim_fastq.nf'
include { map_hic } from '../modules/map_hic.nf'
include { yahs_scaffold } from '../modules/yahs.nf'
include { pairs_to_hic; generate_assembly_agp } from '../modules/juicer.nf'

workflow SCAFFOLDING {
    take:
    fasta_ch        // Channel: [sample, fasta_file]
    hic_info_ch     // Channel: [sample, hic_r1, hic_r2]

    main:
    // Join fasta with HiC info: [sample, fasta, hic_r1, hic_r2]
    fasta_with_hic_ch = fasta_ch
        .combine(hic_info_ch, by: 0)

    // Branch into HiC and non-HiC samples
    fasta_with_hic_ch.branch {
        hic: it[2] != null      // Has HiC reads (hic_r1 is not null)
        no_hic: it[2] == null   // No HiC reads
    }.set { branched_ch }

    // Index assemblies for HiC samples
    indexed_assemblies = bwa_index(
        branched_ch.hic.map { sample, fasta, hic_r1, hic_r2 ->
            tuple(sample, fasta)
        }
    )

    // Trim HiC reads once per sample
    hic_reads_unique = branched_ch.hic
        .map { sample, fasta, hic_r1, hic_r2 ->
            tuple(sample, hic_r1, hic_r2)
        }
        .unique()

    trimmed_hic = fastp_trim_hic(hic_reads_unique)

    // Combine indexed assemblies with trimmed reads for HiC mapping
    hic_mapping_input = indexed_assemblies.indexed_genome
        .combine(trimmed_hic.trimmed_reads, by: 0)
        .map { sample, assembly, index_files, reads1, reads2 ->
            tuple(sample, assembly, index_files, reads1, reads2)
        }

    // Map HiC reads to assemblies
    mapped_hic = map_hic(hic_mapping_input)

    // Run YaHS scaffolding
    scaffolded = yahs_scaffold(mapped_hic.hic_bam)

    // Generate Juicer heatmaps for visualization
    juicer_input = scaffolded.scaffolded_assembly
        .join(mapped_hic.hic_pairs, by: [0,1])
        .map { sample, assembly, pairs -> [sample, assembly, pairs] }

    pairs_to_hic(juicer_input)
    generate_assembly_agp(scaffolded.scaffolded_assembly)

    // Combine scaffolded and non-scaffolded assemblies
    hic_assemblies = scaffolded.scaffolded_assembly
    non_hic_assemblies = branched_ch.no_hic
        .map { sample, fasta, hic_r1, hic_r2 -> tuple(sample, fasta) }

    all_assemblies = hic_assemblies.mix(non_hic_assemblies)

    emit:
    assemblies = all_assemblies  // Channel: [sample, fasta_file]
    hic_pairs = mapped_hic.hic_pairs.ifEmpty([])  // Channel: [sample, assembly, pairs_file]
}
