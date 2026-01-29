#!/usr/bin/env nextflow

// QC and statistics workflow - runs assembly QC tools

include { asmStats } from '../modules/asm_stats.nf'
include { merqury } from '../modules/merqury.nf'
include { busco; busco_download_db } from '../modules/busco.nf'
include { repeatMasker; create_replib } from '../modules/repeatmasker.nf'
include { merylCounts } from '../modules/meryl.nf'

workflow QC_AND_STATS {
    take:
    assemblies_ch   // Channel: [sample, fasta_file]
    reads_ch        // Channel: [sample, reads_file]

    main:
    // Run assembly stats in parallel with meryl
    asmStats(assemblies_ch)
    
    // Run meryl database creation in parallel
    meryl_db_ch = merylCounts(reads_ch).meryl_db

    // Run Merqury once per sample with all assemblies grouped
    fasta_grouped_ch = assemblies_ch
        .groupTuple(by: 0)  // [sample, [fa1, fa2, fa3]]

    merqury_in_ch = fasta_grouped_ch
        .join(meryl_db_ch, by: 0)  // Join by sample name

    merqury(merqury_in_ch)

    // Run BUSCO completeness assessment
    db_ch = busco_download_db()
    busco(assemblies_ch, db_ch)

    // Run RepeatMasker
    // Configure repeat library
    if (params.repeatMasker_library_path) {
        repeat_lib_ch = Channel.fromPath(params.repeatMasker_library_path)
    } else {
        create_replib(params.repeatMasker_species)
        repeat_lib_ch = create_replib.out
    }

    repeatMasker(assemblies_ch, repeat_lib_ch)

    emit:
    stats = asmStats.out.fasta          // Channel: [sample, stats_file]
    merqury = merqury.out               // Channel: [files]
    busco = busco.out.busco_summary     // Channel: [sample, busco_file]
    repeats = repeatMasker.out.repeats_tbl  // Channel: [sample, tbl_file]
}
