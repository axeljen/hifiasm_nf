#!/usr/bin/env nextflow

// Include modules
include { runHifiasm } from './modules/hifiasm.nf'
include { gfa2fa } from './modules/gfa2fa.nf'
include { asmStats } from './modules/asm_stats.nf'
include { merylCounts } from './modules/meryl.nf'
include { merqury } from './modules/merqury.nf'
include { busco } from './modules/busco.nf'

workflow {
    // Define input channel
    inreads_ch = Channel.fromPath(params.input_reads)
        // remove any lines starting with a hash
        .splitCsv(header: false, sep: ',')
        .filter { row -> !row[0].startsWith('#') }
        .map { row -> 
            def sample = row[0]
            def reads = file(row[1])
            def mode = row[2]
            return [sample, reads, mode]
        }

    // run hifiasm to generate assemblies
    assemblies_ch = runHifiasm(inreads_ch)

    // prepare a meryl db in parallel for all input reads
    merylCounts(inreads_ch)

    gfas_ch = assemblies_ch.files
    .flatMap { item ->
        def (sample, files) = item
        files.collect { f -> tuple(sample, f) }
    }

    // convert gfa to fasta
    fasta_ch = gfa2fa(gfas_ch)

    // run the assembly stats
    asmStats(fasta_ch)

    // run merqury once per sample, a combined run with the diploid fasta files and a single run with the primaries will be performed
    // group fasta files by sample
    fasta_grouped_ch = fasta_ch
        .groupTuple(by: 0)  // [ sample, [fa1, fa2, fa3] ]

    // join fasta files with meryl dbs
    merqury_in_ch = fasta_grouped_ch
        .join(merylCounts.out.meryl_db, by: 0)  // join by sample name
    // run merqury
    merqury(merqury_in_ch)

    // run busco completeness assessment
    busco(fasta_ch)

}