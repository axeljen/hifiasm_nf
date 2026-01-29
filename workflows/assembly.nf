#!/usr/bin/env nextflow

// Assembly workflow - handles Hifiasm assembly and GFA to FASTA conversion

include { runHifiasm } from '../modules/hifiasm.nf'
include { gfa2fa } from '../modules/gfa2fa.nf'

workflow ASSEMBLY {
    take:
    inreads_ch  // Channel: [sample, reads, mode, hic_r1, hic_r2]

    main:
    // Run hifiasm to generate assemblies
    assemblies_ch = runHifiasm(inreads_ch)

    // For ONT samples, use ec reads; for HiFi samples, use original reads
    reads_ch = inreads_ch
        .map { sample, reads, mode, hic_r1, hic_r2 ->
            tuple(sample, reads, mode)
        }
        .join(assemblies_ch.files.map { sample, gfas, ec_reads -> tuple(sample, ec_reads) })
        .map { sample, orig_reads, mode, ec_reads ->
            def reads_to_use = (mode == 'ont') ? ec_reads : orig_reads
            tuple(sample, reads_to_use)
        }

    // Flatten GFA files for conversion
    gfas_ch = assemblies_ch.files
        .flatMap { sample, gfas, ec_reads ->
            gfas.collect { f -> tuple(sample, f) }
        }

    // Convert GFA to FASTA
    fasta_ch = gfa2fa(gfas_ch).fasta

    emit:
    fasta = fasta_ch       // Channel: [sample, fasta_file]
    reads = reads_ch       // Channel: [sample, reads] - ec for ONT, original for HiFi
}
