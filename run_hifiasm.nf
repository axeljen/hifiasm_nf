#!/usr/bin/env nextflow

// Include modules
include { runHifiasm } from './modules/hifiasm.nf'
include { gfa2fa } from './modules/gfa2fa.nf'
include { asmStats } from './modules/asm_stats.nf'
include { runCompleasm; downloadOdb } from './modules/compleasm.nf'

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

   assemblies_ch = runHifiasm(inreads_ch)

    // there are lots of redundant files among those assemblies, 
   // based on the file path endings
   def suffixes = [
            '.bp.p_ctg.gfa',
            '.bp.hap1.p_ctg.gfa',
            '.bp.hap2.p_ctg.gfa'
        ]
    gfas_ch = assemblies_ch.files
        .collect()
        .map { sample, files ->
            def keep = files.findAll { f ->
                def name = f.getFileName().toString()
                suffixes.any { s -> name.endsWith(s) }
            }
            tuple(sample,keep)
        }
        .flatMap { sample, files ->
            files.collect { f -> tuple(sample, f) }
        }

    fasta_ch = gfa2fa(gfas_ch)

    asmStats(fasta_ch)

    // run the compleasm analysis
    runCompleasm(fasta_ch)

}