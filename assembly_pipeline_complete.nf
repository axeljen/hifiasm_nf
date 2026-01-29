#!/usr/bin/env nextflow

// Include modules
include { runHifiasm } from './modules/hifiasm.nf'
include { gfa2fa } from './modules/gfa2fa.nf'
include { asmStats } from './modules/asm_stats.nf'
include { merylCounts } from './modules/meryl.nf'
include { merqury } from './modules/merqury.nf'
include { busco } from './modules/busco.nf'
include { busco_download_db } from './modules/busco.nf'
include { create_replib } from './modules/repeatmasker.nf'
include { repeatMasker } from './modules/repeatmasker.nf'
include { download_dfam } from './modules/repeatmasker.nf'
include { generateAssemblyReport } from './modules/assembly_report.nf'
include { bwa_index } from './modules/bwa_index.nf'
include { fastp_trim_hic } from './modules/trim_fastq.nf'
include { map_hic } from './modules/map_hic.nf'
include { yahs_scaffold } from './modules/yahs.nf'
include { pairs_to_hic; generate_assembly_agp } from './modules/juicer.nf'

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
            def hic_r1 = row[3] ? file(row[3]) : null
            def hic_r2 = row[4] ? file(row[4]) : null
            return [sample, reads, mode, hic_r1, hic_r2]
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
    fasta_ch = gfa2fa(gfas_ch).fasta

    // Join fasta channel with original input to get HiC info
    // Keep only sample and hic_r1, hic_r2 from inreads_ch
    hic_info_ch = inreads_ch
        .map { sample, reads, mode, hic_r1, hic_r2 ->
            tuple(sample, hic_r1, hic_r2)
        }
        .unique()  // One entry per sample

    // Join fasta with HiC info: [sample, fasta, hic_r1, hic_r2]
    fasta_with_hic_ch = fasta_ch
        .map { sample, fasta -> tuple(sample, fasta) }
        .combine(hic_info_ch, by: 0)

    // Branch into HiC and non-HiC samples
    fasta_with_hic_ch.branch {
        hic: it[2] != null      // Has HiC reads (hic_r1 is not null)
        no_hic: it[2] == null   // No HiC reads
    }.set { branched_ch }

    // For HiC samples: run bwa_index and fastp_trim_hic in parallel
    // Index assemblies for HiC samples (runs once per assembly)
    indexed_assemblies = bwa_index(
        branched_ch.hic.map { sample, fasta, hic_r1, hic_r2 ->
            tuple(sample, fasta)
        }
    )

    // Trim HiC reads once per sample (not once per assembly)
    // Get unique sample+reads combinations
    hic_reads_unique = branched_ch.hic
        .map { sample, fasta, hic_r1, hic_r2 ->
            tuple(sample, hic_r1, hic_r2)
        }
        .unique()

    trimmed_hic = fastp_trim_hic(hic_reads_unique)

    // Combine indexed assemblies with trimmed reads for HiC mapping
    // indexed_assemblies: [sample, assembly, index_files]
    // trimmed_hic: [sample, R1_trimmed, R2_trimmed]
    hic_mapping_input = indexed_assemblies.indexed_genome
        .combine(trimmed_hic.trimmed_reads, by: 0)
        .map { sample, assembly, index_files, reads1, reads2 ->
            tuple(sample, assembly, index_files, reads1, reads2)
        }

    // Map HiC reads to assemblies
    mapped_hic = map_hic(hic_mapping_input)

    // Run YaHS scaffolding
    // mapped_hic.hic_bam already has the format: [sample, assembly, bam]
    scaffolded = yahs_scaffold(mapped_hic.hic_bam)

    // Generate Juicer heatmaps for visualization
    // Combine scaffolded assemblies with pairs files for heatmap generation
    juicer_input = scaffolded.scaffolded_assembly
        .join(mapped_hic.hic_pairs, by: [0,1])  // Join by sample and assembly
        .map { sample, assembly, pairs -> [sample, assembly, pairs] }

    pairs_to_hic(juicer_input)
    generate_assembly_agp(scaffolded.scaffolded_assembly)

    // Use scaffolded assemblies for HiC samples in QC
    hic_assemblies_for_qc = scaffolded.scaffolded_assembly

    // For non-HiC samples: use assemblies as-is
    non_hic_assemblies = branched_ch.no_hic
        .map { sample, fasta, hic_r1, hic_r2 -> tuple(sample, fasta) }

    // Combine both paths for QC
    all_assemblies_for_qc = hic_assemblies_for_qc.mix(non_hic_assemblies)

    // run the assembly stats
    asmStats(all_assemblies_for_qc)

    // run merqury once per sample, a combined run with the diploid fasta files and a single run with the primaries will be performed
    // group fasta files by sample
    fasta_grouped_ch = all_assemblies_for_qc
        .groupTuple(by: 0)  // [ sample, [fa1, fa2, fa3] ]

    // join fasta files with meryl dbs
    merqury_in_ch = fasta_grouped_ch
        .join(merylCounts.out.meryl_db, by: 0)  // join by sample name
    // run merqury
    merqury(merqury_in_ch)

    // run busco completeness assessment
    // start with downloading the busco lineage dataset, as we'll only need to do this once
    db_ch = busco_download_db()
    // and then run busco on all fasta files
    busco(all_assemblies_for_qc, db_ch)

    // if no path to dfam db is provided, download it first
    // dfam_db_path = params.repeatMasker_dfam_path ?
    //     Channel.value(params.repeatMasker_dfam_path) :
    //     Channel.value("${workflow.launchDir}/results/repeatmasker_db/dfam")

    // // download dfam if needed
    // if (!params.repeatMasker_dfam_path) {
    //     download_dfam()
    // }

    // Configure famdb once before running repeatmasker, if it's not already given in the params 
    if (params.repeatMasker_library_path) {
    repeat_lib_ch = Channel.fromPath(params.repeatMasker_library_path)
    } else {
        create_replib(params.repeatMasker_species)
        repeat_lib_ch = create_replib.out
    }

    // run repeatmasker on all fasta files
    // combine fasta channel with famdb_ready to ensure famdb is configured before any repeatMasker runs
    repeatMasker(all_assemblies_for_qc, repeat_lib_ch)

    // Generate comprehensive assembly quality report per sample
    // Group stats files by sample: [sample, stats_file]
    stats_by_sample = asmStats.out.fasta
        .map { sample, stats -> tuple(sample, stats) }
        .groupTuple(by: 0)

    // Group BUSCO summary files by sample: [sample, busco_file]
    busco_by_sample = busco.out.busco_summary
        .map { sample, busco -> tuple(sample, busco) }
        .groupTuple(by: 0)

    // Group merqury files by sample
    // Filter PNG files and QV files separately
    merqury_split = merqury.out
    .flatMap { files -> files.collect { f -> f } }
    .map { file ->
        def sample = file.name.split('_merq_')[0]
        tuple(sample, file)
    }
    .branch {
        png: it[1].name.endsWith('.png')
        qv: it[1].name.endsWith('.qv')
    }

    merqury_png_by_sample = merqury_split.png.groupTuple(by: 0)
    merqury_qv_by_sample = merqury_split.qv.groupTuple(by: 0)

    merqury_by_sample = merqury_png_by_sample
        .join(merqury_qv_by_sample)
        .view { "Merqury by sample: $it" }
        
    // merqury_by_sample = merqury.out
    //     .flatMap { files ->
    //         files.collect { f -> f }
    //     }
    //     .map { file ->
    //         // Extract sample name from filename (e.g., "U_arctos_18_pr_137_001_merq_primary.qv")
    //         def sample = file.name.split('_merq_')[0]
    //         tuple(sample, file)
    //     }
    //     .groupTuple(by: 0)
    //     .map { sample, files ->
    //         def png_files = files.findAll { it.name.endsWith('.png') }
    //         def qv_files = files.findAll { it.name.endsWith('.qv') }
    //         tuple(sample, png_files, qv_files)
    //     }

    // Group RepeatMasker .tbl files by sample: [sample, tbl_file]
    repeat_by_sample = repeatMasker.out.repeats_tbl
        .map { sample, tbl -> tuple(sample, tbl) }
        .groupTuple(by: 0)

    // Combine all channels by sample into a single tuple
    // [sample, [stats_files], [busco_files], [merqury_png_files], [merqury_qv_files], [repeat_files]]
   report_input = stats_by_sample
    .join(busco_by_sample, by: 0)
    .join(repeat_by_sample, by: 0)
    .join(merqury_by_sample, by: 0)
    .view { "After all joins: size=${it.size()}, content=$it" }
    .map { items ->
        def sample = items[0]
        def stats = items[1]
        def busco = items[2]
        def repeats = items[3]
        def png_files = items[4]  // Already a list!
        def qv_files = items[5]   // Already a list!
        tuple(sample, stats, busco, png_files, qv_files, repeats)
    }
    .view { "Report input: $it" }
    .view { sample, stats, busco, png, qv, repeats ->
        """
        Sample: $sample
        Stats count: ${stats instanceof List ? stats.size() : 1}
        Busco count: ${busco instanceof List ? busco.size() : 1}
        PNG count: ${png instanceof List ? png.size() : 1}
        QV count: ${qv instanceof List ? qv.size() : 1}
        Repeats count: ${repeats instanceof List ? repeats.size() : 1}
        PNG files: $png
        QV files: $qv
        """
    }
    //.view { "Report input after map: $it" }
    // report_input = stats_by_sample
    //     .join(busco_by_sample, by: 0)
    //     .join(merqury_by_sample, by: 0)
    //     .view { "Stats by sample: $it" }
    //     .join(repeat_by_sample, by: 0)
    //     .map { sample, stats, busco, merqury_tuple, repeats ->
    //         def (png_files, qv_files) = merqury_tuple
    //         tuple(sample, stats, busco, png_files, qv_files, repeats)
    //     }
    //     .view()

    // Generate the report for each sample
    generateAssemblyReport(report_input)

}