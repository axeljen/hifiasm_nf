#!/usr/bin/env nextflow

/*
 * RepeatMasker and Report Pipeline
 *
 * Start from existing assemblies, run only RepeatMasker, and regenerate reports.
 * Use this when you want to update repeat annotations without re-running other QC tools.
 *
 * This workflow expects that other QC results (stats, BUSCO, Merqury) already exist
 * in the expected locations, as it will gather them for report generation.
 *
 * Input:
 *   --input_reads: CSV file with sample,reads,mode,hic_r1,hic_r2 (same as main pipeline)
 *   --results_dir: Directory where existing assemblies and QC results are located (default: results)
 *
 * Usage:
 *   nextflow run repeatmasker_report.nf --input_reads samples.csv -resume
 *
 * The pipeline will:
 * - Look for assemblies in: {results_dir}/{sample}/assemblies/*.fa or {results_dir}/{sample}/yahs/*.fa
 * - Look for existing stats in: {results_dir}/{sample}/stats/*.stats.txt
 * - Look for existing BUSCO in: {results_dir}/{sample}/busco/short_summary*.txt
 * - Look for existing Merqury in: {results_dir}/{sample}/merqury/*
 */

// Include modules and workflows
include { repeatMasker; create_replib } from './modules/repeatmasker.nf'
include { ASSEMBLY_REPORT } from './workflows/reporting.nf'

workflow {
    // Parse input CSV (only first column used for sample names)
    // Format: sample,reads,mode,hic_r1,hic_r2
    input_ch = Channel.fromPath(params.input_reads)
        .splitCsv(header: false, sep: ',')
        .filter { row -> !row[0].startsWith('#') }
        .map { row -> row[0] }  // Extract only sample name
        .unique()

    // Define results directory
    results_dir = params.results_dir ?: 'results'

    // Find existing assemblies in results directory
    // Priority: scaffolded_assemblies > assemblies
    assembly_ch = input_ch
        .flatMap { sample ->
            // Check scaffolded_assemblies first (takes precedence)
            def scaffolded_dir = file("${results_dir}/${sample}/scaffolded_assemblies")
            def assembly_dir = file("${results_dir}/${sample}/assemblies")

            def assemblies = []

            // Try scaffolded assemblies first
            if (scaffolded_dir.exists()) {
                assemblies = scaffolded_dir.listFiles()
                    .findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') ||
                              it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }
                    .collect { fasta -> tuple(sample, fasta) }
            }

            // If no scaffolded assemblies, use original assemblies
            if (assemblies.isEmpty() && assembly_dir.exists()) {
                assemblies = assembly_dir.listFiles()
                    .findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') ||
                              it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }
                    .collect { fasta -> tuple(sample, fasta) }
            }

            if (assemblies.isEmpty()) {
                log.warn "No assemblies found for sample ${sample} in ${scaffolded_dir} or ${assembly_dir}"
            }

            assemblies
        }

    // Configure repeat library
    if (params.repeatMasker_library_path) {
        repeat_lib_ch = Channel.fromPath(params.repeatMasker_library_path)
    } else {
        create_replib(params.repeatMasker_species)
        repeat_lib_ch = create_replib.out
    }

    // Run RepeatMasker
    repeatMasker(assembly_ch, repeat_lib_ch)

    // Gather existing QC results
    // Stats files
    stats_ch = input_ch
        .flatMap { sample ->
            def stats_dir = file("${results_dir}/${sample}/assembly_stats")
            if (stats_dir.exists()) {
                stats_dir.listFiles()
                    .findAll { it.name.endsWith('.stats.txt') }
                    .collect { f -> tuple(sample, f) }
            } else {
                log.warn "Stats directory not found for sample ${sample}: ${stats_dir}"
                []
            }
        }

    // BUSCO files
    busco_ch = input_ch
        .flatMap { sample ->
            def busco_dir = file("${results_dir}/${sample}/busco")
            if (busco_dir.exists()) {
                busco_dir.listFiles()
                    .findAll { it.name.startsWith('short_summary') && it.name.endsWith('.txt') }
                    .collect { f -> tuple(sample, f) }
            } else {
                log.warn "BUSCO directory not found for sample ${sample}: ${busco_dir}"
                []
            }
        }

    // Merqury files
    merqury_ch = input_ch
        .flatMap { sample ->
            def merqury_dir = file("${results_dir}/${sample}/merqury")
            if (merqury_dir.exists()) {
                merqury_dir.listFiles()
                    .findAll { it.name.endsWith('.png') || it.name.endsWith('.qv') }
            } else {
                log.warn "Merqury directory not found for sample ${sample}: ${merqury_dir}"
                []
            }
        }
        .collect()

    // Generate updated reports with new RepeatMasker results
    ASSEMBLY_REPORT(
        stats_ch,
        busco_ch,
        merqury_ch,
        repeatMasker.out.repeats_tbl
    )
}

workflow.onComplete {
    log.info """
    ========================================
    RepeatMasker analysis completed!
    ========================================
    Updated reports: ${params.outdir ?: 'results'}/*/assembly_report.html
    RepeatMasker results: ${params.outdir ?: 'results'}/*/repeatmasker/
    """
}
