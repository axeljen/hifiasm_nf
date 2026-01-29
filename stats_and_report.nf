#!/usr/bin/env nextflow

/*
 * Stats and Report Pipeline
 *
 * Start from existing assemblies and run QC/stats and generate reports.
 * Use this when you already have finished assemblies and want to regenerate stats and reports.
 *
 * Input:
 *   --input_reads: CSV file with sample,reads,mode,hic_r1,hic_r2 (same as main pipeline)
 *   --results_dir: Directory where existing assemblies are located (default: results)
 *
 * Usage:
 *   nextflow run stats_and_report.nf --input_reads samples.csv -resume
 *
 * The pipeline will look for assemblies in: {results_dir}/{sample}/assemblies/*.fa
 * Raw reads (for Merqury) are taken from column 2 of the CSV file.
 * If no reads are provided in the CSV, Merqury will be skipped.
 */

// Include workflows
include { QC_AND_STATS } from './workflows/qc_and_stats.nf'
include { ASSEMBLY_REPORT } from './workflows/reporting.nf'
include { merylCounts } from './modules/meryl.nf'

workflow {
    // Parse input CSV (same format as main pipeline)
    // Format: sample,reads,mode,hic_r1,hic_r2
    input_ch = Channel.fromPath(params.input_reads)
        .splitCsv(header: false, sep: ',')
        .filter { row -> !row[0].startsWith('#') }
        .map { row ->
            def sample = row[0]
            def reads = row.size() > 1 && row[1] ? file(row[1]) : null
            def mode = row.size() > 2 ? row[2] : 'hifi'
            def hic_r1 = row.size() > 3 && row[3] ? file(row[3]) : null
            def hic_r2 = row.size() > 4 && row[4] ? file(row[4]) : null
            return [sample, reads, mode, hic_r1, hic_r2]
        }

    // Define results directory
    results_dir = params.results_dir ?: 'results'

    // Find existing assemblies in results directory
    // Priority: scaffolded_assemblies > assemblies
    assembly_ch = input_ch
        .map { sample, reads, mode, hic_r1, hic_r2 -> sample }
        .unique()
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

    // Handle Meryl databases
    reads_ch = input_ch.filter { sample, reads, mode, hic_r1, hic_r2 -> reads != null }

    if (reads_ch.count().val > 0) {
        merylCounts(reads_ch)
        meryl_db_ch = merylCounts.out.meryl_db
    } else {
        meryl_db_ch = Channel.empty()
        log.warn "No reads provided in CSV. Merqury analysis will be skipped."
    }

    // Run QC and stats
    QC_AND_STATS(assembly_ch, meryl_db_ch)

    // Generate reports
    ASSEMBLY_REPORT(
        QC_AND_STATS.out.stats,
        QC_AND_STATS.out.busco,
        QC_AND_STATS.out.merqury,
        QC_AND_STATS.out.repeats
    )
}

workflow.onComplete {
    log.info """
    ========================================
    Stats and report generation completed!
    ========================================
    Assembly reports: ${params.outdir ?: 'results'}/*/assembly_report.html
    """
}
