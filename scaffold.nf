#!/usr/bin/env nextflow

/*
 * Scaffolding Pipeline
 *
 * Start from existing assemblies and add HiC scaffolding, then run stats and reporting.
 * Use this when you already have finished Hifiasm assemblies and want to scaffold them.
 *
 * Input:
 *   --input_reads: CSV file with sample,reads,mode,hic_r1,hic_r2 (same as main pipeline)
 *   --results_dir: Directory where existing assemblies are located (default: results)
 *
 * Usage:
 *   nextflow run scaffold.nf --input_reads samples.csv -resume
 *
 * The pipeline will look for assemblies in: {results_dir}/{sample}/assemblies/*.fa
 * HiC reads and raw reads (for Merqury) are taken from the CSV file.
 */

// Include workflows
include { SCAFFOLDING } from './workflows/scaffolding.nf'
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

    // Extract HiC info from input
    hic_info_ch = input_ch
        .map { sample, reads, mode, hic_r1, hic_r2 ->
            tuple(sample, hic_r1, hic_r2)
        }
        .unique()

    // Run scaffolding
    SCAFFOLDING(assembly_ch, hic_info_ch)

    // Create Meryl databases if reads are provided
    reads_with_hic = input_ch.filter { sample, reads, mode, hic_r1, hic_r2 -> reads != null }

    if (reads_with_hic) {
        merylCounts(reads_with_hic)
        meryl_db_ch = merylCounts.out.meryl_db
    } else {
        // Try to find existing Meryl databases
        log.info "No reads provided in CSV. Looking for existing Meryl databases..."
        meryl_db_ch = Channel.empty()
    }

    // Run QC and stats
    QC_AND_STATS(SCAFFOLDING.out.assemblies, meryl_db_ch)

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
    Scaffolding pipeline completed!
    ========================================
    Scaffolded assemblies: ${params.outdir ?: 'results'}/*/yahs/
    Assembly reports: ${params.outdir ?: 'results'}/*/assembly_report.html
    """
}
