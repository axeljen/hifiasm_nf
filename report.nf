#!/usr/bin/env nextflow

/*
 * Report Generation Pipeline
 *
 * Generate HTML reports from existing QC results without re-running any analyses.
 * Use this when all QC results exist and you just want to regenerate/update the reports.
 *
 * This is the fastest option - it only reads existing files and generates the HTML report.
 * Perfect for:
 * - Updating report formatting/styling
 * - Regenerating reports after manual edits to QC files
 * - Creating reports with different settings
 *
 * Input:
 *   --input_reads: CSV file with sample names (only column 1 is used)
 *   --results_dir: Directory where QC results are located (default: results)
 *
 * Usage:
 *   nextflow run report.nf --input_reads samples.csv -resume
 *
 * The pipeline expects to find existing QC results in:
 *   {results_dir}/{sample}/stats/*.stats.txt
 *   {results_dir}/{sample}/busco/short_summary*.txt
 *   {results_dir}/{sample}/merqury/*.png and *.qv
 *   {results_dir}/{sample}/repeatmasker/*.tbl
 */

// Include workflow
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

    // Gather existing stats files
    stats_ch = input_ch
        .flatMap { sample ->
            def stats_dir = file("${results_dir}/${sample}/assembly_stats")
            if (stats_dir.exists()) {
                stats_dir.listFiles()
                    .findAll { it.name.endsWith('.stats.txt') }
                    .collect { f -> tuple(sample, f) }
            } else {
                log.error "Stats directory not found for sample ${sample}: ${stats_dir}"
                []
            }
        }

    // Gather existing BUSCO files
    busco_ch = input_ch
        .flatMap { sample ->
            def busco_dir = file("${results_dir}/${sample}/busco")
            if (busco_dir.exists()) {
                busco_dir.listFiles()
                    .findAll { it.name.startsWith('short_summary') && it.name.endsWith('.txt') }
                    .collect { f -> tuple(sample, f) }
            } else {
                log.error "BUSCO directory not found for sample ${sample}: ${busco_dir}"
                []
            }
        }

    // Gather existing Merqury files
    merqury_ch = input_ch
        .flatMap { sample ->
            def merqury_dir = file("${results_dir}/${sample}/merqury")
            if (merqury_dir.exists()) {
                merqury_dir.listFiles()
                    .findAll { it.name.endsWith('.png') || it.name.endsWith('.qv') }
            } else {
                log.error "Merqury directory not found for sample ${sample}: ${merqury_dir}"
                []
            }
        }
        .collect()

    // Gather existing RepeatMasker files
    repeats_ch = input_ch
        .flatMap { sample ->
            def repeatmasker_dir = file("${results_dir}/${sample}/repeatmasker")
            if (repeatmasker_dir.exists()) {
                repeatmasker_dir.listFiles()
                    .findAll { it.name.endsWith('.repeats.tbl') }
                    .collect { f -> tuple(sample, f) }
            } else {
                log.error "RepeatMasker directory not found for sample ${sample}: ${repeatmasker_dir}"
                []
            }
        }

    // Generate reports
    ASSEMBLY_REPORT(
        stats_ch,
        busco_ch,
        merqury_ch,
        repeats_ch
    )
}

workflow.onComplete {
    log.info """
    ========================================
    Report generation completed!
    ========================================
    Assembly reports: ${params.outdir ?: 'results'}/*/assembly_report.html

    To view the reports, open them in a web browser:
      firefox ${params.outdir ?: 'results'}/*/assembly_report.html

    Or download them to your local machine:
      scp -r dardel:${workflow.launchDir}/${params.outdir ?: 'results'}/*/assembly_report.html ./
    """
}
