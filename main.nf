#!/usr/bin/env nextflow

/*
 * Full Hifiasm Assembly Pipeline
 *
 * This is the complete pipeline that runs:
 * 1. Hifiasm assembly
 * 2. HiC scaffolding (if HiC reads provided)
 * 3. QC and statistics (assembly stats, Merqury, BUSCO, RepeatMasker)
 * 4. Comprehensive HTML report generation
 *
 * Usage:
 *   nextflow run main.nf --input_reads samples.csv
 */

// Include workflows
include { ASSEMBLY } from './workflows/assembly.nf'
include { SCAFFOLDING } from './workflows/scaffolding.nf'
include { QC_AND_STATS } from './workflows/qc_and_stats.nf'
include { ASSEMBLY_REPORT } from './workflows/reporting.nf'

workflow {
    // Parse input CSV
    // Format: sample,reads,mode,hic_r1,hic_r2
    inreads_ch = Channel.fromPath(params.input_reads)
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

    // Step 1: Run assembly
    ASSEMBLY(inreads_ch)

    // Step 2: Run scaffolding with HiC (if available)
    // Extract HiC info from input
    hic_info_ch = inreads_ch
        .map { sample, reads, mode, hic_r1, hic_r2 ->
            tuple(sample, hic_r1, hic_r2)
        }
        .unique()

    SCAFFOLDING(ASSEMBLY.out.fasta, hic_info_ch)

    // Step 3: Run QC and stats on final assemblies
    QC_AND_STATS(SCAFFOLDING.out.assemblies, ASSEMBLY.out.reads)

    // Step 4: Generate comprehensive reports
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
    Pipeline completed successfully!
    ========================================
    Results directory: ${params.outdir ?: 'results'}
    Assembly reports: ${params.outdir ?: 'results'}/*/assembly_report.html
    """
}
