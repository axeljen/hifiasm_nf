#!/usr/bin/env nextflow

// Reporting workflow - generates comprehensive assembly quality reports

include { generateAssemblyReport } from '../modules/assembly_report.nf'

workflow ASSEMBLY_REPORT {
    take:
    stats_ch       // Channel: [sample, stats_file]
    busco_ch       // Channel: [sample, busco_file]
    merqury_ch     // Channel: [files]
    repeats_ch     // Channel: [sample, tbl_file]

    main:
    // Group stats files by sample
    stats_by_sample = stats_ch
        .map { sample, stats -> tuple(sample, stats) }
        .groupTuple(by: 0)

    // Group BUSCO files by sample
    busco_by_sample = busco_ch
        .map { sample, busco -> tuple(sample, busco) }
        .groupTuple(by: 0)

    // Group merqury files by sample and split PNG/QV files
    merqury_split = merqury_ch
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

    // Group RepeatMasker files by sample
    repeat_by_sample = repeats_ch
        .map { sample, tbl -> tuple(sample, tbl) }
        .groupTuple(by: 0)

    // Combine all channels by sample
    // [sample, [stats_files], [busco_files], [merqury_png_files], [merqury_qv_files], [repeat_files]]
    report_input = stats_by_sample
        .join(busco_by_sample, by: 0)
        .join(repeat_by_sample, by: 0)
        .join(merqury_by_sample, by: 0)
        .map { items ->
            def sample = items[0]
            def stats = items[1]
            def busco = items[2]
            def repeats = items[3]
            def png_files = items[4]
            def qv_files = items[5]
            tuple(sample, stats, busco, png_files, qv_files, repeats)
        }

    // Generate the report for each sample
    generateAssemblyReport(report_input)

    emit:
    reports = generateAssemblyReport.out
}
