#!/usr/bin/env nextflow

/*
Module for generating comprehensive assembly quality report per sample
*/

process generateAssemblyReport {
    conda params.assembly_report_conda_env
    cpus 1
    time '1h'
    queue params.assembly_report_queue

    publishDir "results/${sample}/", mode: 'copy'

    input:
        tuple val(sample), path(stats_files), path(busco_files), path(merqury_png_files), path(merqury_qv_files), path(repeat_files)

    output:
        tuple val(sample), path("${sample}_assembly_quality_report.html"), emit: report

    script:
    def stats_list = stats_files instanceof List ? stats_files.join(' ') : stats_files
    def busco_list = busco_files instanceof List ? busco_files.join(' ') : busco_files
    def merqury_png_list = merqury_png_files instanceof List ? merqury_png_files.join(' ') : merqury_png_files
    def merqury_qv_list = merqury_qv_files instanceof List ? merqury_qv_files.join(' ') : merqury_qv_files
    def repeat_list = repeat_files instanceof List ? repeat_files.join(' ') : repeat_files
    """
    # Prevent Python from using user site-packages
    export PYTHONNOUSERSITE=1

    generate_assembly_report.py \
        --sample ${sample} \
        --stats ${stats_files.join(' ')} \
        --busco ${busco_files.join(' ')} \
        --merqury-png ${merqury_png_files.join(' ')} \
        --merqury-qv ${merqury_qv_files.join(' ')} \
        --repeats ${repeat_files.join(' ')} \
        --output ${sample}_assembly_quality_report.html
    """
}
