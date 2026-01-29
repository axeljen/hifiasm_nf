#!/usr/bin/env nextflow

/*
Module for assessing assembly completeness using BUSCO.
*/

process busco {
    cpus params.busco_threads
    time params.busco_time
    queue params.busco_queue
//    conda params.busco_conda

    publishDir "results/$sample/busco", mode: 'copy', pattern: "*.{txt,json,tsv}"

    input:
       // using the reads channel here, we can ignore the mode though
       tuple val(sample), path(fasta)
       path(busco_db)

    output:
        tuple val(sample), path("*.txt"), emit: busco_summary
        tuple val(sample), path("*.json"), emit: busco_json, optional: true
        tuple val(sample), path("*.tsv"), emit: busco_table, optional: true

    script:

    """
    # Set up scratch directory for BUSCO
    BUSCO_TMPDIR=${params.scratch_dir}/${fasta.baseName}_busco_\${SLURM_JOB_ID:-\$\$}
    mkdir -p \${BUSCO_TMPDIR}

    # Run BUSCO in scratch storage
    busco -i ${fasta} \\
          -o ${fasta.baseName}.busco \\
          -l ${params.busco_lineage} \\
          -m genome \\
          -c ${params.busco_threads} \\
          --offline \\
          --download_path ${busco_db} \\
          --out_path \${BUSCO_TMPDIR}

    # Copy only the important summary files to work directory
    cp \${BUSCO_TMPDIR}/${fasta.baseName}.busco/short_summary.*.txt . || true
    cp \${BUSCO_TMPDIR}/${fasta.baseName}.busco/short_summary.json . || true

    # Copy the full table (small text file with detailed results)
    cp \${BUSCO_TMPDIR}/${fasta.baseName}.busco/run_*/full_table.tsv ${fasta.baseName}.busco_full_table.tsv || true

    # Clean up scratch storage
    # rm -rf \${BUSCO_TMPDIR}
    """
}

process busco_download_db {
    cpus 1
    time '1h'
    queue params.busco_download_queue

    publishDir "results/busco_db/lineages/${params.busco_lineage}/", mode: 'copy'

    output:
        path "busco_db/", emit: busco_db
    script:
    """
    # Download BUSCO lineage dataset
    busco --download ${params.busco_lineage} --download_path busco_db
    """
}