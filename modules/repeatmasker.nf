#!/usr/bin/env nextflow

/*
Module for running repeatmasker on the assemblies
*/

process create_replib {
    conda params.repeatMasker_conda_env
    cpus 1
    time '30m'
    queue params.repeatMasker_queue

    publishDir "results/repeatmasker_db/famdb/", mode: 'copy'

    input:
        val species

    output:
        path "${species}_famdb.fa", emit: famdb_ready

    script:
    """
    # Prevent Python from using user site-packages (avoid version conflicts)
    export PYTHONNOUSERSITE=1

    # Path to famdb.py script in the repeatmasker environment
    FAMDB_SCRIPT=\${CONDA_PREFIX}/share/RepeatMasker/famdb.py
    FAMDB_DIR=\${CONDA_PREFIX}/share/RepeatMasker/Libraries/famdb

    # Check if famdb.py exists
    if [ ! -f "\${FAMDB_SCRIPT}" ]; then
        echo "Error: famdb.py not found at \${FAMDB_SCRIPT}"
        exit 1
    fi

    # Run famdb.py info to ensure the database is accessible and properly configured
    # This will fail if the h5 files are not properly set up
    python \${FAMDB_SCRIPT} -i \${FAMDB_DIR} families ${species} --curated --ancestors --descendants --include-class-in-name --stage 0 -f fasta_name > ${species}_famdb.fa
    """
}

process repeatMasker {
    conda params.repeatMasker_conda_env
    cpus params.repeatMasker_cpus
    time params.repeatMasker_time
    queue params.repeatMasker_queue

    publishDir "results/${sample}/repeatmasker", mode: 'copy'

    input:
       tuple val(sample), path(fasta)
       path famdb_ready

    output:
        tuple val(sample), path("${fasta.baseName}.repeats.gff"), emit: gff
        tuple val(sample), path("${fasta.baseName}.repeats.out"), emit: repeats_out
        tuple val(sample), path("${fasta.baseName}.repeats.tbl"), emit: repeats_tbl
        tuple val(sample), path("${fasta.baseName}.repeats.masked.fa.gz"), emit: masked_fasta, optional: true

    script:

    """
    # Prevent Python from using user site-packages (avoid version conflicts)
    export PYTHONNOUSERSITE=1

    REPMASKER_TMPDIR=${params.scratch_dir}/${fasta.baseName}_repeatmasker_\${SLURM_JOB_ID:-\$\$}
    mkdir -p \${REPMASKER_TMPDIR}

    # Check if input is gzipped and decompress if necessary
    if [[ ${fasta} == *.gz ]]; then
        echo "Input is gzipped, decompressing..."
        gunzip -c ${fasta} > ${fasta.baseName}
        INPUT_FASTA=${fasta.baseName}
    else
        INPUT_FASTA=${fasta}
    fi

    RepeatMasker -pa ${params.repeatMasker_cpus} \${INPUT_FASTA} -dir \${REPMASKER_TMPDIR} -gff -lib ${famdb_ready}

    # Copy GFF file from scratch to working directory with appropriate naming
    cp \${REPMASKER_TMPDIR}/\${INPUT_FASTA}.out.gff ${fasta.baseName}.repeats.gff
    cp \${REPMASKER_TMPDIR}/\${INPUT_FASTA}.out ${fasta.baseName}.repeats.out
    cp \${REPMASKER_TMPDIR}/\${INPUT_FASTA}.tbl ${fasta.baseName}.repeats.tbl

    # Optionally generate masked fasta
    if [ "${params.repeatMasker_generate_masked_fasta}" = "true" ]; then
        gzip -c \${REPMASKER_TMPDIR}/\${INPUT_FASTA}.masked > ${fasta.baseName}.repeats.masked.fa.gz
    fi

    # Clean up scratch directory
    # rm -rf \${REPMASKER_TMPDIR}

    """
}

process download_dfam {
    cpus params.repeatMasker_download_cpus
    time params.repeatMasker_download_time
    queue params.repeatMasker_download_queue
    publishDir "results/repeatmasker_db/dfam/", mode: 'copy'

    output:
        path "dfam/*", emit: dfam_db
    script:
    """
    mkdir -p dfam
    wget https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.1.h5.gz -O dfam/dfam39_full.1.h5.gz
    gunzip dfam/dfam39_full.1.h5.gz
    """
}