// Module for generating HiC heatmaps for Juicebox visualization

/*
Convert pairs file to .hic format for Juicebox
Uses juicer_tools to generate contact maps at multiple resolutions
*/
process pairs_to_hic {
    tag "${sample}_${assembly.baseName}"
    publishDir "results/${sample}/juicer", mode: 'copy'

    cpus params.juicer_cpus
    memory params.juicer_memory
    time params.juicer_time
    queue params.juicer_queue

    input:
        tuple val(sample), path(assembly), path(pairs_file)

    output:
        tuple val(sample), path("${assembly.baseName}.hic"), emit: hic_file
        tuple val(sample), path("${assembly.baseName}.hic"), path(assembly), emit: hic_with_assembly

    script:
    def mem_gb = task.memory.toGiga() - 2
    """
    # Create chromosome sizes file from assembly
    if [[ ${assembly} == *.gz ]]; then
        gunzip -c ${assembly} > assembly_temp.fa
        samtools faidx assembly_temp.fa
        cut -f1,2 assembly_temp.fa.fai > chrom.sizes
        rm assembly_temp.fa assembly_temp.fa.fai
    else
        samtools faidx ${assembly}
        cut -f1,2 ${assembly}.fai > chrom.sizes
    fi

    # Sort pairs file if not already sorted
    # Use LC_ALL=C for faster sorting
    LC_ALL=C sort --parallel=${task.cpus} -S ${mem_gb}G -k2,2 -k4,4 -k3,3n -k5,5n ${pairs_file} > sorted.pairs

    # Convert sorted pairs to .hic format using juicer_tools from conda
    # Resolutions: 2.5M, 1M, 500K, 250K, 100K, 50K, 25K, 10K, 5K
    juicer_tools pre \\
        -j ${task.cpus} \\
        -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 \\
        sorted.pairs \\
        ${assembly.baseName}.hic \\
        chrom.sizes

    # Clean up temp files
    rm sorted.pairs
    """
}

/*
Generate contact matrix in text format for custom visualization
Useful for creating custom plots or additional analysis
*/
process pairs_to_matrix {
    tag "${sample}_${assembly.baseName}"
    publishDir "results/${sample}/contact_matrices", mode: 'copy'

    cpus params.matrix_cpus
    memory params.matrix_memory
    time params.matrix_time
    queue params.matrix_queue

    input:
        tuple val(sample), path(assembly), path(pairs_file)

    output:
        tuple val(sample), path("${assembly.baseName}_matrix.txt"), emit: matrix

    script:
    """
    # Use cooler or custom script to generate contact matrix
    # This is a placeholder - you can customize based on your needs

    # Create chromosome sizes
    if [[ ${assembly} == *.gz ]]; then
        gunzip -c ${assembly} | samtools faidx - > temp.fai
        cut -f1,2 temp.fai > chrom.sizes
    else
        samtools faidx ${assembly}
        cut -f1,2 ${assembly}.fai > chrom.sizes
    fi

    # Convert pairs to bedpe format and bin
    awk 'BEGIN{OFS="\\t"} !/^#/ {print \$2,\$3,\$3+1,\$4,\$5,\$5+1}' ${pairs_file} | \\
        bedtools intersect -a stdin -b stdin -wa -wb | \\
        awk '{bin1=int(\$2/1000000); bin2=int(\$5/1000000); print \$1"_"bin1"\\t"\$4"_"bin2}' | \\
        sort | uniq -c | \\
        awk '{print \$2"\\t"\$3"\\t"\$1}' > ${assembly.baseName}_matrix.txt
    """
}

/*
Generate assembly AGP file for Juicebox Assembly Tools (optional)
Useful if you want to manually curate scaffolds in Juicebox
*/
process generate_assembly_agp {
    tag "${sample}_${assembly.baseName}"
    publishDir "results/${sample}/juicer", mode: 'copy'

    cpus 1
    memory '4GB'
    time '1h'
    queue params.juicer_queue

    input:
        tuple val(sample), path(assembly)

    output:
        tuple val(sample), path("${assembly.baseName}.assembly"), emit: assembly_file

    script:
    """
    # Create assembly file for Juicebox
    if [[ ${assembly} == *.gz ]]; then
        gunzip -c ${assembly} > temp.fa
        FASTA=temp.fa
    else
        FASTA=${assembly}
    fi

    # Generate simple assembly format (scaffold_name length)
    awk '/^>/ {if(name) print name"\\t"len; name=substr(\$1,2); len=0; next} {len+=length(\$0)} END {print name"\\t"len}' \$FASTA > ${assembly.baseName}.assembly

    [[ ${assembly} == *.gz ]] && rm temp.fa || true
    """
}
