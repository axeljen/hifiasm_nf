# Nextflow pipeline for genome assembly using Hifiasm

## !! Under development !!

Pipeline to assemble pacbio or ONT reads using hifiasm. The following steps are currently included in the pipeline:

1. Assemble reads using hifiasm. 
2. Convert relevant GFA output to FASTA format using the custom gfa2fa.py script.
3. Run the custom 'asm_stats.py' script to generate some standard assembly statistics (N50/L50, numbed of contigs, etc.).
4. Run Compleasm to assess assembly quality.

All results will be collected inside the ./results directory.

## Dependencies

- The pipeline is written in Nextflow (https://www.nextflow.io/) and requires a working installation of Nextflow.
- The custom scripts 'gfa2fa.py' and 'asm_stats.py' require Python (developed and tested with v.3.12.1)
- Requires hifiasm >0.25 and compleasm https://github.com/huangnengCSU/compleasm.
- The recommended way is to run the pipeline with Singularity/Apptainer. After installing apptainer/singularity, pull the container with the dependencies preinstalled:
```
apptainer pull oras://community.wave.seqera.io/library/compleasm_hifiasm:953a436935c830be
```



## Running the pipeline

To run the pipeline, clone this directory:

```
git clone https://github.com/axeljen/hifiasm_nf.git

```

Change into the directory

```
cd hifiasm_nf
```

If running the pipeline with the singularity/apptainer container, pull the container to this directory:
```
apptainer pull oras://community.wave.seqera.io/library/compleasm_hifiasm:953a436935c830be
```

Create a csv file that lists all the samples to be assembled, with three comma-separated columns: 1) sample_id, 2) path/to/reads.fq, 3) type (either "hifi" or "ont"). See 'input_reads.csv' for an example. The pipeline has only been tested with hifi reads so far.

Edit the parameters in the nextflow.config file to fit your environment and needs.

The pipeline has been developed for a slurm cluster, to run on slurm:

```

sbatch slurm_wrapper.sh

```


