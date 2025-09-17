# Nextflow pipeline for genome assembly using Hifiasm

## !! Under development !!

Pipeline to assemble pacbio or ONT reads using hifiasm (only tested with hifiasm as of yet). The following steps are currently included in the pipeline. In addition to assemblying the reads with hifiasm, the pipeline also converts the output gfa files to fasta, and performes a number of qc/completeness steps using merqury, busco and asm_stats.py.

## Dependencies

- The pipeline is written in Nextflow (https://www.nextflow.io/) and requires a working installation of Nextflow.
- The custom scripts "gfa2fa.py" and "asm_stats.py" require Python (developed and tested with v.3.12.1)
- [hifiasm >0.25](https://github.com/chhylp123/hifiasm).
- [merqury](https://github.com/marbl/merqury)
- [meryl](https://github.com/marbl/meryl) (dependency of merqury)
- [busco](https://busco.ezlab.org/)

## Getting started

Clone the repository:

```bash
git clone https://github.com/axeljen/hifiasm_nf.git
cd hifiasm_nf
```

### Setting up the environment

If the dependencies are not installed on your system, the environment can be set up using conda.

Make sure you have [conda](https://docs.conda.io/en/latest/) installed.  

Enable conda in the workflow by setting `conda.enabled = true` in the nextflow.config file.

```bash

# create environment from yml file
conda env create -f assembly_env.yaml

# find the path to the main conda environment
conda env list | awk ' $1 == "assembly_env" {print $2}'
# copy the path and paste it in the nextflow.config file as the value for process.conda

```

## Running the pipeline

1) Create a csv file that lists all the samples to be assembled, with three comma-separated columns: sample_id,path/to/reads.fq,type(either "hifi" or "ont"). Any lines starting with "#" will be ignored. See 'input_reads.csv' for an example.

Edit the parameters in the nextflow.config file to fit your environment and needs.

The pipeline has been developed for a slurm cluster, to run on slurm:

```bash

# edit slurm headers in slurm_wrapper.sh to fit your system, then submit the job with: 
sbatch slurm_wrapper.sh

```
