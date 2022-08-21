#!/bin/sh   
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mail-user=XXXX
#SBATCH --job-name=snakemake
#SBATCH -o snakemake_stdout_stderr.txt
#SBATCH --mail-type=ALL

source fastqc-0.11.9
source package 65df873c-d601-44ba-ac61-64644b55dfbb #quast
source package /tgac/software/testing/bin/centrifuge-1.0.4_6cc874e

snakemake -j 225 --snakefile /hpc-home/mcgowan/workflows/SkimSeqQC/GenerateSamplesheet.smk --cluster-config /hpc-home/mcgowan/workflows/SkimSeqQC/cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time}"
 
snakemake -j 230 --retries 2 --snakefile /hpc-home/mcgowan/workflows/SkimSeqQC/SkimSeqQC.smk --cluster-config /hpc-home/mcgowan/workflows/SkimSeqQC/cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time}"
