#!/bin/sh   
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mail-user=XXXXXX
#SBATCH --job-name=snakemake
#SBATCH -o snakemake_stdout_stderr.txt
#SBATCH --mail-type=END

source fastqc-0.11.9 #fastqc
source package 65df873c-d601-44ba-ac61-64644b55dfbb #quast
source package /tgac/software/testing/bin/centrifuge-1.0.4_6cc874e #centrifuge

snakemake -j 225 -p --snakefile GenerateSamplesheet.smk --cluster-config cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time}"
 
snakemake -j 225 -p --snakefile SkimSeqQC.smk --cluster-config cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time}"
