## SkimSeqQC Pipeline
Snakemake pipeline for QC of low-coverage skim-sequencing of G&amp;T seq data

Expects reads to be in a directory called `READS` with filenames in the format: `{sample}_cDNA_R[12].fastq.gz` and `{sample}_gDNA_R[12].fastq.gz`



Modify the config file `config.yaml` to set run_name and choose min number of gDNA and cDNA reads (see example below).

First run the `GenerateSamplesheet.smk` pipeline which just runs fastp on the raw reads to count reads and generate a "samplesheet":



```
snakemake -j 100 --snakefile GenerateSamplesheet.smk --cluster-config cluter.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude}"
```

Then run the full pipeline:

```
snakemake -j 100 --snakefile ../SkimSeqQC.smk --cluster-config ../cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude}"
```


Example sbatch script:

```
#!/bin/sh   
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mail-user=NNNNNN
#SBATCH --job-name=snakemake
#SBATCH -o snakemake_stdout_stderr.txt
#SBATCH --mail-type=END

source snakemake-5.9.1_CBG
source fastqc-0.11.9
source package 65df873c-d601-44ba-ac61-64644b55dfbb #quast
source package /tgac/software/testing/bin/centrifuge-1.0.4_6cc874e

snakemake -j 100 --snakefile GenerateSamplesheet.smk --cluster-config cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude}"
 
snakemake -j 150 --snakefile SkimSeqQC.smk --cluster-config cluster.json --latency-wait 60 --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude}"

```


Example `config.yaml` file:

```
run_name: PLATE_NAME
samplesheet: samplesheet.tsv
fastq_dir: "READS"

# Minimum gDNA read count to be passed to SPAdes
min_gDNA_read_count: 70000

# Minimum cDNA read count to be passed to Trinity
min_cDNA_read_count: 50000

# Minimum rRNA length to report in top blast hits vs pr2 DB
# This is the length of the blast alignment not the annotated gene
min_rRNA_blast_length: 500

# Interval to use for calculating random kmer uniquess (not currently reported in summary)
rku_gDNA_interval: 25000
rku_cDNA_interval: 25000

# Identity threshold for clustering transcriptome assembly
cdhit_identity_threshold: 0.98

# How many reads to classify using centrifuge [-u/--upto <int>   stop after first <int> reads/pairs (no limit)]
centrifuge_upto: 1000000

# Paths to databases
bbduk_adapters: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/bbmap/resources/adapters.fa
pr2_database: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/databases/pr2/4.14.0/pr2_version_4.14.0_SSU_taxo_long.fasta
diamond_database: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/databases/diamond/reference_proteomes.dmnd
taxonomy_database: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/databases/taxdump/fullnamelineage.dmp
kraken2_database: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/databases/kraken2/
centrifuge_database: /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Reference/databases/centrifuge/centrifugedb
centrifuge_NT_database: /ei/public/databases/centrifuge/ncbi/nt_20200707/nt
```
