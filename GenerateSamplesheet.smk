# Snakemake pipeline to generate a samplesheet with read counts from G&T-Seq skim sequencing data
# Expects to be run in a directory containing a "READS" directory containing gDNA and cDNA reads in the format:
# {sample}_cDNA_R[12].fastq.gz and {sample}_gDNA_R[12].fastq.gz

import sys
from os.path import join
from os import listdir

configfile: "config.yaml"

FASTQ_DIR = config["fastq_dir"]

ALL_SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,.+[^/]+}_cDNA_R1.fastq.gz'))

cDNA_SAMPLES = [f.split("_cDNA_R1.fastq.gz")[0] for f in listdir(FASTQ_DIR) if f.endswith("_cDNA_R1.fastq.gz")]
gDNA_SAMPLES = [f.split("_gDNA_R1.fastq.gz")[0] for f in listdir(FASTQ_DIR) if f.endswith("_gDNA_R1.fastq.gz")]

ALL_SAMPLES = set(cDNA_SAMPLES + gDNA_SAMPLES)

rule generate_samplesheet:
    input:
        fastp_cDNA = expand(join(FASTQ_DIR, '{sample}_{type}.fastp.json'), sample = cDNA_SAMPLES, type = ["cDNA"]),
        fastp_gDNA = expand(join(FASTQ_DIR, '{sample}_{type}.fastp.json'), sample = gDNA_SAMPLES, type = ["gDNA"])
    output:
        samplesheet = config["samplesheet"]
    params:
        samples = ALL_SAMPLES,
        cDNA_samples = cDNA_SAMPLES,
        gDNA_samples = gDNA_SAMPLES,
        fastq_dir = FASTQ_DIR
    script:
        "scripts/generate_samplesheet.py"

rule fastp:
    input:
        raw_read1 = join(FASTQ_DIR, "{sample}_{type}_R1.fastq.gz"),
        raw_read2 = join(FASTQ_DIR, "{sample}_{type}_R2.fastq.gz")
    output:
        fastp_json = join(FASTQ_DIR, "{sample}_{type}.fastp.json"),
        fastp_html = join(FASTQ_DIR, "{sample}_{type}.fastp.html")
    threads: 1
    log: "logs/fastp_{sample}_{type}.log"
    benchmark: "benchmarks/fastp_{sample}_{type}.tsv"
    shell: "fastp --in1 {input.raw_read1} --in2 {input.raw_read2} --thread {threads} -j {output.fastp_json} -h {output.fastp_html} > {log} 2>&1"
