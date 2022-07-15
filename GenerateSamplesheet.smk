# Snakemake pipeline to generate a samplesheet with read counts from G&T-Seq skim sequencing data
# Expects to be run in a directory containing a "READS" directory containing gDNA and cDNA reads in the format:
# {sample}_cDNA_R[12].fastq.gz and {sample}_gDNA_R[12].fastq.gz

import sys
from os.path import join

configfile: "config.yaml"

FASTQ_DIR = config["fastq_dir"]

ALL_SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,.+[^/]+}_cDNA_R1.fastq.gz'))

rule generate_samplesheet:
    input:
        fastp_json_files = expand(join(FASTQ_DIR, '{sample}_{type}.fastp.json'),sample=ALL_SAMPLES, type=["gDNA", "cDNA"])
    output:
        samplesheet = config["samplesheet"]
    params:
        samples = expand("{sample}", sample=ALL_SAMPLES),
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
    shell: "source fastp-0.23.2 && fastp --in1 {input.raw_read1} --in2 {input.raw_read2} --thread {threads} -j {output.fastp_json} -h {output.fastp_html} > {log} 2>&1"
