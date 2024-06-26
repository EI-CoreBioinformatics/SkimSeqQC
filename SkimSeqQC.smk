# Snakemake pipeline to QC G&T-Seq skim-sequencing data.
# Expects to be run in a directory containing a "READS" directory containing gDNA and cDNA reads in the format:
# {sample}_cDNA_R[12].fastq.gz and {sample}_gDNA_R[12].fastq.gz
# Working directory should include a samplesheet generated by GenerateSamplesheet.smk pipeline.

import sys
from os.path import join, basename
from os import getcwd
from re import split

configfile: "config.yaml"

FASTQ_DIR = config["fastq_dir"]

ALL_SAMPLES = []
gDNA_SAMPLES = []
cDNA_SAMPLES = []

with open(config["samplesheet"], "r") as f:
    f.readline() # Ignore header line
    for line in f:
        line = line.strip().split("\t")
        sample = line[0]

        ALL_SAMPLES.append(sample)

        gDNA_read_count = line[1]
        cDNA_read_count = line[2]

        if gDNA_read_count != "NA":
            gDNA_read_count = int(gDNA_read_count)

            if gDNA_read_count >= config["min_gDNA_read_count"]:
                gDNA_SAMPLES.append(sample)

        if cDNA_read_count != "NA":
            cDNA_read_count = int(cDNA_read_count)

            if cDNA_read_count >= config["min_cDNA_read_count"]:
                cDNA_SAMPLES.append(sample)

summary_filename = config["run_name"] + "_summary.tsv"
readme_filename = config["run_name"] + "_readme.txt"
rRNA_summary_filename = config["run_name"] + "_rRNA_summary.tsv"
# multiqc_report_filename = config["run_name"] + "_multiqc_report.html" # TODO


# Write info and config params to a readme file, only if not run in dry mode
if "-n" not in sys.argv and "--dry-run" not in sys.argv:
    fo = open(readme_filename, "w")
    fo.write("wd = " + getcwd() + "\n")
    fo.write("\nAll samples (" + str(len(ALL_SAMPLES)) + "): " + ", ".join(map(str, sorted(ALL_SAMPLES))) + "\n")
    fo.write("\ngDNA samples (" + str(len(gDNA_SAMPLES)) + "): " + ", ".join(map(str, sorted(gDNA_SAMPLES))) + "\n")
    fo.write("\ncDNA samples (" + str(len(cDNA_SAMPLES)) + "): " + ", ".join(map(str, sorted(cDNA_SAMPLES))) + "\n")
    fo.write("\n")

    excluded_gDNA_samples = list(set(ALL_SAMPLES) - set(gDNA_SAMPLES))
    excluded_cDNA_samples = list(set(ALL_SAMPLES) - set(cDNA_SAMPLES))

    fo.write("gDNA samples excluded due to having < " + str(config["min_gDNA_read_count"]) + " reads (" + str(len(excluded_gDNA_samples)) + "): " + ", ".join(sorted(excluded_gDNA_samples)) + "\n")
    fo.write("cDNA samples excluded due to having < " + str(config["min_cDNA_read_count"]) + " reads (" + str(len(excluded_cDNA_samples)) + "): " + ", ".join(sorted(excluded_cDNA_samples)) + "\n")
    fo.write("\n======================\n")
    fo.write("Config file parameters")
    fo.write("\n======================\n")
    for c in config:
        fo.write(c + ": " + str(config[c]) + "\n")
    fo.close()


TARGETS = [summary_filename]

if config["map_gDNA"]:
    TARGETS.append(expand(join("{sample}", "{sample}_gDNA_alignment.flagstat"), sample=gDNA_SAMPLES))

# Can only run hisat2 on samples that meet min cDNA read count AND min gDNA read count
if config["map_cDNA"]:
    for sample in cDNA_SAMPLES:
        if sample in gDNA_SAMPLES:
            TARGETS.append(expand(join("{sample}", "{sample}.hisat2_summary.txt"), sample=[sample]))

rule all:
    input:
        TARGETS

rule bbduk_gDNA:
    input: 
        r1 = join(FASTQ_DIR, "{sample}_gDNA_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_gDNA_R2.fastq.gz")
    output:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R2.trimmed.fastq.gz")
    params:
        adapters = config["bbduk_adapters"]
    threads: 2
    log: "logs/bbduk_gDNA_{sample}.log"
    benchmark: "benchmarks/bbduk_gDNA_{sample}.tsv"
    shell: "bbduk.sh ref={params.adapters} ktrim=r k=21 mink=7 hdist=1 qtrim=lr trimq=10 maq=20 minlength=50 tpe tbo in1={input.r1} in2={input.r2} out1={output.trimmed_r1} out2={output.trimmed_r2} t={threads} > {log} 2>&1"

rule bbduk_cDNA:
    input: 
        r1 = join(FASTQ_DIR, "{sample}_cDNA_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_cDNA_R2.fastq.gz")
    output:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R2.trimmed.fastq.gz")
    params:
        adapters = config["bbduk_adapters"]
    threads: 2
    log: "logs/bbduk_cDNA_{sample}.log"
    benchmark: "benchmarks/bbduk_cDNA_{sample}.tsv"
    shell: "bbduk.sh t={threads} in1={input.r1} in2={input.r2} out=stdout.fq minlength=50 ktrim=l k=14 mink=5 hdist=1 hdist2=0 rcomp=f tbo literal=AAGCAGTGGTATCAACGCAGAGT,TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG,GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG,GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG 2> {log} | bbduk.sh t={threads} in=stdin.fq int=t out1={output.trimmed_r1} out2={output.trimmed_r2} minlength=50 ktrim=r k=14 mink=5 hdist=1 hdist2=0 rcomp=f qtrim=lr trimq=10 maq=20 tbo tpe trimpolygright=5 trimpolya=5 literal=CTGTCTCTTATACACATCTGACGCTGCCGACGA,CTGTCTCTTATACACATCTCCGAGCCCACGAGAC,CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC,CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC,ATCTCGTATGCCGTCTTCTGCTTG,ACTCTGCGTTGATACCACTGCTT,AGATCGGAAGAGCACACG,TGGAATTCTCGGGTGCCAAGG 2>>{log}"

rule fastp_trimmed_reads:
    input:
        r1_trimmed = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),
        r2_trimmed = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz")
    output:
        fastp_json = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_trimmed.fastp.json"),
        fastp_html = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_trimmed.fastp.html")
    threads: 1
    log: "logs/fastp_trimmed_{sample}_{sampletype}.log"
    benchmark: "benchmarks/fastp_trimmed_{sample}_{sampletype}.tsv"
    shell: "fastp --in1 {input.r1_trimmed} --in2 {input.r2_trimmed} --thread {threads} -j {output.fastp_json} -h {output.fastp_html} > {log} 2>&1"

rule make_centrifuge_samplesheets:
    input:
        trimmed_r1_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        trimmed_r2_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        trimmed_r1_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        trimmed_r2_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz"),sample=gDNA_SAMPLES, sampletype=["gDNA"])
    output:
        centrifuge_gDNA_samplesheet = "samplesheet_centrifuge_gDNA.tsv",
        centrifuge_cDNA_samplesheet = "samplesheet_centrifuge_cDNA.tsv",
    run:    
        fo_cDNA = open(output.centrifuge_cDNA_samplesheet, "w")

        for r1, r2 in zip(input.trimmed_r1_cDNA, input.trimmed_r2_cDNA):
            file_basename = basename(r1)
            sample = split("_[cg]DNA", file_basename)[0]
            sample_type = "cDNA"

            centrifuge_output = "centrifuge/" + sample + "_"+ sample_type + "_centrifuge_classification.out"
            centrifuge_report = "centrifuge/" + sample + "_"+ sample_type + "_centrifuge_report.txt"
            line = ["2", r1, r2, centrifuge_output, centrifuge_report]

            fo_cDNA.write("\t".join(line) + "\n")
        fo_cDNA.close()

        fo_gDNA = open(output.centrifuge_gDNA_samplesheet, "w")
        for r1, r2 in zip(input.trimmed_r1_gDNA, input.trimmed_r2_gDNA):
            file_basename = basename(r1)
            sample = split("_[cg]DNA", file_basename)[0]
            sample_type = "gDNA"

            centrifuge_output = "centrifuge/" + sample + "_"+ sample_type + "_centrifuge_classification.out"
            centrifuge_report = "centrifuge/" + sample + "_"+ sample_type + "_centrifuge_report.txt"
            line = ["2", r1, r2, centrifuge_output, centrifuge_report]

            fo_gDNA.write("\t".join(line) + "\n")
        fo_gDNA.close()

rule centrifuge_gDNA:
    input:
        centrifuge_gDNA_samplesheet = "samplesheet_centrifuge_gDNA.tsv"
    output:
        centrifuge_gDNA_done = join("centrifuge", "centrifuge_gDNA.done"),
        centrifuge_classification_out = expand(join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.out"), sample=gDNA_SAMPLES, sampletype=["gDNA"])
    params:
        centrifuge_database = config["centrifuge_database"],
        centrifuge_upto = config["centrifuge_upto"]
    threads: 16
    log: "logs/centrifuge_gDNA.log"
    benchmark: "benchmarks/centrifuge_gDNA.tsv"
    shell: "centrifuge --threads {threads} -x {params.centrifuge_database} --upto {params.centrifuge_upto} --sample-sheet {input.centrifuge_gDNA_samplesheet} 2> {log} && touch {output.centrifuge_gDNA_done}"

rule centrifuge_cDNA:
    input:
        centrifuge_cDNA_samplesheet = "samplesheet_centrifuge_cDNA.tsv"
    output:
        centrifuge_cDNA_done = join("centrifuge", "centrifuge_cDNA.done"),
        centrifuge_classification_out = expand(join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.out"), sample=cDNA_SAMPLES, sampletype=["cDNA"])
    params:
        centrifuge_database = config["centrifuge_database"],
        centrifuge_upto = config["centrifuge_upto"]
    threads: 16
    log: "logs/centrifuge_cDNA.log"
    benchmark: "benchmarks/centrifuge_cDNA.tsv"
    shell: "centrifuge --threads {threads} -x {params.centrifuge_database} --upto {params.centrifuge_upto} --sample-sheet {input.centrifuge_cDNA_samplesheet} 2> {log} && touch {output.centrifuge_cDNA_done}"

rule centrifuge_kreport:
    input:
        centrifuge_done = join("centrifuge", "centrifuge_{sampletype}.done"),
        centrifuge_classification_out = join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.out")
    output:
        centrifuge_kreport = join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.kreport"),
        centrifuge_summary = join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.summary.txt")
    params:
        centrifuge_database = config["centrifuge_database"]
    log: "logs/centrifuge_kreport_{sample}_{sampletype}.log"
    benchmark: "benchmarks/centrifuge_kreport_{sample}_{sampletype}.tsv"
    threads: 1
    shell: "centrifuge-kreport -x {params.centrifuge_database} {input.centrifuge_classification_out} > {output.centrifuge_kreport} 2> {log} && parse_kreport.py {output.centrifuge_kreport} > {output.centrifuge_summary} 2>> {log}"

rule make_centrifuge_NT_samplesheets:
    input:
        trimmed_r1_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        trimmed_r2_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        trimmed_r1_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        trimmed_r2_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz"),sample=gDNA_SAMPLES, sampletype=["gDNA"])
    output:
        centrifuge_NT_gDNA_samplesheet = "samplesheet_centrifuge_NT_gDNA.tsv",
        centrifuge_NT_cDNA_samplesheet = "samplesheet_centrifuge_NT_cDNA.tsv",
    run:    
        fo_cDNA = open(output.centrifuge_NT_cDNA_samplesheet, "w")

        for r1, r2 in zip(input.trimmed_r1_cDNA, input.trimmed_r2_cDNA):
            file_basename = basename(r1)
            sample = split("_[cg]DNA", file_basename)[0]
            sample_type = "cDNA"

            centrifuge_NT_output = "centrifuge_NT/" + sample + "_"+ sample_type + "_centrifuge_NT_classification.out"
            centrifuge_NT_report = "centrifuge_NT/" + sample + "_"+ sample_type + "_centrifuge_NT_report.txt"
            line = ["2", r1, r2, centrifuge_NT_output, centrifuge_NT_report]

            fo_cDNA.write("\t".join(line) + "\n")
        fo_cDNA.close()

        fo_gDNA = open(output.centrifuge_NT_gDNA_samplesheet, "w")

        for r1, r2 in zip(input.trimmed_r1_gDNA, input.trimmed_r2_gDNA):
            file_basename = basename(r1)
            sample = split("_[cg]DNA", file_basename)[0]
            sample_type = "gDNA"

            centrifuge_NT_output = "centrifuge_NT/" + sample + "_"+ sample_type + "_centrifuge_NT_classification.out"
            centrifuge_NT_report = "centrifuge_NT/" + sample + "_"+ sample_type + "_centrifuge_NT_report.txt"
            line = ["2", r1, r2, centrifuge_NT_output, centrifuge_NT_report]

            fo_gDNA.write("\t".join(line) + "\n")

        fo_gDNA.close()

rule centrifuge_NT_gDNA:
    input:
        centrifuge_NT_gDNA_samplesheet = "samplesheet_centrifuge_NT_gDNA.tsv"
    output:
        centrifuge_NT_gDNA_done = join("centrifuge_NT", "centrifuge_NT_gDNA.done"),
        centrifuge_NT_classification_out = expand(join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.out"), sample=gDNA_SAMPLES, sampletype=["gDNA"])
    params:
        centrifuge_NT_database = config["centrifuge_NT_database"],
        centrifuge_upto = config["centrifuge_upto"]
    threads: 10
    log: "logs/centrifuge_NT_gDNA.log"
    benchmark: "benchmarks/centrifuge_NT_gDNA.tsv"
    shell: "centrifuge --threads {threads} -x {params.centrifuge_NT_database} --upto {params.centrifuge_upto} --sample-sheet {input.centrifuge_NT_gDNA_samplesheet} 2> {log} && touch {output.centrifuge_NT_gDNA_done}"

rule centrifuge_NT_cDNA:
    input:
        centrifuge_NT_cDNA_samplesheet = "samplesheet_centrifuge_NT_cDNA.tsv"
    output:
        centrifuge_NT_cDNA_done = join("centrifuge_NT", "centrifuge_NT_cDNA.done"),
        centrifuge_NT_classification_out = expand(join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.out"), sample=cDNA_SAMPLES, sampletype=["cDNA"])
    params:
        centrifuge_NT_database = config["centrifuge_NT_database"],
        centrifuge_upto = config["centrifuge_upto"]
    threads: 10
    log: "logs/centrifuge_NT_cDNA.log"
    benchmark: "benchmarks/centrifuge_NT_cDNA.tsv"
    shell: "centrifuge --threads {threads} -x {params.centrifuge_NT_database} --upto {params.centrifuge_upto} --sample-sheet {input.centrifuge_NT_cDNA_samplesheet} 2> {log} && touch {output.centrifuge_NT_cDNA_done}"

rule centrifuge_NT_kreport:
    input:
        centrifuge_NT_done = join("centrifuge_NT", "centrifuge_NT_{sampletype}.done"),
        centrifuge_NT_classification_out = join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.out")
    output:
        centrifuge_NT_kreport = join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.kreport"),
        centrifuge_NT_summary = join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.summary.txt")
    params:
        centrifuge_NT_database = config["centrifuge_NT_database"]
    log: "logs/centrifuge_NT_kreport_{sample}_{sampletype}.log"
    benchmark: "benchmarks/centrifuge_NT_kreport_{sample}_{sampletype}.tsv"
    threads: 1
    shell: "centrifuge-kreport -x {params.centrifuge_NT_database} {input.centrifuge_NT_classification_out} > {output.centrifuge_NT_kreport} 2> {log} && parse_kreport.py {output.centrifuge_NT_kreport} > {output.centrifuge_NT_summary} 2>> {log}"

rule rku_gDNA:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R2.trimmed.fastq.gz")
    output:
        rku_gDNA = join("{sample}", "{sample}_rku_gDNA.tsv")
    params:
        interval = config["rku_gDNA_interval"]
    log: "logs/bbcountunique_gDNA_{sample}.log"
    benchmark: "benchmarks/bbcountunique_gDNA_{sample}.tsv"
    shell: "bbcountunique.sh in1={input.trimmed_r1} in2={input.trimmed_r2} out={output.rku_gDNA} interval={params.interval} percent=t > {log} 2>&1"

rule rku_cDNA:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R2.trimmed.fastq.gz")
    output:
        rku_cDNA = join("{sample}", "{sample}_rku_cDNA.tsv")
    params:
        interval = config["rku_cDNA_interval"]
    log: "logs/bbcountunique_cDNA_{sample}.log"
    benchmark: "benchmarks/bbcountunique_cDNA_{sample}.tsv"
    shell: "bbcountunique.sh in1={input.trimmed_r1} in2={input.trimmed_r2} out={output.rku_cDNA} interval={params.interval} percent=t > {log} 2>&1"

rule fastqc:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz")
    output:
        trimmed_r1_qc = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed_fastqc.html"),
        trimmed_r2_qc = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed_fastqc.html")
    threads: 2
    log: "logs/fastqc_{sampletype}_{sample}.log"
    benchmark: "benchmarks/fastqc_{sampletype}_{sample}.tsv"
    shell: "fastqc -t {threads} {input} > {log} 2>&1"

rule bbmerge:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz")
    output:
        merged_reads = join("{sample}", "merged_reads", "{sample}_{sampletype}.merged.fastq.gz"),
        ihist = join("{sample}", "merged_reads", "{sample}_{sampletype}.ihist.txt")
    threads: 2
    log: "logs/bbmerge_{sampletype}_{sample}.log"
    benchmark: "benchmarks/bbmerge_{sampletype}_{sample}.tsv"
    shell: "bbmerge.sh in1={input.trimmed_r1} in2={input.trimmed_r2} out={output.merged_reads} ihist={output.ihist} t={threads} > {log} 2>&1"

# TODO UPDATE KRAKEN2 / move parse_kreport.py python script
rule kraken2:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed.fastq.gz")
    output:
        kraken_report = join("kraken2", "{sample}_{sampletype}_kraken.kreport"),
        kraken_classification = join("kraken2", "{sample}_{sampletype}_kraken.krk"),
        kraken_summary = join("kraken2", "{sample}_{sampletype}_kraken.summary")
    params:
        kraken2_database = config["kraken2_database"]
    threads: 8
    log: "logs/kraken2_{sampletype}_{sample}.log"
    benchmark: "benchmarks/kraken2_{sampletype}_{sample}.tsv"
    shell: "kraken2 --threads {threads} --report {output.kraken_report} --db {params.kraken2_database} --paired {input.trimmed_r1} {input.trimmed_r2} > {output.kraken_classification} 2> {log} && parse_kreport.py {output.kraken_report} > {output.kraken_summary}"

rule spades:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R2.trimmed.fastq.gz")
    output:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta")
    params:
        out_dir = directory(join("{sample}", "assembly")),
        original_contigs = join("{sample}", "assembly", "contigs.fasta"),
        original_scaffolds = join("{sample}", "assembly", "scaffolds.fasta"),
        new_contigs = join("{sample}", "assembly", "{sample}.contigs.fasta"),
        new_scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta")
    threads: 8
    log: "logs/spades_{sample}.log"
    benchmark: "benchmarks/spades_{sample}.tsv"
    shell: "spades.py --sc -1 {input.trimmed_r1} -2 {input.trimmed_r2} --threads {threads} -o {params.out_dir} > {log} 2>&1 && mv {params.original_contigs} {params.new_contigs} && mv {params.original_scaffolds} {params.new_scaffolds}"

rule seqtk_seq:
    input:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta")
    output:
        filtered_scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.1000.fasta")
    params:
        min_length = 1000
    threads: 1
    benchmark: "benchmarks/seqtk_seq_{sample}.tsv"
    shell: "seqtk seq -l 70 -L {params.min_length} {input.scaffolds} > {output.filtered_scaffolds}"

rule quast:
    input:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta"),
        filtered_scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.1000.fasta")
    output:
        quast_report = join("{sample}", "quast", "report.txt")
    params:
        out_dir = directory(join("{sample}", "quast"))
    threads: 1
    log: "logs/quast_{sample}.log"
    benchmark: "benchmarks/quast_{sample}.tsv"
    shell: "quast -t {threads} -o {params.out_dir} {input.scaffolds} {input.filtered_scaffolds} > {log} 2>&1"

rule minimap2_gDNA:
    input:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta"),
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_gDNA_R2.trimmed.fastq.gz")
    output:
        gDNA_alignment = join("{sample}", "{sample}_gDNA_alignment.bam"),
        gDNA_alignment_flagstat = join("{sample}", "{sample}_gDNA_alignment.flagstat"),
    threads: 4
    log: "logs/minimap2_gDNA_{sample}.log"
    benchmark: "benchmarks/minimap2_gDNA_{sample}.tsv"
    shell: "minimap2 -t {threads} -ax sr {input.scaffolds} {input.trimmed_r1} {input.trimmed_r2} | samtools sort -@ {threads} -o {output.gDNA_alignment} && samtools index -@ {threads} {output.gDNA_alignment} && samtools flagstat -@ {threads} {output.gDNA_alignment} > {output.gDNA_alignment_flagstat}"

rule hisat2_build:
    input:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta")
    output:
        hisat2_index_done = join("{sample}", "assembly", "hisat2_build.done")
    threads: 1
    log: "logs/hisat2_build_{sample}.log"
    benchmark: "benchmarks/hisat2_build_{sample}.tsv"
    shell: "hisat2-build -p {threads} {input.scaffolds} {input.scaffolds} && touch {output.hisat2_index_done}"

rule hisat2:
    input:
        hisat2_index_done = join("{sample}", "assembly", "hisat2_build.done"),
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta"),
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R2.trimmed.fastq.gz")
    output:
        cDNA_alignment = join("{sample}", "{sample}_cDNA_alignment.bam"),
        hisat2_summary = join("{sample}", "{sample}.hisat2_summary.txt")
    threads: 4
    log: "logs/hisat2_{sample}.log"
    benchmark: "benchmarks/hisat2_{sample}.tsv"
    shell: "hisat2 -q -p {threads} --dta --summary-file {output.hisat2_summary} -t -x {input.scaffolds} -1 {input.trimmed_r1} -2 {input.trimmed_r2} | samtools sort -@ {threads} -O BAM -o {output.cDNA_alignment} && samtools index -@ {threads} {output.cDNA_alignment}"

rule barrnap_gDNA:
    input:
        scaffolds = join("{sample}", "assembly", "{sample}.scaffolds.fasta")
    output:
        rrna = join("{sample}", "rrna", "{sample}_gDNA.rrna.fasta")
    threads: 2
    log: "logs/barrnap_gDNA_{sample}.log"
    benchmark: "benchmarks/barrnap_gDNA_{sample}.tsv"
    shell: "barrnap --kingdom euk --threads {threads} --outseq {output.rrna} {input.scaffolds} > {log} 2>&1"

# TODO move print_best_blast_hits.py python script
rule pr2_blast_gDNA:
    input:
        rrna = join("{sample}", "rrna", "{sample}_gDNA.rrna.fasta")
    output:
        hits = join("{sample}", "rrna", "{sample}_gDNA.rrna.blast.tsv"),
        top_hits = join("{sample}", "rrna", "{sample}_gDNA.rrna.blast.top.tsv")
    params:
        db = config["pr2_database"]
    threads: 2
    benchmark: "benchmarks/pr2_blast_gDNA_{sample}.tsv"
    shell: "blastn -query {input.rrna} -db {params.db} -outfmt 6 -out {output.hits} -num_threads {threads} && print_best_blast_hits.py {output.hits} > {output.top_hits}"

rule trinity:
    input:
        trimmed_r1 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R1.trimmed.fastq.gz"),
        trimmed_r2 = join("{sample}", "trimmed_reads", "{sample}_cDNA_R2.trimmed.fastq.gz")
    output:
        transcriptome = join("{sample}", "trinity", "{sample}.Trinity.fasta")
    params:
        output_dir = directory(join("{sample}", "trinity")),
        original_output = join("{sample}", "trinity.Trinity.fasta"),
        new_output = join("{sample}", "trinity", "{sample}.Trinity.fasta")
    threads: 8
    log: "logs/trinity_{sample}.log"
    benchmark: "benchmarks/trinity_{sample}.tsv"
    shell:
        """
        Trinity --full_cleanup --seqType fq --max_memory 70G --left {input.trimmed_r1} --right {input.trimmed_r2} --CPU {threads} --output {params.output_dir} > {log} 2>&1
        rm -rf {params.output_dir}
        mkdir -p {params.output_dir}
        mv {params.original_output} {params.new_output}
        """

rule cdhit:
    input:
        transcriptome = join("{sample}", "trinity", "{sample}.Trinity.fasta")
    output:
        transcriptome_cdhit = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta")
    params:
        identity_threshold = config["cdhit_identity_threshold"]
    threads: 4
    log: "logs/cdhit_{sample}.log"
    benchmark: "benchmarks/cdhit_{sample}.log"
    shell: "cd-hit-est -o {output.transcriptome_cdhit} -c {params.identity_threshold} -i {input.transcriptome} -p 1 -d 0 -b 3 -T {threads}  > {log} 2>&1"    
    
rule barrnap_cDNA:
    input:
        transcriptome_cdhit = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta")
    output:
        rrna = join("{sample}", "rrna", "{sample}_cDNA.rrna.fasta")
    threads: 2
    log: "logs/barrnap_cDNA_{sample}.log"
    benchmark: "benchmarks/barrnap_cDNA_{sample}.tsv"
    shell: "barrnap --kingdom euk --threads {threads} --outseq {output.rrna} {input.transcriptome_cdhit} > {log} 2>&1"

# TODO move print_best_blast_hits.py python script
rule pr2_blast_cDNA:
    input:
        rrna = join("{sample}", "rrna", "{sample}_cDNA.rrna.fasta")
    output:
        hits = join("{sample}", "rrna", "{sample}_cDNA.rrna.blast.tsv"),
        top_hits = join("{sample}", "rrna", "{sample}_cDNA.rrna.blast.top.tsv")
    params:
        db = config["pr2_database"]
    threads: 2
    benchmark: "benchmarks/pr2_blast_cDNA_{sample}.tsv"
    shell: "blastn -query {input.rrna} -db {params.db} -outfmt 6 -out {output.hits} -num_threads {threads} && print_best_blast_hits.py {output.hits} > {output.top_hits}"

rule transdecoder_LongOrfs:
    input:
        transcriptome_cdhit = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta")
    output:
        orfs = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder_dir", "longest_orfs.pep")
    params:
        output_dir = directory(join("{sample}", "trinity"))
    log: "logs/transdecoder_LongOrfs_{sample}.log"
    benchmark: "benchmarks/transdecoder_LongOrfs_{sample}.tsv"
    shell: "TransDecoder.LongOrfs -t {input.transcriptome_cdhit} -m 60 --output_dir {params.output_dir} > {log} 2>&1"

rule transdecoder_Predict:
    input:
        transcriptome_cdhit = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta"),
	orfs = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder_dir", "longest_orfs.pep")
    output:
        proteins = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.pep")
    params:
        output_dir = directory(join("{sample}", "trinity")),
        original_output = "{sample}.Trinity.cdhit.fasta.transdecoder.*",
        new_output = join("{sample}", "trinity")
    log: "logs/transdecoder_Predict_{sample}.log"
    benchmark: "benchmarks/transdecoder_Predict_{sample}.tsv"
    shell: "TransDecoder.Predict -t {input.transcriptome_cdhit} --output_dir {params.output_dir} > {log}"

rule diamond_transcriptome:
    input:
        proteins = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.pep")
    output:
        classifications = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.diamond.out")
    params:
        diamond_database = config["diamond_database"]
    threads: 8
    log: "logs/diamond_transcriptome_{sample}.log"
    benchmark: "benchmarks/diamond_transcriptome_{sample}.log"
    shell: "diamond blastp --query {input.proteins} --outfmt 102 --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 --db {params.diamond_database} > {output.classifications}"

# TODO move map_taxids.py python script
rule map_taxids:
    input:
        classifications = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.diamond.out")
    output:
        mapped_classifications = join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.diamond.mapped.out")
    params:
        taxonomy_db = config["taxonomy_database"]
    threads: 1
    log: "logs/map_taxids_{sample}.log"
    benchmark: "benchmarks/map_taxids_{sample}.tsv"
    shell: "map_taxids.py {params.taxonomy_db} {input.classifications} > {output.mapped_classifications}"

# This rule summarises all 18S rRNA hits in either the genome or transcriptome, reporting which samples each hit is 
# found in. Limitation to this is that we are only looking at 18S genes (barrnap --kingdom euk & pr2 database) but 
# there will still be some prokaryotic genes annotated and only looking at the best hit. Recovered genes are likely
# to be partial whereas in a coassembly might be more complete. Organellar rRNAs are poorly represented in databases.
# To report summary of ALL samples (gDNA+cDNA), gDNA only, and cDNA only.
rule rRNA_summary:
    input:
        top_rrna_hits_gDNA = expand(join("{sample}", "rrna", "{sample}_gDNA.rrna.blast.top.tsv"), sample=gDNA_SAMPLES),
        top_rrna_hits_cDNA = expand(join("{sample}", "rrna", "{sample}_cDNA.rrna.blast.top.tsv"), sample=cDNA_SAMPLES)
    output:
        rRNA_summary = rRNA_summary_filename
    params:
        min_rRNA_blast_length = int(config["min_rRNA_blast_length"])
    run:
        fo = open(output.rRNA_summary, "w")
        samples_per_top_hit = {}
        samples_per_top_hit_gDNA = {}
        samples_per_top_hit_cDNA = {}

        samples = []

        for i in input.top_rrna_hits_gDNA:
            sample = basename(i).split("_gDNA")[0]
            samples.append(sample + "_gDNA")

            with open(i, "r") as f:
                for line in f:
                    line = line.strip().split("\t")
                    query, hit, identity, aln_len = line[0], line[1], line[2], line[3]

                    if "18S" in query or "16S" in query:
                        if int(aln_len) > params.min_rRNA_blast_length:
                            if hit not in samples_per_top_hit:
                                samples_per_top_hit[hit] = set()
                            samples_per_top_hit[hit].add(sample)

                            if hit not in samples_per_top_hit_gDNA:
                                samples_per_top_hit_gDNA[hit] = set()
                            samples_per_top_hit_gDNA[hit].add(sample)

        for i in input.top_rrna_hits_cDNA:
            sample = basename(i).split("_cDNA")[0]
            samples.append(sample + "_cDNA")

            with open(i, "r") as f:
                for line in f:
                    line = line.strip().split("\t")
                    query, hit, identity, aln_len = line[0], line[1], line[2], line[3]

                    if "18S" in query or "16S" in query:
                        if int(aln_len) > int(params.min_rRNA_blast_length):
                            if hit not in samples_per_top_hit:
                                samples_per_top_hit[hit] = set()
                            samples_per_top_hit[hit].add(sample)

                            if hit not in samples_per_top_hit_cDNA:
                                samples_per_top_hit_cDNA[hit] = set()
                            samples_per_top_hit_cDNA[hit].add(sample)

        fo.write("Summary of top blast hits of SSU rRNA genes against the pr2 database\n")
        fo.write("Only considering sequences annotated by barrnap as 18S rRNA genes and blast alignments with alignment length >= " + str(params.min_rRNA_blast_length) + "\n\n")
        fo.write("\n")

        fo.write("All gDNA and cDNA samples combined:\n")
        fo.write(", ".join(sorted(samples)))

        fo.write("\n\n")
        fo.write("pr2 gene\tNumber of samples with top hits\tSamples with top hits\n")
        for top_hit in samples_per_top_hit:
            fo.write("\t".join([top_hit, str(len(samples_per_top_hit[top_hit])), "; ".join(sorted(samples_per_top_hit[top_hit]))]) + "\n")

        fo.write("\n\n")
        fo.write("gDNA samples only:\n")
        gDNA_samples = [basename(s).replace("_gDNA.rrna.blast.top.tsv", "") for s in input.top_rrna_hits_gDNA]
        fo.write(", ".join(sorted(gDNA_samples)))
        fo.write("\n\n")
        fo.write("pr2 gene\tNumber of samples with top hits\tSamples with top hits\n")
        for top_hit in samples_per_top_hit_gDNA:
            fo.write("\t".join([top_hit, str(len(samples_per_top_hit_gDNA[top_hit])), "; ".join(sorted(samples_per_top_hit_gDNA[top_hit]))]) + "\n")

        fo.write("\n\n")
        fo.write("cDNA samples only:\n")
        cDNA_samples = [basename(s).replace("_cDNA.rrna.blast.top.tsv", "") for s in input.top_rrna_hits_cDNA]
        fo.write(", ".join(sorted(cDNA_samples)))
        fo.write("\n\n")
        fo.write("pr2 gene\tNumber of samples with top hits\tSamples with top hits\n")
        for top_hit in samples_per_top_hit_cDNA:
            fo.write("\t".join([top_hit, str(len(samples_per_top_hit_cDNA[top_hit])), "; ".join(sorted(samples_per_top_hit_cDNA[top_hit]))]) + "\n")

        fo.close()

rule generate_summary:
    input:
        trimmed_gDNA_R1 = expand(join("{sample}", "trimmed_reads", "{sample}_gDNA_R1.trimmed.fastq.gz"), sample=gDNA_SAMPLES),
        trimmed_gDNA_R2 = expand(join("{sample}", "trimmed_reads", "{sample}_gDNA_R2.trimmed.fastq.gz"), sample=gDNA_SAMPLES),
        trimmed_cDNA_R1 = expand(join("{sample}", "trimmed_reads", "{sample}_cDNA_R1.trimmed.fastq.gz"), sample=cDNA_SAMPLES),
        trimmed_cDNA_R2 = expand(join("{sample}", "trimmed_reads", "{sample}_cDNA_R2.trimmed.fastq.gz"), sample=cDNA_SAMPLES),
        rku_gDNA = expand(join("{sample}", "{sample}_rku_gDNA.tsv"), sample=gDNA_SAMPLES),
        rku_cDNA = expand(join("{sample}", "{sample}_rku_cDNA.tsv"), sample=cDNA_SAMPLES),
        trimmed_r1_qc_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed_fastqc.html"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        trimmed_r2_qc_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed_fastqc.html"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        trimmed_r1_qc_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R1.trimmed_fastqc.html"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        trimmed_r2_qc_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_R2.trimmed_fastqc.html"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        ihist_gDNA = expand(join("{sample}", "merged_reads", "{sample}_{sampletype}.ihist.txt"), sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        ihist_cDNA = expand(join("{sample}", "merged_reads", "{sample}_{sampletype}.ihist.txt"), sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        quast_report = expand(join("{sample}", "quast", "report.txt"), sample=gDNA_SAMPLES),
        top_rrna_hits_gDNA = expand(join("{sample}", "rrna", "{sample}_gDNA.rrna.blast.top.tsv"), sample=gDNA_SAMPLES),
        top_rrna_hits_cDNA = expand(join("{sample}", "rrna", "{sample}_cDNA.rrna.blast.top.tsv"), sample=cDNA_SAMPLES),
        rRNA_summary = rRNA_summary_filename,
        mapped_diamond_classifications = expand(join("{sample}", "trinity", "{sample}.Trinity.cdhit.fasta.transdecoder.diamond.mapped.out"), sample=cDNA_SAMPLES),
        centrifuge_kreport_gDNA = expand(join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.kreport"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        centrifuge_kreport_cDNA = expand(join("centrifuge", "{sample}_{sampletype}_centrifuge_classification.kreport"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        centrifuge_NT_kreport_gDNA = expand(join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.kreport"),sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        centrifuge_NT_kreport_cDNA = expand(join("centrifuge_NT", "{sample}_{sampletype}_centrifuge_NT_classification.kreport"),sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        kraken2_reports_gDNA = expand(join("kraken2", "{sample}_{sampletype}_kraken.kreport"), sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        kraken2_reports_cDNA = expand(join("kraken2", "{sample}_{sampletype}_kraken.kreport"), sample=cDNA_SAMPLES, sampletype=["cDNA"]),
        fastp_json_gDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_trimmed.fastp.json"), sample=gDNA_SAMPLES, sampletype=["gDNA"]),
        fastp_json_cDNA = expand(join("{sample}", "trimmed_reads", "{sample}_{sampletype}_trimmed.fastp.json"), sample=cDNA_SAMPLES, sampletype=["cDNA"])
    output:
        summary = summary_filename
    params:
        list_all_samples = ALL_SAMPLES,
        list_gDNA_samples = gDNA_SAMPLES,
        list_cDNA_samples = cDNA_SAMPLES,
        min_rRNA_blast_length = config["min_rRNA_blast_length"]
    script:
        "scripts/generate_summary.py"
