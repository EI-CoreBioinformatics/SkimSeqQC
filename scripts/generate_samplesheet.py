import json
import sys
from os import path

# Check if a samplesheet already exists (snakemake will check so shouldn't be the case)
if path.exists(snakemake.output.samplesheet):
    print("Error,", snakemake.output.samplesheet, "already exists")
    sys.exit(-1)

FASTQ_DIR = snakemake.params.fastq_dir
ALL_SAMPLES = sorted(snakemake.params.samples)

fo = open(snakemake.output.samplesheet, "w")

samples = {}

# Read in fastp json files to extract number of reads
for sample in ALL_SAMPLES:
    samples[sample] = {}

    # Get gDNA read count
    if sample in snakemake.params.gDNA_samples:
        data = ""

        with open(path.join(FASTQ_DIR, sample + "_gDNA.fastp.json"), "r") as f:
            data = json.load(f)

        samples[sample]["gDNA_read_count"] = data["summary"]["before_filtering"]["total_reads"]

    else:
        samples[sample]["gDNA_read_count"] = "NA"

    # Get cDNA read count
    if sample in snakemake.params.cDNA_samples:
        data = ""

        with open(path.join(FASTQ_DIR, sample + "_cDNA.fastp.json"), "r") as f:
            data = json.load(f)

        samples[sample]["cDNA_read_count"] = data["summary"]["before_filtering"]["total_reads"]

    else:
        samples[sample]["cDNA_read_count"] = "NA"

fo.write("\t".join(["Sample", "gDNA_Raw_Reads", "cDNA_Raw_Reads"]) + "\n")

for sample in ALL_SAMPLES:
    fo.write("\t".join(map(str, [sample, samples[sample]["gDNA_read_count"], samples[sample]["cDNA_read_count"]])) + "\n")

fo.close()
