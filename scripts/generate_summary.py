import json
import sys
from os import path

ALL_SAMPLES = sorted(snakemake.params.list_all_samples)
gDNA_SAMPLES = sorted(snakemake.params.list_gDNA_samples)
cDNA_SAMPLES = sorted(snakemake.params.list_cDNA_samples)

output_file = snakemake.output.summary

samples = {}
for sample in ALL_SAMPLES:
    samples[sample] = {}

# Get raw read counts from fastp json in ./READS directory
for sample in ALL_SAMPLES:
    data = ""
    
    # gDNA raw reads
    if sample in gDNA_SAMPLES:
        with open(path.join("READS", sample + "_gDNA.fastp.json")) as f:
            data = json.load(f)

        total_reads = data["summary"]["before_filtering"]["total_reads"]
        total_bases = data["summary"]["before_filtering"]["total_bases"]

        samples[sample]["gDNA_reads"] = total_reads
        samples[sample]["gDNA_bases"] = total_bases
    else:
        samples[sample]["gDNA_reads"] = "NA"
        samples[sample]["gDNA_bases"] = "NA"

    # cDNA raw reads
    if sample in cDNA_SAMPLES:
        with open(path.join("READS", sample + "_cDNA.fastp.json")) as f:
            data = json.load(f)

        total_reads = data["summary"]["before_filtering"]["total_reads"]
        total_bases = data["summary"]["before_filtering"]["total_bases"]

        samples[sample]["cDNA_reads"] = total_reads
        samples[sample]["cDNA_bases"] = total_bases
    else:
        samples[sample]["cDNA_reads"] = "NA"
        samples[sample]["cDNA_bases"] = "NA"

    # gDNA trimmed reads
    if sample in gDNA_SAMPLES:
        with open(path.join(sample, "trimmed_reads", sample + "_gDNA_trimmed.fastp.json")) as f:
            data = json.load(f)

        total_reads = data["summary"]["before_filtering"]["total_reads"]
        total_bases = data["summary"]["before_filtering"]["total_bases"]
        duplication_rate = data["duplication"]["rate"]

        samples[sample]["gDNA_trimmed_reads"] = total_reads
        samples[sample]["gDNA_trimmed_bases"] = total_bases

        if samples[sample]["gDNA_bases"] == 0:
            samples[sample]["gDNA_trimmed_bases_retained"] = "NA"
        else:
            samples[sample]["gDNA_trimmed_bases_retained"] = (total_bases / samples[sample]["gDNA_bases"]) * 100
        samples[sample]["gDNA_trimmed_duplication_rate"] = duplication_rate

        samples[sample]["gDNA_insert_size_fastp"] = data["insert_size"]["peak"]

        with open(path.join(sample, "merged_reads", sample + "_gDNA.ihist.txt")) as f:
            f.readline() # Ignore first line (mean insert size)
            median_insert_size = f.readline().strip().split("\t")[1]
            samples[sample]["gDNA_insert_size_bbmerge"] = median_insert_size
    else:
        samples[sample]["gDNA_trimmed_reads"] = "NA"
        samples[sample]["gDNA_trimmed_bases"] = "NA"
        samples[sample]["gDNA_trimmed_bases_retained"] = "NA"
        samples[sample]["gDNA_trimmed_duplication_rate"] = "NA"
        samples[sample]["gDNA_insert_size_fastp"] = "NA"
        samples[sample]["gDNA_insert_size_bbmerge"] = "NA"

    # cDNA trimmed reads
    if sample in cDNA_SAMPLES:
        with open(path.join(sample, "trimmed_reads", sample + "_cDNA_trimmed.fastp.json")) as f:
            data = json.load(f)

        total_reads = data["summary"]["before_filtering"]["total_reads"]
        total_bases = data["summary"]["before_filtering"]["total_bases"]
        duplication_rate = data["duplication"]["rate"]

        samples[sample]["cDNA_trimmed_reads"] = total_reads
        samples[sample]["cDNA_trimmed_bases"] = total_bases

        if samples[sample]["cDNA_bases"] == 0:
            samples[sample]["cDNA_trimmed_bases_retained"] = "NA"
        else:
            samples[sample]["cDNA_trimmed_bases_retained"] = (total_bases / samples[sample]["cDNA_bases"]) * 100
        samples[sample]["cDNA_trimmed_duplication_rate"] = duplication_rate

        samples[sample]["cDNA_insert_size_fastp"] = data["insert_size"]["peak"]

        with open(path.join(sample, "merged_reads", sample + "_cDNA.ihist.txt")) as f:
            f.readline() # Ignore first line (mean insert size)
            median_insert_size = f.readline().strip().split("\t")[1]
            samples[sample]["cDNA_insert_size_bbmerge"] = median_insert_size
    else:
        samples[sample]["cDNA_trimmed_reads"] = "NA"
        samples[sample]["cDNA_trimmed_bases"] = "NA"
        samples[sample]["cDNA_trimmed_bases_retained"] = "NA"
        samples[sample]["cDNA_trimmed_duplication_rate"] = "NA"
        samples[sample]["cDNA_insert_size_fastp"] = "NA"
        samples[sample]["cDNA_insert_size_bbmerge"] = "NA"

    # kraken2 gDNA
    if sample in gDNA_SAMPLES:
        with open(path.join("kraken2", sample + "_gDNA_kraken.summary"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["gDNA_kraken2_unclassified"] = unclassified
            samples[sample]["gDNA_kraken2_bacterial"] = bacterial
            samples[sample]["gDNA_kraken2_eukaryotic"] = eukaryotic
            samples[sample]["gDNA_kraken2_top_orders"] = top_orders
    else:
            samples[sample]["gDNA_kraken2_unclassified"] = "NA"
            samples[sample]["gDNA_kraken2_bacterial"] = "NA"
            samples[sample]["gDNA_kraken2_eukaryotic"] = "NA"
            samples[sample]["gDNA_kraken2_top_orders"] = "NA"

    # kraken2 cDNA
    if sample in cDNA_SAMPLES:
        with open(path.join("kraken2", sample + "_cDNA_kraken.summary"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["cDNA_kraken2_unclassified"] = unclassified
            samples[sample]["cDNA_kraken2_bacterial"] = bacterial
            samples[sample]["cDNA_kraken2_eukaryotic"] = eukaryotic
            samples[sample]["cDNA_kraken2_top_orders"] = top_orders
    else:
            samples[sample]["cDNA_kraken2_unclassified"] = "NA"
            samples[sample]["cDNA_kraken2_bacterial"] = "NA"
            samples[sample]["cDNA_kraken2_eukaryotic"] = "NA"
            samples[sample]["cDNA_kraken2_top_orders"] = "NA"

    # Centrifuge gDNA
    if sample in gDNA_SAMPLES:
        with open(path.join("centrifuge", sample + "_gDNA_centrifuge_classification.summary.txt"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["gDNA_centrifuge_unclassified"] = unclassified
            samples[sample]["gDNA_centrifuge_bacterial"] = bacterial
            samples[sample]["gDNA_centrifuge_eukaryotic"] = eukaryotic
            samples[sample]["gDNA_centrifuge_top_orders"] = top_orders
    else:
            samples[sample]["gDNA_centrifuge_unclassified"] = "NA"
            samples[sample]["gDNA_centrifuge_bacterial"] = "NA"
            samples[sample]["gDNA_centrifuge_eukaryotic"] = "NA"
            samples[sample]["gDNA_centrifuge_top_orders"] = "NA"

    # Centrifuge cDNA
    if sample in cDNA_SAMPLES:
        with open(path.join("centrifuge", sample + "_cDNA_centrifuge_classification.summary.txt"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["cDNA_centrifuge_unclassified"] = unclassified
            samples[sample]["cDNA_centrifuge_bacterial"] = bacterial
            samples[sample]["cDNA_centrifuge_eukaryotic"] = eukaryotic
            samples[sample]["cDNA_centrifuge_top_orders"] = top_orders
    else:
            samples[sample]["cDNA_centrifuge_unclassified"] = "NA"
            samples[sample]["cDNA_centrifuge_bacterial"] = "NA"
            samples[sample]["cDNA_centrifuge_eukaryotic"] = "NA"
            samples[sample]["cDNA_centrifuge_top_orders"] = "NA"

    # Centrifuge NT gDNA
    if sample in gDNA_SAMPLES:
        with open(path.join("centrifuge_NT", sample + "_gDNA_centrifuge_NT_classification.summary.txt"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["gDNA_centrifuge_NT_unclassified"] = unclassified
            samples[sample]["gDNA_centrifuge_NT_bacterial"] = bacterial
            samples[sample]["gDNA_centrifuge_NT_eukaryotic"] = eukaryotic
            samples[sample]["gDNA_centrifuge_NT_top_orders"] = top_orders
    else:
            samples[sample]["gDNA_centrifuge_NT_unclassified"] = "NA"
            samples[sample]["gDNA_centrifuge_NT_bacterial"] = "NA"
            samples[sample]["gDNA_centrifuge_NT_eukaryotic"] = "NA"
            samples[sample]["gDNA_centrifuge_NT_top_orders"] = "NA"

    # Centrifuge NT cDNA
    if sample in cDNA_SAMPLES:
        with open(path.join("centrifuge_NT", sample + "_cDNA_centrifuge_NT_classification.summary.txt"), "r") as f:
            line = f.readline() # ignore first line header
            line = f.readline().strip().split("\t")

            unclassified = line[0]
            bacterial = line[1]
            eukaryotic = line[2]

            if len(line) == 4:
                top_orders = line[3]
            else:
                top_orders = "NA"

            samples[sample]["cDNA_centrifuge_NT_unclassified"] = unclassified
            samples[sample]["cDNA_centrifuge_NT_bacterial"] = bacterial
            samples[sample]["cDNA_centrifuge_NT_eukaryotic"] = eukaryotic
            samples[sample]["cDNA_centrifuge_NT_top_orders"] = top_orders
    else:
            samples[sample]["cDNA_centrifuge_NT_unclassified"] = "NA"
            samples[sample]["cDNA_centrifuge_NT_bacterial"] = "NA"
            samples[sample]["cDNA_centrifuge_NT_eukaryotic"] = "NA"
            samples[sample]["cDNA_centrifuge_NT_top_orders"] = "NA"

for sample in gDNA_SAMPLES:
    # rRNA genes from gDNA
    with open(path.join(sample, "rrna", sample + "_gDNA.rrna.blast.top.tsv"), "r") as f:
        hits = []
        for line in f:
            query, hit, identity, aln_length, *_ = line.strip().split("\t")

            # Only report hit as most detailed taxonomic level
            hit = hit.split("|")[-1]

            if "18S" in query:
                if int(aln_length) > int(snakemake.params.min_rRNA_blast_length):
                    hits.append(hit + " (" + identity + "% " + aln_length + " bp)")
        
        samples[sample]["gDNA_rRNA_top_hits"] = ", ".join(sorted(hits))

for sample in cDNA_SAMPLES:
    # rRNA genes from cDNA
    with open(path.join(sample, "rrna", sample + "_cDNA.rrna.blast.top.tsv"), "r") as f:
        hits = []
        for line in f:
            query, hit, identity, aln_length, *_ = line.strip().split("\t")

            # Only report hit as most detailed taxonomic level
            hit = hit.split("|")[-1]

            if "18S" in query:
                if int(aln_length) > int(snakemake.params.min_rRNA_blast_length):
                    hits.append(hit + " (" + identity + "% " + aln_length + " bp)")
        
        samples[sample]["cDNA_rRNA_top_hits"] = ", ".join(sorted(hits))

    # Diamond/blast transcript classification
    with open(path.join(sample, "trinity", sample + ".Trinity.cdhit.fasta.transdecoder.diamond.mapped.out"), "r") as f:
        hits = {}
        for line in f:
            print(line)
            classification = line.strip().split("\t")[1]

            if classification not in hits:
                hits[classification] = 1
            else:
                hits[classification] += 1

        hits = sorted(hits.items(), key = lambda x: x[1], reverse = True)

        formatted_hits = []
        for hit in hits:
            formatted_hits.append(hit[0] + " (" + str(hit[1]) + ")")

        samples[sample]["diamond_hits"] = ", ".join(formatted_hits)

fo = open(output_file, "w")

header = [""]
header += ["Raw gDNA reads", "Raw gDNA bases", "Trimmed gDNA reads", "Trimmed gDNA bases", "gDNA % bp retained", "gDNA duplication rate (fastp)"]
header += ["gDNA insert size peak (fastp)", "gDNA insert size median (bbmerge)"]
header += ["Raw cDNA reads", "Raw cDNA bases", "Trimmed cDNA reads", "Trimmed cDNA bases", "cDNA % bp retained", "cDNA duplication rate (fastp)"]
header += ["cDNA insert size peak (fastp)", "cDNA insert size peak (bbmerge)"]
header += ["gDNA Kraken2 % Unclassified", "gDNA Kraken2 % Bacterial", "gDNA Kraken2 % Eukaryotic", "gDNA Kraken2 Top Orders"]
header += ["cDNA Kraken2 % Unclassified", "cDNA Kraken2 % Bacterial", "cDNA Kraken2 % Eukaryotic", "cDNA Kraken2 Top Orders"]
header += ["gDNA Centrifuge % Unclassified", "gDNA Centrifuge % Bacterial", "gDNA Centrifuge % Eukaryotic", "gDNA Centrifuge Top Orders"]
header += ["cDNA Centrifuge % Unclassified", "cDNA Centrifuge % Bacterial", "cDNA Centrifuge % Eukaryotic", "cDNA Centrifuge Top Orders"]
header += ["gDNA Centrifuge NT % Unclassified", "gDNA Centrifuge NT % Bacterial", "gDNA Centrifuge NT % Eukaryotic", "gDNA Centrifuge NT Top Orders"]
header += ["cDNA Centrifuge NT % Unclassified", "cDNA Centrifuge NT % Bacterial", "cDNA Centrifuge NT % Eukaryotic", "cDNA Centrifuge NT Top Orders"]
header += ["gDNA Top 18S rRNA genes"]
header += ["cDNA Top 18S rRNA genes"]
header += ["Transcripts best hits"]

print("\t".join(header))
fo.write("\t".join(header) + "\n")

# TODO decide here what to include for samples where read count was lower than cutoff
for sample in ALL_SAMPLES:
    line = [sample]
    line.append(samples[sample]["gDNA_reads"])
    line.append(samples[sample]["gDNA_bases"])
    line.append(samples[sample]["gDNA_trimmed_reads"])
    line.append(samples[sample]["gDNA_trimmed_bases"])
    line.append(samples[sample]["gDNA_trimmed_bases_retained"])
    line.append(samples[sample]["gDNA_trimmed_duplication_rate"])
    line.append(samples[sample]["gDNA_insert_size_fastp"])
    line.append(samples[sample]["gDNA_insert_size_bbmerge"])    

    line.append(samples[sample]["cDNA_reads"])
    line.append(samples[sample]["cDNA_bases"])
    line.append(samples[sample]["cDNA_trimmed_reads"])
    line.append(samples[sample]["cDNA_trimmed_bases"])
    line.append(samples[sample]["cDNA_trimmed_bases_retained"])
    line.append(samples[sample]["cDNA_trimmed_duplication_rate"])
    line.append(samples[sample]["cDNA_insert_size_fastp"])
    line.append(samples[sample]["cDNA_insert_size_bbmerge"])

    line.append(samples[sample]["gDNA_kraken2_unclassified"])
    line.append(samples[sample]["gDNA_kraken2_bacterial"])
    line.append(samples[sample]["gDNA_kraken2_eukaryotic"])
    line.append(samples[sample]["gDNA_kraken2_top_orders"])

    line.append(samples[sample]["cDNA_kraken2_unclassified"])
    line.append(samples[sample]["cDNA_kraken2_bacterial"])
    line.append(samples[sample]["cDNA_kraken2_eukaryotic"])
    line.append(samples[sample]["cDNA_kraken2_top_orders"])

    line.append(samples[sample]["gDNA_centrifuge_unclassified"])
    line.append(samples[sample]["gDNA_centrifuge_bacterial"])
    line.append(samples[sample]["gDNA_centrifuge_eukaryotic"])
    line.append(samples[sample]["gDNA_centrifuge_top_orders"])

    line.append(samples[sample]["cDNA_centrifuge_unclassified"])
    line.append(samples[sample]["cDNA_centrifuge_bacterial"])
    line.append(samples[sample]["cDNA_centrifuge_eukaryotic"])
    line.append(samples[sample]["cDNA_centrifuge_top_orders"])

    line.append(samples[sample]["gDNA_centrifuge_NT_unclassified"])
    line.append(samples[sample]["gDNA_centrifuge_NT_bacterial"])
    line.append(samples[sample]["gDNA_centrifuge_NT_eukaryotic"])
    line.append(samples[sample]["gDNA_centrifuge_NT_top_orders"])

    line.append(samples[sample]["cDNA_centrifuge_NT_unclassified"])
    line.append(samples[sample]["cDNA_centrifuge_NT_bacterial"])
    line.append(samples[sample]["cDNA_centrifuge_NT_eukaryotic"])
    line.append(samples[sample]["cDNA_centrifuge_NT_top_orders"])

    if sample in gDNA_SAMPLES:
        if len(samples[sample]["gDNA_rRNA_top_hits"]) == 0:
            line.append("No genes or no  hits")
        else:    
            line.append(samples[sample]["gDNA_rRNA_top_hits"])
    else:
        line.append("NA")

    if sample in cDNA_SAMPLES:
        if len(samples[sample]["cDNA_rRNA_top_hits"]) == 0:
            line.append("No genes or no hits")
        else:
            line.append(samples[sample]["cDNA_rRNA_top_hits"])
        
        line.append(samples[sample]["diamond_hits"])
    else:
        line.append("NA")
        line.append("NA")

    print("\t".join(map(str, line)))
    fo.write("\t".join(map(str, line)) + "\n")

fo.close()
