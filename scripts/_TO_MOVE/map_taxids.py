#!/usr/bin/env python

# Maps taxonomy IDs to species names and lineages
# Usage: python map_taxids.py fullnamelineage.dmp query_ids.txt
# Prints to stdout: classification\tfull lineage

import sys

if len(sys.argv) < 2:
	print("Usage: python map_taxids.py fullnamelineage.dmp query_ids.txt")
	sys.exit(0)

mapping = {}

with open(sys.argv[1], "r") as f:
	for line in f:
		l = line.strip().split("|")

		mapping[l[0].strip()] = [l[1].strip(), l[2].strip()]

with open(sys.argv[2], "r") as f:
	for line in f:
		l = line.strip().split("\t")
		tax_id = l[1]
		
		if tax_id == "0":
			print("Unclassified")
		else:
			try:
				print("\t".join(mapping[tax_id]))
			except:
				print("MISSING ID FROM DUMPFILE ->", tax_id)
