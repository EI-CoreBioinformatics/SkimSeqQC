#!/usr/bin/env python

# reports top BLAST hits from BLAST tabular output

import sys

seen = []

with open(sys.argv[1], "r") as f:
	for line in f:
		query = line.strip().split("\t")[0]
		
		if query not in seen:
			seen.append(query)
			print(line.strip())
