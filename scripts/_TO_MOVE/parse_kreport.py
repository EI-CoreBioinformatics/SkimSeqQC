#!/usr/bin/env python

# Parses Kraken2/Centrifuge kreport to report % Unclassified, % Bacteria and % Eukaryotic, and top 5 orders

import sys
from os.path import exists

if not exists(sys.argv[1]):
	print("File missing: " +  sys.argv[1])
	sys.exit(0)

bacterial = "0"
eukaryotic = "0"
unclassified = "0"

orders = []

with open(sys.argv[1], "r") as f:
	for line in f:
		l = line.strip().split("\t")
		for i in range(0, len(l)):
			l[i] = l[i].strip()

		if l[3] == "U" and l[5] == "unclassified":
			unclassified = l[0]
		elif l[3] == "D" and l[5] == "Bacteria":
			bacterial = l[0]
		elif l[3] == "D" and l[5] == "Eukaryota":
			eukaryotic = l[0]
		elif l[3] == "O":
			order = [l[5], float(l[0])]
			orders.append(order)

sorted_orders = sorted(orders, key = lambda x: x[1], reverse = True)

top_5_orders =  sorted_orders[0:5]
top_5_orders_formatted = []
for i in top_5_orders:
	top_5_orders_formatted.append(i[0] + " (" + str(i[1]) + "%)")

print("#" + "\t".join(["%Unclassified", "%Bacterial", "%Eukaryotic", "Top_Orders"]))
print("\t".join([unclassified, bacterial, eukaryotic, ";".join(top_5_orders_formatted)]))
