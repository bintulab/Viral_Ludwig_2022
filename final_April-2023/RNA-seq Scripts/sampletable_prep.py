#Made by ART, 1/10/22

import os
import sys
import pandas as pd
import csv 

# specify directory with SAM files after alignment
samDir = sys.argv[1]
parentDir = os.path.dirname(os.path.dirname(samDir))

names = []
for file in os.listdir(samDir):
	names.append(file.split('.')[0])

header = ['', 'SampleName', 'cell', 'plasmid', 'dox', 'Rep', 'Run']
cells = []
plasmids = []
dox = []
Rep = []

for name in names:
	name_components = name.split('_')
	cells.append(name_components[0])
	plasmids.append(name_components[1])
	dox.append(name_components[2])
	Rep.append(name_components[3])

files = [names, names, cells, plasmids, dox, Rep, names]
prep = zip(*files)

sample_table = os.path.join(parentDir, 'sample_table.csv')
with open(sample_table, 'w', encoding = 'UTF8') as f:
	writer = csv.writer(f)
	writer.writerow(header)
	writer.writerows(prep)
