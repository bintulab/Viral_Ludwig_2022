# Name: domains_to_codon_opt_oligos.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

import os
import sys
import pandas as pd
import math
from dnachisel import *

# define optimizeOligo function:
#   |_ failed = list of oligos that fail codon optimization (input and output - updates)
#   |_ ID = tile ID (input)
#   |_ dna_sequence = starting DNA sequence to optimize (input; may be output if codon optimization fails)
#   |_ GCglobMax = upper limit for enforcement of GC content (input)
#   |_ optimizedDNA = codon-opimized DNA sequence (output)

def optimizeOligo(failed, ID, dna_sequence, GCglobMax):
	problem = DnaOptimizationProblem(
		sequence=dna_sequence,
		constraints=[
			AvoidPattern('BsmBI_site'),
			AvoidPattern('7xC'),
			AvoidRareCodons(species='h_sapiens', min_frequency=0.1),
		    EnforceGCContent(mini=0.25, maxi=GCglobMax),
		    EnforceGCContent(mini=0.20, maxi=0.70, window=50),
		    EnforceTranslation(),
		    AvoidStopCodons()],
		objectives=[CodonOptimize(species='h_sapiens', method='match_codon_usage')]
		)
	try:
		problem.resolve_constraints()
	except:
		if ID not in failed:
			failed.append(ID)
		GCglobMax += 0.01
		print('++++++++++ GOING UP TO %s ++++++++++' % str(GCglobMax))
		optimizedDNA, failed = optimizeOligo(failed, ID, dna_sequence, GCglobMax)
		return optimizedDNA, failed
	problem.optimize()
	optimizedDNA = str(problem.sequence)
	return optimizedDNA, failed


# sys.argv[1] = input csv with unique tiles
# sys.argv[2] = library name

df = pd.read_csv(sys.argv[1])
libName = sys.argv[2]

# store tile ID and protein sequence column values as lists
tileID = df['Tile ID'].values.tolist()
tileAAseq = df['Tile Sequence'].values.tolist()

# libName = [sys.argv[2]] * len(tileID)

# initialize arrays to store information from codon optimization
tileDNAseq = []
failedList = []

# for each tile protein sequence, reverse-translate and attempt to codon-optimize with a default max global GC content of 65%
# return the IDs of tiles that only pass codon optimization with relaxed global GC content constraints
for i in range(len(tileID)):
	initialDNAseq = reverse_translate(tileAAseq[i])
	coDNAseq, failedList = optimizeOligo(failedList, tileID[i], initialDNAseq, GCglobMax=0.65)
	tileDNAseq.append(coDNAseq)

# if any oligos required relaxed global GC content constraints for optimization, report these IDs and save a text file with them
if len(failedList) != 0:
	print('Oligos with global GC content > 65%: ', failedList)
	failedFile = sys.argv[1][:-4] + '_oligos-globalGC-gt65-IDs.txt'
	with open(failedFile, 'w') as f:
		for failedID in failedList:
			f.write("%s\n" % failedID)
else:
	print('All oligos passed codon optimization with specified constraints')

# add codon-optimized DNA sequence for each oligo and library tag to the dataframe and save
df['DNA Sequence'] = tileDNAseq
df['Library'] = libName
savename = sys.argv[1][:-4] + '_codon-opt-oligos.csv'

df.to_csv(savename, index=False)

