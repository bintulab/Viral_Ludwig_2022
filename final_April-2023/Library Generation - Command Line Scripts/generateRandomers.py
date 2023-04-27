'''
Use DNA chisel to generate random negative control sequences

09-29-2020
'''

import pandas as pd
import re
import argparse
from dnachisel import *

parser = argparse.ArgumentParser(description='Make sequences for random negative controls')
# parser.add_argument('residues', type=int,
#                     help='Number of residues per element')

parser.add_argument('number', type=int,
                    help='Number of random control sequences to include',
                    default = 500)
parser.add_argument('length', type=int,
                    help='Length of oligo',
                    default = 240)
parser.add_argument('out_file', type=str,
                    help='Output filename start')

#########################################################################################################

# Edited by Connor 10/05/2020 to match domains_to_codon_opt_oligos_v2.py
def optimizeOligo(dna_sequence, pattern):
	problem = DnaOptimizationProblem(
		sequence=dna_sequence,
		constraints=[AvoidPattern('BsmBI_site'),
					 AvoidPattern('7xC'),
					 AvoidRareCodons(species='h_sapiens', min_frequency=0.1),
					 # newly added
					 EnforceGCContent(mini=0.25, maxi=0.65),
					 EnforceGCContent(mini=0.20, maxi=0.70, window=50),
					 AvoidStopCodons(),
					 # original enforcement in Nicole's script (toggled off here)
					 # EnforceGCContent(mini=0.20, maxi=0.75, window=50),
					 EnforceTranslation()],
		objectives=[CodonOptimize(species='h_sapiens', method='match_codon_usage')]
	)
	try:
		problem.resolve_constraints()
	except (TypeError, KeyError):
		print('optimization failed')
		return dna_sequence
	problem.optimize()
	optDNA = str(problem.sequence)
	return optDNA

#########################################################################################################
args = parser.parse_args()
BsmBIpattern = EnzymeSitePattern("BsmBI")
print('Avoiding BsmBI sites')
numControls = args.number
lenControls = args.length

rows_list = [] #list of sequences

# Make random negative controls
print('Generating '+str(numControls)+' random controls')
for ctrlNumber in range(0,numControls):
	# get a random protein with no stop codons
	dna_sequence = random_dna_sequence(lenControls)
	while '*' in translate(dna_sequence):
		dna_sequence = random_dna_sequence(lenControls)
	# optimize to remove RE sites	
	variantDNA = optimizeOligo(dna_sequence, BsmBIpattern)
	variant = translate(variantDNA)
	# original
	# rows_list.append({'label':'Random_control;'+str(ctrlNumber), 'Variant protein': variant, 'Variant DNA': variantDNA, 'Category': 'Random control'})
	# newly added 10/05/202
	rows_list.append({'Tile ID':'random_' + str(ctrlNumber), 'Tile Sequence': variant, 'DNA Sequence': variantDNA, 'Library': 'random_control'})


# Make df
df = pd.DataFrame(rows_list, columns = ['Tile ID', 'Tile Sequence','DNA Sequence', 'Library'])

print('Saving randomers. # Oligos:', df.shape[0])
df.to_csv(args.out_file, index = False)
# df.to_csv(args.out_file + '.csv', index = False)

