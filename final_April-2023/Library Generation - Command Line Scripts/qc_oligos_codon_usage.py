# Name: qc_oligos_codon_usage.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 9/29/2020

import os
import sys
import pandas as pd
import math
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
# from Bio import Entrez, SeqIO
import seaborn as sns

df = pd.read_csv(sys.argv[1])
DNAseqList = df['DNA Sequence'].values.tolist()
libName = df['Library'][1]

# Codon Usage Table from GenScript: https://www.genscript.com/tools/codon-frequency-table
humanUsage = {'TTT':['F', 0.45, 16.9], 'TTC':['F', 0.55, 20.4], 'TTA':['L', 0.07,  7.2], 'TTG':['L', 0.13, 12.6],
			  'TAT':['Y', 0.43, 12.0], 'TAC':['Y', 0.57, 15.6], 'TAA':['*', 0.28,  0.7], 'TAG':['*', 0.20,  0.5],
			  'CTT':['L', 0.13, 12.8], 'CTC':['L', 0.20, 19.4], 'CTA':['L', 0.07,  6.9], 'CTG':['L', 0.41, 40.3],
			  'CAT':['H', 0.41, 10.4], 'CAC':['H', 0.59, 14.9], 'CAA':['Q', 0.25, 11.8], 'CAG':['Q', 0.75, 34.6],
			  'ATT':['I', 0.36, 15.7], 'ATC':['I', 0.48, 21.4], 'ATA':['I', 0.16,  7.1], 'ATG':['M', 1.00, 22.3],
			  'AAT':['N', 0.46, 16.7], 'AAC':['N', 0.54, 19.5], 'AAA':['K', 0.42, 24.0], 'AAG':['K', 0.58, 32.9],
			  'GTT':['V', 0.18, 10.9], 'GTC':['V', 0.24, 14.6], 'GTA':['V', 0.11,  7.0], 'GTG':['V', 0.47, 28.9],
			  'GAT':['D', 0.46, 22.3], 'GAC':['D', 0.54, 26.0], 'GAA':['E', 0.42, 29.0], 'GAG':['E', 0.58, 40.8],
			  'TCT':['S', 0.18, 14.6], 'TCC':['S', 0.22, 17.4], 'TCA':['S', 0.15, 11.7], 'TCG':['S', 0.06,  4.5],
			  'TGT':['C', 0.45,  9.9], 'TGC':['C', 0.55, 12.2], 'TGA':['*', 0.52,  1.3], 'TGG':['W', 1.00, 12.8],
			  'CCT':['P', 0.28, 17.3], 'CCC':['P', 0.33, 20.0], 'CCA':['P', 0.27, 16.7], 'CCG':['P', 0.11,  7.0],
			  'CGT':['R', 0.08,  4.7], 'CGC':['R', 0.19, 10.9], 'CGA':['R', 0.11,  6.3], 'CGG':['R', 0.21, 11.9],
			  'ACT':['T', 0.24, 12.8], 'ACC':['T', 0.36, 19.2], 'ACA':['T', 0.28, 14.8], 'ACG':['T', 0.12,  6.2],
			  'AGT':['S', 0.15, 11.9], 'AGC':['S', 0.24, 19.4], 'AGA':['R', 0.20, 11.5], 'AGG':['R', 0.20, 11.4],
			  'GCT':['A', 0.26, 18.6], 'GCC':['A', 0.40, 28.5], 'GCA':['A', 0.23, 16.0], 'GCG':['A', 0.11,  7.6],
			  'GGT':['G', 0.16, 10.8], 'GGC':['G', 0.34, 22.8], 'GGA':['G', 0.25, 16.3], 'GGG':['G', 0.25, 16.4]}

oligoUsage = {'TTT':0, 'TTC':0, 'TTA':0, 'TTG':0,
			  'TAT':0, 'TAC':0, 'TAA':0, 'TAG':0,
			  'CTT':0, 'CTC':0, 'CTA':0, 'CTG':0,
			  'CAT':0, 'CAC':0, 'CAA':0, 'CAG':0,
			  'ATT':0, 'ATC':0, 'ATA':0, 'ATG':0,
			  'AAT':0, 'AAC':0, 'AAA':0, 'AAG':0,
			  'GTT':0, 'GTC':0, 'GTA':0, 'GTG':0,
			  'GAT':0, 'GAC':0, 'GAA':0, 'GAG':0,
			  'TCT':0, 'TCC':0, 'TCA':0, 'TCG':0,
			  'TGT':0, 'TGC':0, 'TGA':0, 'TGG':0,
			  'CCT':0, 'CCC':0, 'CCA':0, 'CCG':0,
			  'CGT':0, 'CGC':0, 'CGA':0, 'CGG':0,
			  'ACT':0, 'ACC':0, 'ACA':0, 'ACG':0,
			  'AGT':0, 'AGC':0, 'AGA':0, 'AGG':0,
			  'GCT':0, 'GCC':0, 'GCA':0, 'GCG':0,
			  'GGT':0, 'GGC':0, 'GGA':0, 'GGG':0}


for N in range(len(DNAseqList)):
	codons = [DNAseqList[N][k:k+3] for k in range(0, len(DNAseqList[N]), 3)]
	for c in codons:
		oligoUsage[c] += 1

usageDF = pd.DataFrame(humanUsage)
usageDF = usageDF.T
usageDF.columns = ['Residue', 'Fraction', 'Frequency/1000']
usageDF['Calculation'] = 'Expected'

oligoDF = pd.DataFrame(oligoUsage, index=[0])
oligoDF = oligoDF.T
oligoDF.columns = ['Frequency/1000']
total_codons = oligoDF['Frequency/1000'].sum()
oligoDF['Frequency/1000'] = 1000 * oligoDF['Frequency/1000']/total_codons
oligoDF['Residue'] = usageDF['Residue']
oligoDF['Fraction'] = usageDF['Fraction']
oligoDF['Calculation'] = 'Designed'

combinedDF = pd.concat([usageDF, oligoDF], ignore_index=False)
combinedDF['Codon'] = combinedDF.index
combinedDF['Codon & Residue'] = combinedDF['Codon'] + ' (' + combinedDF['Residue'] + ')'

# combinedDF = pd.concat([usageDF, oligoDF], ignore_index=True)
# combinedDF = combinedDF.T
# combinedDF.columns = ['Residue', 'Fraction', 'Expected', 'Designed']
# combinedDF['Designed'] = combinedDF['Designed']/1000

print(combinedDF)

plt.figure(figsize=(20,5))
g = sns.barplot(x='Codon & Residue', y='Frequency/1000', hue='Calculation', data=combinedDF, palette="Paired")
plt.ylabel('Frequency/1000')
plt.xlabel('Codon')
plt.ylim(0,50)
plt.xticks(rotation=90)
plt.title('Codon Usage for %s Library' % libName)
plt.tight_layout()


figsaveFile = sys.argv[1][:-4] + '_codon-usage.png'
fig_to_save = g.get_figure()
fig_to_save.savefig(figsaveFile)



# rows = 8
# cols = 8
# fig, axs = plt.subplots(rows, cols, figsize=(16,16), dpi=300)
# counter = 0
# for r in np.arange(0, rows):
# 	for c in np.arange(0, cols):
# 		sns.b

# oligoFract = {'TTT':(oligoUsage['TTT']/(oligoUsage['TTT'] + oligoUsage['TTC'])),
# 			  'TTC':(oligoUsage['TTC']/(oligoUsage['TTT'] + oligoUsage['TTC'])),
# 			  'TTA':(oligoUsage['TTA']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'TTG':(oligoUsage['TTG']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'TAT':(oligoUsage['TAT']/(oligoUsage['TAT'] + oligoUsage['TAC'])),
# 			  'TAC':(oligoUsage['TAC']/(oligoUsage['TAT'] + oligoUsage['TAC'])),
# 			  'CTT':(oligoUsage['CTT']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'CTC':(oligoUsage['CTC']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'CTA':(oligoUsage['CTA']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'CTG':(oligoUsage['CTG']/(oligoUsage['TTA'] + oligoUsage['TTG'] + oligoUsage['CTT'] + oligoUsage['CTC'] + oligoUsage['CTA'] + oligoUsage['CTG'])),
# 			  'CAT':(oligoUsage['CAT']/(oligoUsage['CAT'] + oligoUsage['CAC'])),
# 			  'CAC':(oligoUsage['CAC']/(oligoUsage['CAT'] + oligoUsage['CAC'])),
# 			  'CAA':(oligoUsage['CAA']/(oligoUsage['CAA'] + oligoUsage['CAG'])),
# 			  'CAG':(oligoUsage['CAG']/(oligoUsage['CAA'] + oligoUsage['CAG'])),
# 			  'ATT':(oligoUsage['ATT']/(oligoUsage['ATT'] + oligoUsage['ATC'] + oligoUsage['ATA'])),
# 			  'ATC':(oligoUsage['ATC']/(oligoUsage['ATT'] + oligoUsage['ATC'] + oligoUsage['ATA'])),
# 			  'ATA':(oligoUsage['ATA']/(oligoUsage['ATT'] + oligoUsage['ATC'] + oligoUsage['ATA'])),
# 			  'AAT':(oligoUsage['AAT']/(oligoUsage['AAT'] + oligoUsage['AAC'])),
# 			  'AAC':(oligoUsage['AAC']/(oligoUsage['AAT'] + oligoUsage['AAC'])),
# 			  'AAA':(oligoUsage['AAA']/(oligoUsage['AAA'] + oligoUsage['AAG'])),
# 			  'AAG':(oligoUsage['AAG']/(oligoUsage['AAA'] + oligoUsage['AAG'])),
# 			  'GTT':(oligoUsage['GTT']/(oligoUsage['GTT'] + oligoUsage['GTC'] + oligoUsage['GTA'] + oligoUsage['GTG'])),
# 			  'GTC':(oligoUsage['GTC']/(oligoUsage['GTT'] + oligoUsage['GTC'] + oligoUsage['GTA'] + oligoUsage['GTG'])),
# 			  'GTA':(oligoUsage['GTA']/(oligoUsage['GTT'] + oligoUsage['GTC'] + oligoUsage['GTA'] + oligoUsage['GTG'])),
# 			  'GTG':(oligoUsage['GTG']/(oligoUsage['GTT'] + oligoUsage['GTC'] + oligoUsage['GTA'] + oligoUsage['GTG'])),
# 			  'GAT':(oligoUsage['GAT']/(oligoUsage['GAT'] + oligoUsage['GAC'])),
# 			  'GAC':(oligoUsage['GAC']/(oligoUsage['GAT'] + oligoUsage['GAC'])),
# 			  'GAA':(oligoUsage['GAA']/(oligoUsage['GAA'] + oligoUsage['GAG'])),
# 			  'GAG':(oligoUsage['GAG']/(oligoUsage['GAA'] + oligoUsage['GAG'])),
# 			  'TCT':(oligoUsage['TCT']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'TCC':(oligoUsage['TCC']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'TCA':(oligoUsage['TCA']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'TCG':(oligoUsage['TCG']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'TGT':(oligoUsage['TGT']/(oligoUsage['TGT'] + oligoUsage['TGC'])),
# 			  'TGC':(oligoUsage['TGC']/(oligoUsage['TGT'] + oligoUsage['TGC'])),
# 			  'CCT':(oligoUsage['CCT']/(oligoUsage['CCT'] + oligoUsage['CCC'] + oligoUsage['CCA'] + oligoUsage['CCG'])),
# 			  'CCC':(oligoUsage['CCC']/(oligoUsage['CCT'] + oligoUsage['CCC'] + oligoUsage['CCA'] + oligoUsage['CCG'])),
# 			  'CCA':(oligoUsage['CCA']/(oligoUsage['CCT'] + oligoUsage['CCC'] + oligoUsage['CCA'] + oligoUsage['CCG'])),
# 			  'CCG':(oligoUsage['CCG']/(oligoUsage['CCT'] + oligoUsage['CCC'] + oligoUsage['CCA'] + oligoUsage['CCG'])),
# 			  'CGT':(oligoUsage['CGT']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'CGC':(oligoUsage['CGC']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'CGA':(oligoUsage['CGA']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'CGG':(oligoUsage['CGG']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'ACT':(oligoUsage['ACT']/(oligoUsage['ACT'] + oligoUsage['ACC'] + oligoUsage['ACA'] + oligoUsage['ACG'])),
# 			  'ACC':(oligoUsage['ACC']/(oligoUsage['ACT'] + oligoUsage['ACC'] + oligoUsage['ACA'] + oligoUsage['ACG'])),
# 			  'ACA':(oligoUsage['ACA']/(oligoUsage['ACT'] + oligoUsage['ACC'] + oligoUsage['ACA'] + oligoUsage['ACG'])),
# 			  'ACG':(oligoUsage['ACG']/(oligoUsage['ACT'] + oligoUsage['ACC'] + oligoUsage['ACA'] + oligoUsage['ACG'])),
# 			  'AGT':(oligoUsage['AGT']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'AGC':(oligoUsage['AGC']/(oligoUsage['TCT'] + oligoUsage['TCC'] + oligoUsage['TCA'] + oligoUsage['TCG'] + oligoUsage['AGT'] + oligoUsage['AGC'])),
# 			  'AGA':(oligoUsage['AGA']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'AGG':(oligoUsage['AGG']/(oligoUsage['CGT'] + oligoUsage['CGC'] + oligoUsage['CGA'] + oligoUsage['CGG'] + oligoUsage['AGA'] + oligoUsage['AGG'])),
# 			  'GCT':(oligoUsage['GCT']/(oligoUsage['GCT'] + oligoUsage['GCC'] + oligoUsage['GCA'] + oligoUsage['GCG'])),
# 			  'GCC':(oligoUsage['GCC']/(oligoUsage['GCT'] + oligoUsage['GCC'] + oligoUsage['GCA'] + oligoUsage['GCG'])),
# 			  'GCA':(oligoUsage['GCA']/(oligoUsage['GCT'] + oligoUsage['GCC'] + oligoUsage['GCA'] + oligoUsage['GCG'])),
# 			  'GCG':(oligoUsage['GCG']/(oligoUsage['GCT'] + oligoUsage['GCC'] + oligoUsage['GCA'] + oligoUsage['GCG'])),
# 			  'GGT':(oligoUsage['GGT']/(oligoUsage['GGT'] + oligoUsage['GGC'] + oligoUsage['GGA'] + oligoUsage['GGG'])),
# 			  'GGC':(oligoUsage['GGC']/(oligoUsage['GGT'] + oligoUsage['GGC'] + oligoUsage['GGA'] + oligoUsage['GGG'])),
# 			  'GGA':(oligoUsage['GGA']/(oligoUsage['GGT'] + oligoUsage['GGC'] + oligoUsage['GGA'] + oligoUsage['GGG'])),
# 			  'GGG':(oligoUsage['GGG']/(oligoUsage['GGT'] + oligoUsage['GGC'] + oligoUsage['GGA'] + oligoUsage['GGG']))}



# print(oligoUsage)
# print(oligoFract)
# print(oligoCount)



# oligoUsage = {'TTT':[], 'TTC':[], 'TTA':[], 'TTG':[],
# 			  'TAT':[], 'TAC':[], 'TAA':[], 'TAG':[],
# 			  'CTT':[], 'CTC':[], 'CTA':[], 'CTG':[],
# 			  'CAT':[], 'CAC':[], 'CAA':[], 'CAG':[],
# 			  'ATT':[], 'ATC':[], 'ATA':[], 'ATG':[],
# 			  'AAT':[], 'AAC':[], 'AAA':[], 'AAG':[],
# 			  'GTT':[], 'GTC':[], 'GTA':[], 'GTG':[],
# 			  'GAT':[], 'GAC':[], 'GAA':[], 'GAG':[],
# 			  'TCT':[], 'TCC':[], 'TCA':[], 'TCG':[],
# 			  'TGT':[], 'TGC':[], 'TGA':[], 'TGG':[],
# 			  'CCT':[], 'CCC':[], 'CCA':[], 'CCG':[],
# 			  'CGT':[], 'CGC':[], 'CGA':[], 'CGG':[],
# 			  'ACT':[], 'ACC':[], 'ACA':[], 'ACG':[],
# 			  'AGT':[], 'AGC':[], 'AGA':[], 'AGG':[],
# 			  'GCT':[], 'GCC':[], 'GCA':[], 'GCG':[],
# 			  'GGT':[], 'GGC':[], 'GGA':[], 'GGG':[]}
