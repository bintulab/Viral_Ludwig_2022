# Name: qc_oligos_GC_content.py
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
tileID = df['Tile ID'].values.tolist()
libName = df['Library'][1]

GCcontent = []

for n in range(len(DNAseqList)):
	GCcount = DNAseqList[n].count('G') + DNAseqList[n].count('C')
	seqLen = len(DNAseqList[n])
	GCcontent.append(100 * GCcount/seqLen)
	if GCcontent[n] > 65:
		print(tileID[n], GCcontent[n])

gte65 = str(round(float(100 * sum(i > 60 for i in GCcontent)/len(GCcontent)), 2))
lte35 = str(round(float(100 * sum(j < 30 for j in GCcontent)/len(GCcontent)), 2))
percent = '%'
print('%s%s of tiles have GC content greater than 60%s' % (gte65, percent, percent))
print('%s%s of tiles have GC content less than 30%s' % (lte35, percent, percent))

g = sns.distplot(GCcontent, bins=100, kde=True)
plt.axvline(x=25, linestyle='--', color='#555555')
plt.axvline(x=65, linestyle='--', color='#555555')
plt.axvline(x=30, linestyle='--', color='#DDDDDD')
plt.axvline(x=60, linestyle='--', color='#DDDDDD')
plt.ylabel('Count')
# plt.yscale('log')
plt.xlabel('GC Content (Percent)')
plt.xlim(20,70)
plt.title('Oligo GC Content Distribution for %s Library' % libName)

lte35_plot = '%s%s' % (lte35, percent)
gte65_plot = '%s%s' % (gte65, percent)

height = 0.1

plt.text(25.5, height, lte35_plot)
plt.text(60.5, height, gte65_plot)

figsaveFile = sys.argv[1][:-4] + '_GC-content_distr.png'
fig_to_save = g.get_figure()
fig_to_save.savefig(figsaveFile)
