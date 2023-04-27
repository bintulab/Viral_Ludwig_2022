# Name: generate_tiles_v2.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

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
tile_len = int(sys.argv[2])
intertile_dist = int(sys.argv[3])
filetype = sys.argv[4]


# for UniProtKB files
if filetype == 'uniprotkb':
	entryIDList = df['Entry'].values.tolist()
	entrynameList = df['Entry name'].values.tolist()
	orgnameList = df['Organism'].values.tolist()
	sequenceList = df['Representative Sequence'].values.tolist()

# for UniRef files
elif filetype == 'uniref':
	entryIDList = df['Cluster ID'].values.tolist()
	entrynameList = df['Cluster name'].values.tolist()
	orgnameList = df['Organisms'].values.tolist()
	sequenceList = df['Representative Sequence'].values.tolist()

# for polyprotein files:
elif filetype == 'polyprotein':
	entryIDList = df['Chain ID'].values.tolist()
	entrynameList = df['Chain Name'].values.tolist()
	orgnameList = df['Virus'].values.tolist()
	sequenceList = df['Chain Sequence'].values.tolist()

else:
	sys.exit('Error: unspecified file type; please enter uniprotkb, uniref, or polyprotein as fourth argument')


tileIDList = []
new_entrynameList = []
startList = []
endList = []
tile_seqList = []
orgnameList2 = []
start = 0

# note: tiles are 1-indexed
for i in range(len(entryIDList)):
	seqlen = len(sequenceList[i])
	# if sequence is shorter than tile length
	if seqlen <= tile_len:
		ntiles = 1
		tileID = entryIDList[i] + '_001'
		tileIDList.append(tileID)
		new_entrynameList.append(entrynameList[i])
		startList.append(1)
		endList.append(seqlen)
		tile_seqList.append(sequenceList[i])
		orgnameList2.append(orgnameList[i])

	# if sequence is longer than tile length
	else:
		ntiles = math.ceil((seqlen - tile_len)/intertile_dist + 1)
		for j in range(ntiles - 1):
			start = j * intertile_dist
			end = start + tile_len
			tileID = entryIDList[i] + '_' + str(j+1).zfill(3)
			tileIDList.append(tileID)
			new_entrynameList.append(entrynameList[i])
			# use 1-indexing for listing amino acids in tile (more intuitive)
			startList.append(start + 1)
			endList.append(end)
			tile_seqList.append(sequenceList[i][start:end])
			orgnameList2.append(orgnameList[i])
		
		# last tile
		tileID = entryIDList[i] + '_' + str(ntiles).zfill(3)
		tileIDList.append(tileID)
		new_entrynameList.append(entrynameList[i])
		startList.append(seqlen - tile_len + 1)
		endList.append(seqlen)
		tile_seqList.append(sequenceList[i][-tile_len:])
		orgnameList2.append(orgnameList[i])

# organize info into dictionary and make dataframe from dictionary
tileDict = {'Tile ID':tileIDList,
			'Tile Info':new_entrynameList,
			'Tile Start':startList,
			'Tile End':endList,
			'Virus':orgnameList2,
			'Tile Sequence':tile_seqList}

tileDF = pd.DataFrame(tileDict)
tileDF.sort_values(by=['Tile ID'], inplace=True)
len_prededup = len(tileDF)

# save all tile info
DFsavename_alltiles = sys.argv[1][:-4] + '_all-tiles.csv'
tileDF.to_csv(DFsavename_alltiles, index=False)


# summarize all tile info
summaryDF = pd.DataFrame()
summaryDF['Entry ID'] = tileDF['Tile ID'].map(lambda x: '_'.join(x.split('_')[:-1]))
summaryDF['Entry Descriptor'] = tileDF['Tile Info']
summaryDF['Virus'] = tileDF['Virus']
summaryDF['Total Number of Tiles'] = 1
summaryDF = summaryDF.groupby(by=['Entry ID', 'Entry Descriptor', 'Virus'], as_index=False)['Total Number of Tiles'].sum()


# de-duplicate tiles and save info
tileDF.drop_duplicates('Tile Sequence', keep = 'first', inplace = True)
len_postdedup = len(tileDF)
DFsavename_uniquetiles = sys.argv[1][:-4] + '_unique-tiles.csv'
tileDF.to_csv(DFsavename_uniquetiles, index=False)


# plot unique tile length distribution
tileDF['Tile Length'] = tileDF['Tile Sequence'].map(lambda x: len(x))
g = sns.distplot(tileDF['Tile Length'], bins=80, kde=False)
plt.ylabel('Count (Log10)')
plt.yscale('log')
plt.xlabel('Tile Length Distribution (Amino Acids)')
plt.xlim(20,70)

figsaveFile = sys.argv[1][:-4] + '_unique-tiles_distr.png'
fig_to_save = g.get_figure()
fig_to_save.savefig(figsaveFile)


# summarize after de-duplication
summaryDF_dedup = pd.DataFrame()
summaryDF_dedup['Entry ID'] = tileDF['Tile ID'].map(lambda x: '_'.join(x.split('_')[:-1]))
summaryDF_dedup['Tiles after De-Duplication'] = 1
summaryDF_dedup = summaryDF_dedup.groupby(by=['Entry ID'], as_index=False)['Tiles after De-Duplication'].sum()


# creat final summary
final_summaryDF = pd.merge(summaryDF, summaryDF_dedup, on='Entry ID', how='outer').fillna(0)
summDFsavename = sys.argv[1][:-4] + '_tile-summary.csv'
final_summaryDF.to_csv(summDFsavename, index=False)

print(final_summaryDF)
print('\nTotal Tile Number before De-Duplication: %d' % len_prededup)
print('\nTotal Tile Number after De-Duplication: %d' % len_postdedup)


