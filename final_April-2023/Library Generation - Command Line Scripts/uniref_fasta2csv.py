# Name: uniref_fasta2csv.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

import os
import sys
import pandas as pd

# sys.argv[1] = UniRef90 fasta file with representative protein sequences for human viral proteins
# sys.argv[2] = UniRef90 metadata tab file for human viral protein clusters
# sys.argv[3] = desired save name for csv file with compiled data

uniref_fasta = sys.argv[1]
uniref_metadata = sys.argv[2]
savename = sys.argv[3]

seqIDList = []
fastaList = []
new_seqID = ''
seqID = ''
seq = ''
first = True

# extract representative fasta sequence for each clusters
with open(uniref_fasta, 'r') as fastaFile:
	for line in fastaFile:
		if line[0] == '>':
			seqID = new_seqID
			new_seqID = line.rstrip('\n').split(' ')[0][1:]

			if first:
				first = False
				continue

			seqIDList.append(seqID)
			fastaList.append(seq)
			seq = ''
			continue
			
		else:
			seq += line.rstrip('\n')

seqIDList.append(new_seqID)
fastaList.append(seq)

# read metadata into dataframe
df = pd.read_csv(uniref_metadata, sep='\t')

# drop columns that are not useful
df.drop(columns=['Unnamed: 1', 'Identity'], inplace=True)

# put fasta sequences from above into dataframe
repseqDict = {'Cluster ID':seqIDList,
			  'Representative Sequence':fastaList}
df2 = pd.DataFrame(repseqDict)

# sort both dataframes by cluster ID in case they are not in the same order
df.sort_values(by=['Cluster ID'], inplace=True)
df2.sort_values(by=['Cluster ID'], inplace=True)

# add new column to the metadata dataframe that contains the representative sequence
df['Representative Sequence'] = df2['Representative Sequence']

# sort values by original index
df.sort_index(inplace=True)
df['Cluster ID'] = df['Cluster ID'].map(lambda x: x.lstrip('UniRef90_'))
df['Cluster name'] = df['Cluster name'].map(lambda x: x.lstrip('Cluster: '))

# check if any sequence contains 'X' and store this info in a different dataframe for manual inspection/correction
df_noX = df[df['Representative Sequence'].str.contains('X') == False]
df_X = df[df['Representative Sequence'].str.contains('X') == True]

# save dataframes as csv files
saveDir = os.path.dirname(uniref_fasta)
savepath = os.path.join(saveDir, savename)

# if proteins with at least one X were detected
if len(df_X) != 0:
	savepath_X = savepath[:-4] + '_contain-X-to-correct.csv'
	df_X.to_csv(savepath_X, index=False)

	savepath_noX = savepath[:-4] + '_no-X.csv'
	df_noX.to_csv(savepath_noX, index=False)

	print('Detected %d protein sequences containing at least one X' % len(df_X))
	print('Creating a separate CSV file for these proteins - please inspect/correct manually ')

# if no Xs detected
else:
	df_noX.to_csv(savepath, index=False)
	print('Data compiled successfully - no Xs detected in protein sequences')
	
