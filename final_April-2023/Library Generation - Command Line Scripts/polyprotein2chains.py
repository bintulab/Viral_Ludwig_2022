# Name: polyprotein2chains.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

import os
import sys
import pandas as pd
import math

# sys.argv[1] = input csv file
# sys.argv[2] = uniprotkb or uniref
df = pd.read_csv(sys.argv[1])
filetype = sys.argv[2]

# for UniProtKB files
if filetype == 'uniprotkb':
	entryID = 'Entry'
	orgname = 'Organism'
	genename = 'Gene names'

# for UniRef files
elif filetype == 'uniref':
	entryID = 'Cluster ID'
	orgname = 'Organisms'
	genename = 'Cluster name'

else:
	sys.exit('Error: unspecified file type; please enter uniprotkb or uniref as second argument')

entryIDList = df[entryID].values.tolist()
chainList = df['Chain'].values.tolist()
orgnameList = df[orgname].values.tolist()
genenameList = df[genename].values.tolist()
sequenceList = df['Representative Sequence'].values.tolist()

chainIDList = []
entrynameList = []
chain_seqList = []
orgnameList2 = []
start = 0

for i in range(len(entryIDList)):
	seqlen = len(sequenceList[i])
	# if no chain info for polyprotein, store original info
	if pd.isna(chainList[i]):
		chainIDList.append(entryIDList[i])
		entrynameList.append(genenameList[i])
		orgnameList2.append(orgnameList[i])
		chain_seqList.append(sequenceList[i])
		continue
	# for polyproteins with chain info, split into constituent proteins
	chains = chainList[i].split('CHAIN ')[1:]
	for j in range(0, len(chains)):
		if ('polyprotein' in chains[j]) or ('Polyprotein' in chains[j]):
			continue
		chainList_elements = chains[j].split('; ')
		indices = chainList_elements[0].split('..')

		if indices[0][0] == '<':
			indices[0] = indices[0][1:]
		start_index = int(indices[0]) - 1
		if indices[1][0] == '>':
			indices[1] = indices[1][1:]
		end_index = int(indices[1])

		if (end_index - start_index) == seqlen:
			continue
		if end_index == seqlen:
			chain_seq = sequenceList[i][start_index:]
		else:
			chain_seq = sequenceList[i][start_index:end_index]

		if '' in chainList_elements:
			chainList_elements.remove('')
		entryname = chainList_elements[1].lstrip(' /note=')[1:-1]
		
		if len(chainList_elements) == 2:
			chainID = entryIDList[i] + '|' + 'none'
		elif len(chainList_elements) == 3:
			chainID = entryIDList[i] + '|' + chainList_elements[2].lstrip(' /id=')[1:-1]
		elif len(chainList_elements) == 4:
			chainID = entryIDList[i] + '|' + chainList_elements[3].lstrip(' /id=')[1:-1]
		else:
			print(chainList_elements)
		chainIDList.append(chainID)
		entrynameList.append(entryname)
		orgnameList2.append(orgnameList[i])
		chain_seqList.append(chain_seq)

# store chain info in dictionary/dataframe
chainDict = {'Chain ID':chainIDList,
			 'Chain Name':entrynameList,
			 'Virus':orgnameList2,
			 'Chain Sequence':chain_seqList}

chainDF = pd.DataFrame(chainDict)
chainDF.sort_values(by=['Chain ID'], inplace=True)
# chainDF.drop_duplicates('Chain Sequence', keep = 'first', inplace = True)

DFsavename = sys.argv[1][:-4] + '_chains.csv'
chainDF.to_csv(DFsavename, index=False)
