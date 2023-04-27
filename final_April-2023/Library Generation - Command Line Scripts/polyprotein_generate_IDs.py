# Name: polyprotein2chains.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/06/2020

import os
import sys
import pandas as pd

# sys.argv[1] = input csv file to filter
# sys.argv[2] = uniprotkb or uniref
df = pd.read_csv(sys.argv[1])
filetype = sys.argv[2]

# for UniProtKB files
# must have a 'Protein names' column for identifying polyproteins from UniProt!
if filetype == 'uniprotkb':
	entryID = 'Entry'
	entryname = 'Protein names'

# for UniRef files
elif filetype == 'uniref':
	entryID = 'Cluster ID'
	entryname = 'Cluster name'

else:
	sys.exit('Error: unspecified file type; please enter uniprotkb or uniref as second argument')

# create two dataframes, one with polyproteins only and one with non-polyproteins only
df_only_poly = df[df[entryname].str.lower().str.contains('polyprotein')==True]
df_wo_poly = df[df[entryname].str.lower().str.contains('polyprotein')==False]

# get IDs of polyproteins and save to text file for uploading to UniProt
repseqID_list = df_only_poly[entryID]
repseqID_list_saveFile = sys.argv[1][:-4] + '_polyprotein-ID-list.txt'

with open(repseqID_list_saveFile, 'w') as outFile:
	for i in repseqID_list:
		line = i + '\n'
		outFile.write(line)

# save dataframe with non-polyproteins only
filtered_df_saveFile = sys.argv[1][:-4] + '_wo-polyproteins.csv'
df_wo_poly.to_csv(filtered_df_saveFile, index=False)
