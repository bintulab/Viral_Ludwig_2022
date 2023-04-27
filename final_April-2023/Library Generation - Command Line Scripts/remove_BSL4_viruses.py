# remove_BSL4_viruses.py
import os
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])

BSL4_filter = []
with open(sys.argv[2], 'r') as virusinFile:
	for line in virusinFile:
		BSL4_filter.append(line.rstrip('\n').lower())
		# BSL4_filter.append(line.rstrip('\n').lower())

# need to change to 'Organisms' for UniRef
df['Organism_lowercase'] = df['Organism'].str.lower()

for keyword in BSL4_filter:
	df = df[df['Organism_lowercase'].str.contains(keyword)==False]

df.drop(columns=['Organism_lowercase'], inplace=True)

print(df)

savename = sys.argv[1][:-4] + '_BSL4-filtered.csv'
df.to_csv(savename, index=False)
