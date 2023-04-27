# Name: uniref_fasta2tsv.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

import os
import sys
import pandas as pd
import math

# sys.argv[1] = input csv file with codon-optimized oligos; must have 'Virus' column with various types/strains
df = pd.read_csv(sys.argv[1])

# define virus and subfamily classification scheme for mapping of tile to virus to subfamily
herpesvirusList = ['HHV-1', 'HHV-2', 'HHV-3', 'HHV-4', 'HHV-5', 'HHV-6', 'HHV-7', 'HHV-8', 'CeHV', 'SuHV']
herpesvirus2subfamily = {'HHV-1':'alpha_herpesvirus',
						 'HHV-2':'alpha_herpesvirus',
						 'HHV-3':'alpha_herpesvirus',
						 'CeHV':'alpha_herpesvirus',
						 'SuHV':'alpha_herpesvirus',
						 'HHV-4':'gamma_herpesvirus',
						 'HHV-5':'beta_herpesvirus',
						 'HHV-6':'beta_herpesvirus',
						 'HHV-7':'beta_herpesvirus',
						 'HHV-8':'gamma_herpesvirus',
						 'other':'alpha_herpesvirus'}

# initialize arrays to store mapping info
herpesvirus = []
subfamily = []

# for each row, count the number of times a herpesvirus type appears in the 'Virus' column
# store the most common as the associated virus type and store the corresponding subfamily info
for row in df['Virus']:
	v_count = 0
	v_id = 'other'
	for v in herpesvirusList:
		if row.count(v) > v_count:
			v_count = row.count(v)
			v_id = v
	herpesvirus.append(v_id)
	subfamily.append(herpesvirus2subfamily[v_id])


# add the virus type and subfamily (sublibrary) to the dataframe
df['Type'] = herpesvirus
df['Library'] = subfamily

# create new dataframes for each sublibrary
df_alpha = df[df['Library'].str.contains('alpha_herpesvirus')]
df_beta = df[df['Library'].str.contains('beta_herpesvirus')]
df_gamma = df[df['Library'].str.contains('gamma_herpesvirus')]

# save dataframes as csv files
savename_parent = sys.argv[1][:-4] + '_full-library.csv'
df.to_csv(savename_parent, index=False)

savename_alpha = sys.argv[1][:-4] + '_alpha-sublibrary.csv'
df_alpha.to_csv(savename_alpha, index=False)

savename_beta = sys.argv[1][:-4] + '_beta-sublibrary.csv'
df_beta.to_csv(savename_beta, index=False)

savename_gamma = sys.argv[1][:-4] + '_gamma-sublibrary.csv'
df_gamma.to_csv(savename_gamma, index=False)

print('Libray/sublibrary sizes before adding controls:')
print('Full herpesvirus library\t%d tiles' % len(df))
print('Alpha herpesvirus sublibrary\t%d tiles' % len(df_alpha))
print('Beta herpesvirus sublibrary\t%d tiles' % len(df_beta))
print('Gamma herpesvirus sublibrary\t%d tiles' % len(df_gamma))
