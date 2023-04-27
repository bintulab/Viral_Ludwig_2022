# Name: compile_dataframes.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Updated: 10/03/2020

import os
import sys
import pandas as pd
import math

# sys.argv[1] = full path and name of csv outfile with compiled dataframe
# sys.argv[2 through N] = full path and name of csv infiles with dataframes that you wish to compile
# NOTE: this script is intended to concatenate dataframes with the same structure (column names)

savepath = sys.argv[1]

# initialize array to store all dataframes
df_list = []

# for each input file, make dataframe from contents, print number/length, and append to df_list
for i in range(len(sys.argv) - 2):
	temp_df = pd.read_csv(sys.argv[i + 2])
	print('Dataframe %d; Length %d' % (i + 1, len(temp_df)))
	df_list.append(temp_df)

# concatenate dataframes, print, and save
df_compiled = pd.concat(df_list, ignore_index=True)

print(df_compiled)

df_compiled.to_csv(savepath, index=False)

# for compiling tiles
# print('\nTiles before size filtering: %d' % len(df_compiled))
# df_compiled = df_compiled[df_compiled['Tile Sequence'].str.len().ge(40)]
# df_compiled.drop_duplicates('Tile Sequence', keep = 'first', inplace = True)
# print('Tiles after size filtering: %d\n' % len(df_compiled))
