# Script Name: batch_hisat2_alignment.py
# Author: Connor Ludwig, with modifications by Abby Thurm
# Organization: Bintu Lab, Stanford University
# Date: 10/1/2021, most recent update 5/6/2022


# NOTE: need to organize paired read files into a subfolder in a parent 'fastazip' folder

# Pipeline to process and analyze CUT&RUN data, specifically mapping E. coli genomic reads
# Start with fastq files of paired-end reads
# Usage: python {path to cutrun_pipeline.py} {path to directory with zipped fasta files} {path to reference genome directory with root of genome name (e.g. hg19)}
# Example: python win_e/Connor/cutrun_pipeline.py win_e/Connor/TestFolder/fastazip win_e/Connor/hg19_nospacer_bt2/hg19_cit-nospacer-mcherry_at_end

import os
import sys
import pandas as pd
from datetime import datetime

# Assign inputs to variable names
inputDir = sys.argv[1]
inputRoot_refgenome = sys.argv[2]

inputDir_parent = os.path.dirname(os.path.dirname(inputDir))
records = open(os.path.join(inputDir_parent, 'analysis_records.txt'), 'a')

# Initialize
init_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Initializing for hisat2 alignment @ %s ~' % init_time)
print(' ~ Initializing for hisat2 alignment @ %s ~' % init_time, file=records)

# Make directories and get paths; store commands in variables
samDir = 'sam'
samDir_path = os.path.join(inputDir_parent, samDir)
os.mkdir(samDir_path)

# updated 10/3/2021 to output unaligned files in case transgene transcripts are mapped there
# also allows identification of contamination
unalDir = 'unaligned'
unalDir_path = os.path.join(inputDir_parent, unalDir)
os.mkdir(unalDir_path)

bamDir = 'bam'
bamDir_path = os.path.join(inputDir_parent, bamDir)
os.mkdir(bamDir_path)

alnstatsDir = 'alignment_stats'
alnstatsDir_path = os.path.join(inputDir_parent, alnstatsDir)
os.mkdir(alnstatsDir_path)

index_cmd = 'samtools index '

for folder in os.listdir(inputDir):
	# Define temporary path for subdirectory with alignment files
	tempDir_path = os.path.join(inputDir, folder)
	files = []
	for file in os.listdir(tempDir_path):
		files.append(os.path.join(tempDir_path, file))

	root = '_'.join(file.split('_')[:-3])
	sam_name = root + '.sam'
	sam_path = os.path.join(samDir_path, sam_name)
	unal_name = root + '_unaligned.fastq'
	unal_path = os.path.join(unalDir_path, unal_name)

	alnstats = root + '_alnstats.txt'
	alnstats_path = os.path.join(alnstatsDir_path, alnstats)

	print(inputRoot_refgenome, files[0], files[1], sam_path)
	hisat2_cmd = 'hisat2 -x %s -1 %s -2 %s -S %s -p 10 --un-conc %s 2> %s' % (inputRoot_refgenome, files[0], files[1], sam_path, unal_path, alnstats_path)
	hisat2_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Aligning file @ %s ~' % hisat2_time)
	print(' ~ Aligning file @ %s ~' % hisat2_time, file=records)
	print(hisat2_cmd)
	print(hisat2_cmd, file=records)
	os.system(hisat2_cmd)

	# Convert sam to bam file
	bam_name = root + '.sorted.bam'
	bam_path = os.path.join(bamDir_path, bam_name)
	# bam_cmd = 'samtools view -S -b ' + sam_path + ' > ' + bam_path
	bam_cmd = 'samtools sort ' + sam_path + ' -o ' + bam_path # updated 10/3/2021 to reduce later sort steps and extra files

	bam_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Converting SAM file to BAM file @ %s ~' % bam_time)
	print(' ~ Converting SAM file to BAM file @ %s ~' % bam_time, file=records)
	print(bam_cmd)
	print(bam_cmd, file=records)
	os.system(bam_cmd)

# Tell user that processing is complete
finish_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Processing of RNA-seq files completed @ %s ~' % finish_time)
print(' ~ Processing of RNA-seq files completed @ %s ~' % finish_time, file=records)



# Extract read information from the alignment stats files
sampleList = []
unaligned = []
single = []
multi = []

for file in os.listdir(alnstatsDir_path):
	fileID_components = file.split('_')[:-1]
	fileID = '_'.join(fileID_components)
	file_path = os.path.join(alnstatsDir_path, file)
	sampleList.append(fileID)

	with open(file_path, 'r') as inFile:
		for line in inFile:
			if '0 times' in line:
				unaligned.append(int(line.lstrip().split(' (')[0]))
			elif 'exactly 1 time' in line:
				single.append(int(line.lstrip().split(' (')[0]))
			elif '>1 times' in line:
				multi.append(int(line.lstrip().split(' (')[0]))
			else:
				continue

df = pd.DataFrame({'Sample':sampleList1,
				   'Unaligned':unaligned,
				   'Single':single,
				   'Multi':multi})
df['Aligned'] = df['Single'] + df['Multi']

summary_path = os.path.join(inputDir_parent, 'analysis_summary.csv')
df.to_csv(summary_path, index=False)

records.close()

# New line

