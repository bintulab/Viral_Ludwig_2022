#Made by ART, 1/10/22
#Use to update genome annotation files for alignment against expressed transgenes (gtf files)

##RUN THIS SCRIPT FROM 'TRANSCRIPTOMES' FOLDER WHERE ALL GTFS ARE STORED
#1st input: fasta file with transgenes
#2nd input: gtf file of entire genome
#3rd input: name of new gtf

import os
import sys
import pandas as pd


namelist = []
seqlist = []
with open(sys.argv[1], 'r') as fasta1:
	for line in fasta1:
		if '>' in line:
			name = line.split('>')[1].strip()
			namelist.append(name)
		else:
			seqlist.append(line.strip())

gtflist = []
for i in range(len(namelist)):
	gtf_cmd = 'echo \'%s\\tunknown\\texon\\t1\\t%s\\t.\\t+\\t.\\tgene_id \"%s\"; transcript_id \"%s\"; gene_biotype \"protein coding\";\'>%s.gtf' %(namelist[i], len(seqlist[i]), namelist[i], namelist[i], namelist[i])
	newgtf = namelist[i]+'.gtf'
	gtflist.append(newgtf)
	print(gtf_cmd)
	os.system(gtf_cmd)

cp_cmd = 'cp %s %s' % (sys.argv[2], sys.argv[3])
print(cp_cmd)
os.system(cp_cmd)

for gtf in gtflist:
	append_cmd = 'cat %s >> %s' % (gtf, sys.argv[3])
	print(append_cmd)
	os.system(append_cmd)

print("New GTF Made")
