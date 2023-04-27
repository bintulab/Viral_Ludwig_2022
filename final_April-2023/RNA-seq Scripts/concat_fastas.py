import sys

fasta1len = 0
fasta2len = 0

with open(sys.argv[3], 'w') as compiled:
	with open(sys.argv[1], 'r') as fasta1:
		for line in fasta1:
			fasta1len += 1
			compiled.write(line)
	with open(sys.argv[2], 'r') as fasta2:
		for line in fasta2:
			fasta2len += 1
			compiled.write(line)

print(fasta1len)
print(fasta2len)
