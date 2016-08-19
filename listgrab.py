#!/usr/bin/env python

from os import sys
from Bio import SeqIO

def parse_genbank_header(header):
	"""Will parse genbank style header
	e.g; jgi|Thaps3|268147|estExt_thaps1_ua_kg.C_chr_20001
	and return the unique gene id, eg 268147"""

	return header.split("|")[2]

# usage: listgrabpy listfile fastafile
# load master FASTA file, list of seq names
SeqListFile = sys.argv[1]
seqlist = open(SeqListFile, 'r')
FASTAfile = sys.argv[2]
fasta = open(FASTAfile, 'r')

# load outfile
OutFileName = SeqListFile + ".fasta"
OutFile = open(OutFileName,'w')

getlist=[]	# declare list
for line in seqlist:
	line = line.strip()
	getlist.append(line)	# append line to list
print getlist

for record in SeqIO.parse(fasta, 'fasta'):
	clean_id = parse_genbank_header(record.id)
	if clean_id in getlist:
		SeqIO.write(record, OutFile, "fasta")

OutFile.close()
