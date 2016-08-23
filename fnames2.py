#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab, University of Washington
# July 2016

"""
fnames2.py : Friendly Names 2.0

Purpose: will parse the headers of FASTA-formatted files to make
human readable, e.g: >Genus_species_ID12345
Currently supports MMETSP and JGI headers. Unlike previous fnames.py, this script
will attempt to auto-detect header style.

Usage:
fnames2.py coolgenes.fasta

Output:
coolgenes.fn.fasta
"""

# load packages
from os import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def usage():
    print "usage: fnames2.py input.fasta"
    print "- Currently supports FASTA files with JGI or MMETSP formatted headers"
    print "- Output is input.fn.fasta"
    sys.exit()

def build_output_handle(infile_path):
    handle_elts = infile_path.split(".")
    handle_elts.insert(-1,"fn")
    outfile_path = ".".join(handle_elts)
    return outfile_path

def determine_type(seq_header):
    """Determines from a fasta header whether it is from the JGI or MMETSP collections.
    JGI example: jgi|Thaps3|23115|estExt_fgenesh1_pg.C_chr_60346
    MMETSP example: CAMPEP_0172341316
    Returns 'jgi' or 'mmetsp'"""

    if seq_header.startswith("jgi"):
        return "jgi"
    elif seq_header.startswith("CAMPEP"):
        return "mmetsp"

def parse_jgi_header(jgi_header):
    """Parses the JGI headers, returning the organism code (ex, Emihu1) and gene/protein ID"""

    header_elts = jgi_header.split("|")
    organism = header_elts[1]
    gene_id = header_elts[2]
    return organism + "_" + str(gene_id)

def parse_mmetsp_header(mmetsp_header):
    """Parses the inconsistent MMETSP headers and attempts to return the ORGANISM string
    along with the unique CAMPEP identifier"""

    header_elts = mmetsp_header.split("/")
    for elt in header_elts:
        elt = elt.lstrip("/")
        if header_elts[2] == 'TAXON_ID=39447 ':
            return "Gymnodinium_catenatum_" + header_elts[0][7:]
        org_index = elt.find("ORGANISM=")
        if org_index != -1:
            new_header = elt[(org_index+10):]
            new_header = new_header.split(",")[0]
            new_header = new_header.replace("\"","")
            new_header = new_header.replace(" ","_")
            new_header = new_header.replace(".","")
            new_header = new_header.strip("_")
            new_header = new_header + "_" + header_elts[0][7:]
            return new_header

if len(sys.argv) != 2:
    usage()
infile_path = sys.argv[1]

# load the input fasta and output fasta
input_fasta = open(infile_path, 'r')
outfile_path = build_output_handle(infile_path)
output_fasta = open(outfile_path, 'w')

# for each header, determine the header type, and parse in source-dependent manner

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    db_type = determine_type(seq_record.id)
    # print seq_record.description
    if db_type == "jgi":
        new_header = parse_jgi_header(seq_record.description)
    elif db_type == "mmetsp":
        new_header = parse_mmetsp_header(seq_record.description)
    # write out seq_record and sequence
    out_record = SeqRecord(seq_record.seq, id= new_header, description="")
    SeqIO.write(out_record, output_fasta, "fasta")

input_fasta.close()
output_fasta.close()
