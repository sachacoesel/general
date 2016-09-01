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
fnames2.py [-tmc] input.fasta

Output:
input.fn.fasta
"""

# load packages
# from os import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input file")
parser.add_argument("-t", "--tax_id", help="Keep TAXON_ID in new headers", action="store_true")
parser.add_argument("-m", "--mmetsp_id", help="Keep MMETSP_ID in new headers", action="store_true")
parser.add_argument("-c", "--campep_id", help="Keep CAMPEP_ID in new headers", action="store_true")
args = parser.parse_args()

global taxon_id
global mmetsp_id
global campep_id

# Make these strings empty if they are not needed
if args.tax_id == False:
    taxon_id = ""
if args.mmetsp_id == False:
    mmetsp_id = ""
if args.campep_id == False:
    campep_id = ""


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

    # This dictionary links JGI genome projects to their NCBI taxonomy ID
    jgi_taxid_dict = {
        "Thaps3":296543,
        "Phatr2":556484,
        "Psemu1":635002,
        "Fracy1":635003,
        "Auran1":44056,
        "Emihu1":2903,
        "Bigna1":227086,
        "Chlre4":3055,
        "OstRCC809_2":385169,
        "Ost9901_3":242159 }

    global taxon_id
    global mmetsp_id
    global campep_id

    header_elts = jgi_header.split("|")
    organism = header_elts[1]
    if args.tax_id == True:
        taxon_id = "_tax" + str(jgi_taxid_dict[organism])
    gene_id = header_elts[2]
    return organism + "_" + str(gene_id) + taxon_id

def parse_mmetsp_header(mmetsp_header):
    """Parses the inconsistent MMETSP headers and attempts to return the ORGANISM string
    along with the unique CAMPEP identifier"""

    global taxon_id
    global mmetsp_id
    global campep_id

    header_elts = mmetsp_header.split("/")

    # Keep the CAMPEP id if requested
    if args.campep_id == True:
        campep_id = "_" + header_elts[0][7:].strip()

    # Look for each desired header element in a non-specific order (they are not precisely ordered in the headers)
    for elt in header_elts:
        elt = elt.lstrip("/") # strip away the fwd-slash because not all fields are prefaced by it

        # Look for the TAXON_ID if asked to
        if args.tax_id == True:
            if elt.find("TAXON_ID=") != -1:
                taxid = elt.split(" ")[0]
                taxon_id = "_tax" + taxid[9:]

        # Look for the MMETSP_ID if asked to
        if args.mmetsp_id == True:
            if elt.find("NCGR_PEP_ID=MMETSP") != -1:
                mmetsp_id = "_" + elt[12:22]

        # special fixes for samples with unusual CAMERA headers
        if header_elts[2] == 'TAXON_ID=39447 ':
            return "Gymnodinium_catenatum_" + campep_id
        elif header_elts[3] == 'TAXON_ID=96639 ':
            return "labyrinthulid_quahog_parasite_QPX" + campep_id

        # always get the binomial name from the ORGANISM field
        org_index = elt.find("ORGANISM=")
        if org_index != -1:
            binomial = elt[(org_index+10):]
            binomial = binomial.split(",")[0]
            binomial = binomial.translate(None, '\".').strip()
            # new_header = new_header.replace("\"","")
            binomial = binomial.replace(" ","_")
            # new_header = new_header.replace(".","")
            # new_header = new_header.strip()
            new_header = binomial + campep_id + mmetsp_id + taxon_id
            return new_header

# load the input fasta and output fasta
input_fasta = open(args.input_fasta, 'r')
outfile_path = build_output_handle(args.input_fasta)
output_fasta = open(outfile_path, 'w')

# for each header, determine the header type, and parse in source-dependent manner

print "Running fnames2.py on", args.input_fasta
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
out_message = "Output: " + outfile_path
if args.campep_id or args.mmetsp_id or args.tax_id == True:
    out_message += " with options: "
    if args.tax_id == True:
        out_message += "--tax_id "
    if args.mmetsp_id == True:
        out_message += "--mmetsp_id "
    if args.campep_id == True:
        out_message += "--campep_id "

print out_message
input_fasta.close()
output_fasta.close()
