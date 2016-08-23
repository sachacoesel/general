#!/usr/bin/env python

# this script will clean-up the final release of the MMETSP translated transcriptomes
# by pooling similar MMETSP samples at the species level.


# import packages
import os
from Bio import SeqIO


# MMETSP_FASTA_PATH = "/Users/rgroussman/data/MMETSP/CAM_P_0001000.pep.fa" # RDG local
# MMETSP_FASTA_PATH = "/share/data/seq/projects/diatom_est/finalized_data/CAM_P_0001000.pep.fa" # on bloom
# MMETSP_FASTA_PATH = "/Users/rgroussman/data/MMETSP/mini_CAM.fa" # for testing purposes
MMETSP_FASTA_PATH = "/home/rgrous83/sandbox/fnames_test/CAM.test.fasta" # for testing on bloom

# OUT_DIR = "//Users/rgroussman/data/MMETSP/clean_mmetsp/" # local
OUT_DIR = "/share/data/seq/projects/diatom_est/finalized_data/cleaned_mmetsp/" # for bloom


def make_fasta_handle(binomial_fasta):
    """Generates a filename from the binomial name.
    Ex: Genus_species -> Genus_species.fa
    """
    OUT_SUFFIX = ".pep.fa"
    return OUT_DIR + binomial_fasta + OUT_SUFFIX

def open_new_fasta(binomial):
    """Opens up a new output file"""
    fasta_handle = make_fasta_handle(binomial_fasta)
    global out_fasta
    out_fasta = open(fasta_handle,"a")

def close_old_fasta(binomial_fasta):

    """If the given old binomial_fasta exists; close it"""
    fasta_handle = make_fasta_handle(binomial_fasta)

    binomial_fasta_exists = os.path.isfile(fasta_handle) # True or False
    if binomial_fasta_exists == True:
        out_fasta.close()

def write_seq_to_fasta(seq, binomial_fasta):
    """NOTHING YET HAHA"""

    print seq
    print binomial_fasta
    print make_fasta_handle(binomial_fasta)
    print
    SeqIO.write(seq, out_fasta, "fasta")

def parse_mmetsp_header(mmetsp_header):
    """Parses the inconsistent MMETSP headers and attempts to return the
    binomial name found in ORGANISM string. Returns Genus_species string"""

    header_elts = mmetsp_header.split("/")
    for elt in header_elts:
        elt = elt.lstrip("/") # not every field is prefaced by fwd-slash
        if header_elts[2] == 'TAXON_ID=39447 ': # special-case fix
            return "Gymnodinium_catenatum"
        org_index = elt.find("ORGANISM=")
        if org_index != -1: # if the ORGANISM field exists
            new_header = elt[(org_index+10):]
            new_header = new_header.split(",")[0]
            new_header = new_header.replace("\"","")
            new_header = new_header.replace(" ","_")
            new_header = new_header.replace(".","")
            new_header = new_header.strip("_")
            return new_header

# open the MMETSP file
mmetsp_fasta = open(MMETSP_FASTA_PATH, 'r')

for seq in SeqIO.parse(mmetsp_fasta, "fasta"):

    # Get the binomial name from the fasta headers
    binomial = parse_mmetsp_header(seq.description)

    # For the first sequence
    binomial_fasta = False

    if binomial == binomial_fasta:
        write_seq_to_fasta(seq, binomial_fasta)

    elif binomial != binomial_fasta:

        if binomial_fasta != False:
            close_old_fasta(binomial_fasta)
            binomial_fasta = binomial
            open_new_fasta(binomial_fasta)
            write_seq_to_fasta(seq, binomial_fasta)

        elif binomial_fasta == False:
            binomial_fasta = binomial
            open_new_fasta(binomial_fasta)
            write_seq_to_fasta(seq, binomial_fasta)

close_old_fasta(binomial_fasta)
