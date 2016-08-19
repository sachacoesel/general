#!/usr/bin/env python

# link_CAMPEP_to_MMETSP_id.py

# This script will process the final MMETSP peptide assemblies,
# linking each CAMERA accession number to one MMETSP id

from Bio import SeqIO


"""
Example header:
>CAMPEP_0114621474 /NCGR_PEP_ID=MMETSP0168-20121206|9247_1 /TAXON_ID=95228 ORGANISM="Vannella sp., Strain DIVA3 517/6/12" /NCGR_SAMPLE_ID=MMETSP0168 /ASSEMBLY_ACC=CAM_ASM_000044 /LENGTH=342 /DNA_ID=CAMNT_0001832673 /DNA_START=127 /DNA_END=1149 /DNA_ORIENTATION=-
We're only interested in the CAMPEP# and the NCGR_SAMPLE_ID... maybe let's get TAXON_ID and organism string too
NOTE that not all MMETSP headers have the NCGR_SAMPLE_ID, so we're going to go with NCGR_PEP_ID instead :(
There are also irregularities in the fasta headers ABOUNDING... this makes it difficult! For example,
the Vannella 'ORGANISM' field is not prefixed by a '/'.
"""

MMETSP_pep_path = "/Users/rgroussman/data/MMETSP/CAM_P_0001000.pep.fa"
# MMETSP_pep_path = "/Users/rgroussman/data/MMETSP/mini_CAM.fa" # for testing...

output_path = "/Users/rgroussman/data/MMETSP/CAMPEP_to_MMETSP.csv"
output_handle = open(output_path, "w")

def get_ncgr_elts(ncgr_elts):
    """Example:
    'NCGR_PEP_ID=MMETSP0168-20121206|4096_1 '
    """
    ncgr_elts = ncgr_elts.split("|")
    mmetsp_id = ncgr_elts[0].split("-")
    mmetsp_id = mmetsp_id[0][12:]
    pep_id = ncgr_elts[1].strip()

    return mmetsp_id, pep_id

for seq_record in SeqIO.parse(MMETSP_pep_path, "fasta"):
    header_elts = seq_record.description.split("/")
    campep = header_elts[0].strip()
    ncgr_elts = header_elts[1]
    mmetsp_id, pep_id = get_ncgr_elts(ncgr_elts)
    taxon_id = header_elts[2][9:]
    organism = header_elts[3][10:]

    outline = campep + "," + mmetsp_id + "," + pep_id + "\n"
    output_handle.write(outline)
