#!/bin/bash


# bash script for retrieving list of MMETSP ids from a FASTA file (raw MMETSP headers)
# Will output a file with sequence count for each unique MMETSP sample
# Ryan Groussman, Armbrust Lab. August 2016.

# usage: extract_mmetsp-id_from_fasta.sh input.fasta

SEARCH_STRING="MMETSP[[:digit:]]{4}"
grep ">" $1 | egrep -o $SEARCH_STRING | sort | uniq -c >> "$1.mmetsp-ids"
