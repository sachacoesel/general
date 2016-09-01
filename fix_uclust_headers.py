#!/usr/bin/env python

# go from this:
>16362|*|CAMPEP_0203799370 /NCGR_PEP_ID=MMETSP0100_2-20121128|9872_1 /ASSEMBLY_ACC=CAM_ASM_000210 /TAXON_ID=96639 /ORGANISM=" , Strain NY0313808BC1" /LENGTH=30 /DNA_ID= /DNA_START= /DNA_END= /DNA_ORIENTATION=

# to this:
>CAMPEP_0203799370 /NCGR_PEP_ID=MMETSP0100_2-20121128|9872_1 /ASSEMBLY_ACC=CAM_ASM_000210 /TAXON_ID=96639 /ORGANISM=" , Strain NY0313808BC1" /LENGTH=30 /DNA_ID= /DNA_START= /DNA_END= /DNA_ORIENTATION=

# by removing this:
16362|*|

cat MMETSP.nr_clust.pep.fa | sed 's/^>[[:digit:]]\+|\*|/>/g' > MMETSP.nr_clust.pep.fa_temp
rm MMETSP.nr_clust.pep.fa
mv MMETSP.nr_clust.pep.fa_temp MMETSP.nr_clust.pep.fa

cat MMETSP_last_seq.fasta | sed 's/^>[[:digit:]]+|\*|/>/g' > MMETSP_last_seq.fasta.fa_temp
rm MMETSP_last_seq.fasta
mv MMETSP_last_seq.fasta.fa_temp MMETSP_last_seq.fasta

on mac nano:
    find: >\d+\|\*\|
    replace: >
