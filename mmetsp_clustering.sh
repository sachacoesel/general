
# pull out the MMETSP id from a file
SEARCH_STRING="NCGR_PEP_ID=MMETSP[[:digit:]]{4}"
for fasta in $(ls -a *.pep.fa); do
  grep ">" $fasta | egrep -o $SEARCH_STRING | sed 's/NCGR_PEP_ID=//g' | sort | uniq -c >> "$fasta.mmetsp-ids"
done

# write out the sequence count for special cases
for special in $(less special_cases.txt); do
  echo $special
  cat $special.mmetsp-ids
done


# write out the tax_ids from a fasta file
SEARCH_TAXON="TAXON_ID=[[:digit:]]*"
for fasta in $(ls -a *.pep.fa); do
  grep ">" $fasta | egrep -o $SEARCH_TAXON | sort | uniq -c >> "$fasta.taxid"
done


# write out the taxid for special cases
SEARCH_TAXON="TAXON_ID=[[:digit:]]*"
for special in $(less special_cases.txt); do
  echo $special >> special_cases.tax-ids
  grep ">" $special | egrep -o $SEARCH_TAXON | sort | uniq -c >> special_cases.tax-ids
done
