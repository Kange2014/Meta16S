#!/bin/bash

## reformat outputs of mothur as R script taxonomic_binning.R input

# Remove the columns "label" and "numOtus", and transpose the table so that OTUs are rows and samples are columns.
awk -v FS="\t" -v OFS="\t" '{printf $2;for(i=4;i<=NF;i=i+1) printf "\t"$i; print ""}' all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared > all.opti_mcc.shared
awk -F'\t' '{for(i=1;i<=NF;i=i+1){a[NR,i]=$i}}END{for(j=1;j<=NF;j++){str=a[1,j];for(i=2;i<=NR;i++){str=str "\t" a[i,j]}print str}}' all.opti_mcc.shared > tr.all.opti_mcc.shared

# append the OTU taxonomy to the above table
colnames=`head -1 tr.all.opti_mcc.shared`
colnames=`echo -e ${colnames} "Taxonomy"`
echo ${colnames/Group/#OTUid} > all.OTU.summary.taxonomy
sed -i 's/ /\t/g' all.OTU.summary.taxonomy
awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$0;next}{if($1 in a) print a[$1],$3}' tr.all.opti_mcc.shared all.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy >> all.OTU.summary.taxonomy

# remove marks of confidence "(100)" 
sed 's/[(][0-9]*[)]//g' all.OTU.summary.taxonomy > all.reformat.OTU.summary.taxonomy
