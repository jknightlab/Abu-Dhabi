##################################
# Resubmit failed permutation jobs
##################################

# list successful tasks
find -name 'permresult_gene*' | perl -lape '~s/.*permresult_gene//g;~s/.rds//g' | sort > summary.results.txt
cat summary.results.txt | sort > summary.succeed.txt

# ID failed
echo 21580 | perl -lane '$total=$F[0]-1; for(1..$total){print "$_"}' | sort > summary.all.txt
comm -2 -3 summary.all.txt summary.succeed.txt | sort -k1n > summary.failed.txt
wc -l summary.failed.txt

## If failures are found, resubmit failed jobs
for i in `cat summary.failed.txt`
do
  qsub -t ${i} ../../../Scripts/eQTL_permutation_25PCs_failed.sh
done
