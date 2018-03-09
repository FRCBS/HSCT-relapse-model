#!/bin/bash

## =================================================================
## Plink GLM association test for relapse status
## Covariates: donor age, diagnosis, graft
## Leaving each sample out in turn
## =================================================================

module load plink/1.90b4.1

# input path for Plink formatted genotype files
INPUT_PATH='/wrkdir/genotypes'
# input path for Plink formatted phenotype files
PHENO_PATH='/wrkdir/phenotypes'

# output folder
cd /wrkdir/recipient_logistic

# extract sample IDs from Plink ped genotype file
awk '{print $1}' $INPUT_PATH/recipients.ped > sample_ids

# run logistic regression association tests omitting each sample in turn
for SAMPLE in `cat sample_ids`; do
        echo "$SAMPLE $SAMPLE" > tmp_exclude # the left-out sample
        plink --file $INPUT_PATH/recipients --remove tmp_exclude --logistic genotypic --covar $PHENO_PATH/covariates \
	--pheno $PHENO_PATH/relapse.pheno --allow-no-sex -out relapse_sans_$SAMPLE
        rm -f tmp_exclude
done

