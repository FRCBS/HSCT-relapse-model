# HSCT-relapse-model

This repository contains scripts used for genomic prediction of relapse after HSCT

## code (./src)

`Plink_feature_selection.sh` Shell script for Plink association test leaving each sample out in turn.

`RF_relapse.R` R script for Random Forest fitting and prediction using the same folds as Plink variant selection (above).

