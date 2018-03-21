# HSCT-relapse-model

This repository contains the scripts for feature selection and RF model fitting for prediction of relapse after HSCT as used in the manuscript:

Ritari J, Hyvärinen K, Koskela S, Itälä-Remes M, Niittyvuopio R, Nihtinen A, Salmenniemi U, Putkonen M, Volin L, Kwan T,  Pastinen T, and Partanen J (2018). Genomic prediction of relapse in recipients of allogeneic haematopoietic stem cell transplantation. _Submitted_.


## code (./src)

`Plink_feature_selection.sh` Shell script for Plink association test leaving each sample out in turn.

`RF_relapse.R` R script for Random Forest fitting and prediction using the same folds as in Plink variant selection (above).

