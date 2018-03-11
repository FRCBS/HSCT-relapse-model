# HSCT-relapse-model

This repository contains the scripts for feature selection and RF model fitting for prediction of relapse after HSCT as used in the manuscript:

J Ritari, K Hyvärinen, S Koskela, M Itälä-Remes, R Niittyvuopio, A Nihtinen, U salmenniemi, M Putkonen, L Volin, T Kwan, T Pastinen, and J Partanen. (2018). Genomic prediction of relapse in recipients of allogeneic haematopoietic stem cell transplantation. _Submitted_.


## code (./src)

`Plink_feature_selection.sh` Shell script for Plink association test leaving each sample out in turn.

`RF_relapse.R` R script for Random Forest fitting and prediction using the same folds as Plink variant selection (above).

