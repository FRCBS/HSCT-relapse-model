
## =================================================================
## - RandomForest fit and prediction for relapse status using 
##   biallelic SNPs
## - LOOCV folds from feature selection by Plink GLM
## =================================================================


# libraries
library(ranger)
library(data.table)


## -----------------------------------------------------------------
# Relapse prediction model 1
# - Fit RandomForest model using SNPs selected by LOOCV 
#   association test
# - Predict relapse outcome for each sample
## -----------------------------------------------------------------

# read phenotype file
pheno.relapse <- fread('/wrkdir/phenotypes/relapse.pheno', data.table=F)[, 2:3]

# read recipient genotype data
ped.file      <- fread('/wrkdir/genotypes/recipients.raw', data.table=F)
ped.meta      <- ped.file[, 1:6] # metadata
ped.file      <- ped.file[, 7:ncol(ped.file)] # genotypes only

# read plink association test outputs
assoc.file.path <- '/wrkdir/CV_logistic/recipient_logistic'
assoc.files     <- paste(assoc.file.path, list.files(assoc.file.path, pattern='relapse2.*logistic$'), sep='/')
assoc.files.smp <- sapply(assoc.files, function(x) strsplit(strsplit(x, '_sans_')[[1]][2], '.assoc')[[1]][1]) # sample IDs

# ensure sample order matching
pheno.relapse   <- pheno.relapse[match(ped.meta[, 1], pheno.relapse[, 1]), ]
assoc.files.smp <- assoc.files.smp[match(ped.meta[, 1], assoc.files.smp)]
assoc.files     <- assoc.files[match(ped.meta[, 1], assoc.files.smp)]

# remove missing
rem.missing     <- is.na(pheno.relapse[, 2])
pheno.relapse   <- pheno.relapse[!rem.missing, ]
ped.file        <- ped.file[!rem.missing, ]
ped.meta        <- ped.meta[!rem.missing, ]
assoc.files     <- assoc.files[!rem.missing]
assoc.files.smp <- assoc.files.smp[!rem.missing] 

# LOOCV: fit random forest, predict relapse in recipients
relapse.rfmodel <- lapply(1:length(assoc.files), function(i) 
{ # looping through association test result files
  
  print(i) # assoc test file index
  
  # read association test results, select SNPs
  assoc.result      <- fread(assoc.files[i], data.table=F) 
  assoc.result      <- assoc.result[!is.na(assoc.result[, 'P']), ]
  assoc.result      <- assoc.result[which(assoc.result[, 'TEST']=='ADD' | assoc.result[, 'TEST']=='DOMDEV'), ]
  p.ind             <- which(assoc.result[, 'P'] < 1e-3) # indices of SNPs under p-value threshold
  assoc.result.sel  <- assoc.result[p.ind, ] # select SNPs under the threshold
  
  # pick selected SNPs from genotype ped file
  tmp.cln <- gsub('_1', '_', colnames(ped.file))
  tmp.cln <- gsub('_HET', '_', tmp.cln)
  ped.snp.ind <- tmp.cln %in% paste(assoc.result.sel[, 'SNP'], '_', sep='')
  ped.subset  <- ped.file[, ped.snp.ind] # select variants
  
  # genotype data for training and predicting
  sans.sample <- unname(assoc.files.smp[i]) # left-out sample ID
  ped.train   <- ped.subset[-match(sans.sample, ped.meta[, 1]), ] # remove left-out sample from genotype file for training
  ped.test    <- ped.subset[ rep(match(sans.sample, ped.meta[, 1]), 2), ] # include left-out sample for prediction
  
  # relapse phenotype data
  pheno.relapse.train <- factor(pheno.relapse[-match(sans.sample, pheno.relapse[, 1]), 2])
  pheno.relapse.test  <- pheno.relapse[ match(sans.sample, pheno.relapse[, 1]), 2]
  
  # fit the random forest model
  c1w       <- 1-(sum(pheno.relapse.train==1) / length(pheno.relapse.train)) # set case weights
  c2w       <- 1-(sum(pheno.relapse.train==2) / length(pheno.relapse.train))
  cweights  <- sapply(pheno.relapse.train, function(x) ifelse(x==1, c1w, c2w))
  rf.model  <- ranger(PHE ~., data=data.frame(ped.train, PHE=pheno.relapse.train), 
                      case.weights=cweights, probability=T, num.trees=2500, importance='permutation')
  
  # write variable importances
  write.table(rf.model$variable.importance, paste('sans', sans.sample, '.RF-importances2', sep=''), quote=F, sep='\t')
  
  # predict the left-out sample using the model; selects prediction prob. for relapse
  out <- matrix(c(as.vector(predict(rf.model, data.frame(ped.test)[1, ])$predictions[2]), pheno.relapse.test), ncol=2) 
  rownames(out) <- sans.sample
  colnames(out) <- c('predicted', 'actual')
  
  out # return prediction and actual outcome for the left-out sample
})
relapse.rfmodel <- do.call(rbind, relapse.rfmodel) 

# write output
write.table(relapse.rfmodel, 'relapse2_rf-model_cv-predictions', sep='\t', quote=F, col.names=T, row.names=F)



## -----------------------------------------------------------------
# Relapse prediction model 2
# - Include all covariates to prediction
# - Fit RandomForest model using SNPs selected by association test
# - Predict relapse outcome for each sample
## -----------------------------------------------------------------

# read phenotype file
pheno.relapse <- fread('/wrkdir/phenotypes/relapse.pheno', data.table=F)[, 2:3]

# read covariate data
# graft NAs replaced with 0.5
covar.relapse <- read.delim('/wrkdir/phenotypes/all.covariates', stringsAsFactors=F)

# read recipient genotype data
ped.file      <- fread('/wrkdir/genotypes/recipients.raw', data.table=F)
ped.meta      <- ped.file[, 1:6] # metadata
ped.file      <- ped.file[, 7:ncol(ped.file)] # genotypes only

# read plink association test outputs
assoc.file.path <- '/wrkdir/CV_logistic/recipient_logistic'
assoc.files     <- paste(assoc.file.path, list.files(assoc.file.path, pattern='relapse2.*logistic$'), sep='/')
assoc.files.smp <- sapply(assoc.files, function(x) strsplit(strsplit(x, '_sans_')[[1]][2], '.assoc')[[1]][1]) # sample IDs

# ensure sample order matching
pheno.relapse   <- pheno.relapse[match(ped.meta[, 1], pheno.relapse[, 1]), ]
assoc.files.smp <- assoc.files.smp[match(ped.meta[, 1], assoc.files.smp)]
assoc.files     <- assoc.files[match(ped.meta[, 1], assoc.files.smp)]
covar.relapse   <- covar.relapse[match(ped.meta[, 1], covar.relapse[, 1]), ]

# remove missing
rem.missing     <- is.na(pheno.relapse[, 2])
pheno.relapse   <- pheno.relapse[!rem.missing, ]
ped.file        <- ped.file[!rem.missing, ]
ped.meta        <- ped.meta[!rem.missing, ]
assoc.files     <- assoc.files[!rem.missing]
assoc.files.smp <- assoc.files.smp[!rem.missing] 
covar.relapse   <- covar.relapse[!rem.missing, -1]
  
# LOOCV: fit random forest, predict relapse in recipients
relapse.rfmodel <- lapply(1:length(assoc.files), function(i) 
{ # looping through association test result files
  
  print(i) # assoc test file index
  
  # read association test results, select SNPs
  assoc.result      <- fread(assoc.files[i], data.table=F) 
  assoc.result      <- assoc.result[!is.na(assoc.result[, 'P']), ]
  assoc.result      <- assoc.result[which(assoc.result[, 'TEST']=='ADD' | assoc.result[, 'TEST']=='DOMDEV'), ]
  p.ind             <- which(assoc.result[, 'P'] < 1e-3) # indices of SNPs under p-value threshold
  assoc.result.sel  <- assoc.result[p.ind, ] # select SNPs under the threshold
  
  # pick selected SNPs from genotype ped file
  tmp.cln <- gsub('_1', '_', colnames(ped.file))
  tmp.cln <- gsub('_HET', '_', tmp.cln)
  ped.snp.ind <- tmp.cln %in% paste(assoc.result.sel[, 'SNP'], '_', sep='')
  ped.subset  <- ped.file[, ped.snp.ind] # select variants
  
  # join covariates to genotype data
  ped.subset  <- cbind(ped.subset, covar.relapse) 
  
  # genotype data for training and predicting
  sans.sample <- unname(assoc.files.smp[i]) # left-out sample ID
  ped.train   <- ped.subset[-match(sans.sample, ped.meta[, 1]), ] # remove left-out sample from genotype file for training
  ped.test    <- ped.subset[ rep(match(sans.sample, ped.meta[, 1]), 2), ] # include left-out sample for prediction
  
  # relapse phenotype data
  pheno.relapse.train <- factor(pheno.relapse[-match(sans.sample, pheno.relapse[, 1]), 2])
  pheno.relapse.test  <- pheno.relapse[ match(sans.sample, pheno.relapse[, 1]), 2]
  
  # fit the random forest model
  c1w       <- 1-(sum(pheno.relapse.train==1) / length(pheno.relapse.train)) # set case weights
  c2w       <- 1-(sum(pheno.relapse.train==2) / length(pheno.relapse.train))
  cweights  <- sapply(pheno.relapse.train, function(x) ifelse(x==1, c1w, c2w))
  rf.model  <- ranger(PHE ~., data=data.frame(ped.train, PHE=pheno.relapse.train), 
                      case.weights=cweights, probability=T, num.trees=2500, importance='permutation')
  
  # write variable importances
  write.table(rf.model$variable.importance, paste('sans', sans.sample, '.RF-importances2cov', sep=''), quote=F, sep='\t')
  
  # predict the left-out sample using the model; selects prediction prob. for relapse
  out <- matrix(c(as.vector(predict(rf.model, data.frame(ped.test)[1, ])$predictions[2]), pheno.relapse.test), ncol=2) 
  rownames(out) <- sans.sample
  colnames(out) <- c('predicted', 'actual')
  
  out # return prediction and actual outcome for the left-out sample
})
relapse.rfmodel <- do.call(rbind, relapse.rfmodel) 

# write output
write.table(relapse.rfmodel, 'relapse2cov_rf-model_cv-predictions', sep='\t', quote=F, col.names=T, row.names=F)



## -----------------------------------------------------------------
# Relapse prediction model 3 
# - Include only covariates to prediction
# - Fit RandomForest model
# - Predict relapse outcome for each sample
## -----------------------------------------------------------------

# read phenotype file
pheno.relapse <- fread('/wrkdir/phenotypes/relapse.pheno', data.table=F)[, 2:3]

# read covariate data
# graft NAs replaced with 0.5
covar.relapse <- read.delim('/wrkdir/phenotypes/2abc.all.covariates', stringsAsFactors=F)

# ensure sample order matching
covar.relapse   <- covar.relapse[match(pheno.relapse[, 1], covar.relapse[, 1]), ]

# remove missing
rem.missing     <- is.na(pheno.relapse[, 2])
pheno.relapse   <- pheno.relapse[!rem.missing, ]
covar.relapse   <- covar.relapse[!rem.missing, ]

# LOOCV: fit random forest, predict relapse in recipients
relapse.rfmodel <- lapply(1:nrow(pheno.relapse), function(i) 
{ # looping through association test result files
  
  print(i) # assoc test file index
  
  # genotype data for training and predicting
  sans.sample <- unname(pheno.relapse[i, 1]) # left-out sample ID
  
  # relapse phenotype data
  pheno.relapse.train <- factor(pheno.relapse[-match(sans.sample, pheno.relapse[, 1]), 2])
  pheno.relapse.test  <- pheno.relapse[ match(sans.sample, pheno.relapse[, 1]), 2]
  
  covar.relapse.train <- covar.relapse[-match(sans.sample, covar.relapse[, 1]), ]
  covar.relapse.test  <- covar.relapse[ match(sans.sample, covar.relapse[, 1]), ]
  covar.relapse.train <- covar.relapse.train[, -1]
  covar.relapse.test  <- covar.relapse.test[, -1]
  
  # fit the random forest model
  c1w       <- 1-(sum(pheno.relapse.train==1) / length(pheno.relapse.train)) # set case weights
  c2w       <- 1-(sum(pheno.relapse.train==2) / length(pheno.relapse.train))
  cweights  <- sapply(pheno.relapse.train, function(x) ifelse(x==1, c1w, c2w))
  rf.model  <- ranger(PHE ~., data=data.frame(covar.relapse.train, PHE=pheno.relapse.train), 
                      case.weights=cweights, probability=T, num.trees=2500, importance='permutation')
  
  # write variable importances
  write.table(rf.model$variable.importance, paste('sans', sans.sample, '.RF-importancesCovOnly', sep=''), quote=F, sep='\t')
  
  # predict the left-out sample using the model; selects prediction prob. for relapse
  out <- matrix(c(as.vector(predict(rf.model, data.frame(covar.relapse.test)[1, ])$predictions[2]), pheno.relapse.test), ncol=2) 
  rownames(out) <- sans.sample
  colnames(out) <- c('predicted', 'actual')
  
  out # return prediction and actual outcome for the sample
})
relapse.rfmodel <- do.call(rbind, relapse.rfmodel) 

# write output
write.table(relapse.rfmodel, 'relapseCovOnly_rf-model_cv-predictions', sep='\t', quote=F, col.names=T, row.names=F)



