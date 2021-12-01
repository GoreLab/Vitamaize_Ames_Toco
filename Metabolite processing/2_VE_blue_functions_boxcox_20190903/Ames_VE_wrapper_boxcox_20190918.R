# wrapper_boxcox_2017_02_16   ### firstly, change the title to uppercase format manually and Check the correspondence between the two inputs!!!
####                                                                        ####
# Created by Elodie Gazave; last modified Feb 16, 2017 with James Chamness edits for augmented design
####                                                                        ####

### The script checks data normality, removes outliers, performs model selection and returns BLUEs or BLUPS
### and heritability of all the traits in an input file.

### See README for more details on usage.

### the current version were updated by Xiaowei Li and Di Wu to calculate BLUEs for Ames VE data.

##############################
######### USER INPUT #########
##############################

dat <- read.table('ve_both_years_new_block.txt', header=TRUE, sep="\t", as.is=TRUE)
colnames(dat)
colnames(dat)[1:8] <- c("GENOTYPE","TIER","New_Block_Range","New_Block_Pass","YEAR","CHECK","IS_EXPERIMENTAL","PLATE")
write.table(dat,'input_test_with_checks.txt',sep='\t',row.names = F,quote=F)

# need to check whether the route changed
path.to.res <- '2_ve_blue_functions_boxcox_20190903\\'  
path.to.infile <- '2_ve_blue_functions_boxcox_20190903\\'
infile <- 'input_test_with_checks.txt'
modeltermfile <- 'model_terms_with_checks.txt'

library(asreml)
setwd(path.to.infile)
source('2_ve_blue_functions_boxcox_20190903\\s1_outlier_removal_function.R')
source('2_ve_blue_functions_boxcox_20190903\\s1_model_fitting_function_from_Di.R')
source('2_ve_blue_functions_boxcox_20190903\\s1_heritability_function_modified_xwl.R')

##############################
#### BELOW, DO NOT TOUCH ####
##############################

myvar <- read.table(modeltermfile, sep='$')
fixed = unlist( strsplit( as.character(myvar[1,]), ' ') )
random = unlist( strsplit( as.character(myvar[2,]), ' ') )
nbvariables = unlist( strsplit( as.character(myvar[3,]), ' ') )
nbtraits = unlist( strsplit( as.character(myvar[4,]), ' ') )
PREDICTION = unlist( strsplit( as.character(myvar[5,]), ' ') )
nested = unlist( strsplit( as.character(myvar[6,]), ' ') )
VCterms = unlist( strsplit( as.character(myvar[7,]), ' ') )
BLUx = unlist( strsplit( as.character(myvar[8,]), ' ') )
genoname = unlist( strsplit( as.character(myvar[9,]), ' ') )
baselinemodel = unlist( strsplit( as.character(myvar[10,]), ' ') )
augmented = unlist( strsplit( as.character(myvar[11,]), ' ') )

# extract fixed terms
if (length(fixed) < 3){
  fixed <- vector()
} else { fixed <- fixed[3:length(fixed)] }

# extract random terms
if (length(random) < 3){
  random <- vector()
} else { random <- random[3:length(random)] }

# extract baseline terms
if (length(baselinemodel) < 3){
  baselinemodel <- vector()
} else { baselinemodel <- baselinemodel[3:length(baselinemodel)] }

# extract heritability var
if (length(VCterms) < 3){
  VCterms <- vector()
} else { VCterms <- VCterms[3:length(VCterms)] }

# extract other model terms
PREDICTION = PREDICTION[3]
BLUx = BLUx[3]
genoname = genoname[3]
# format variable names
nbtraits = as.numeric(nbtraits[3])
nbvariables = as.numeric(nbvariables[3])
augmented = as.logical(augmented[3])

setwd(path.to.infile)
pheno <- read.table(infile, header=TRUE, sep="\t", colClasses=c(rep('factor', nbvariables),rep('numeric', nbtraits)), as.is=TRUE)
if (augmented) {
  ### Included for augmented design with checks; the column IS_EXPERIMENTAL should be the only one NOT to be a factor
  pheno$IS_EXPERIMENTAL <- as.numeric(as.character(pheno$IS_EXPERIMENTAL))
}
# transf.dataset <- pheno[,c(1:nbvariables)]
clean.dataset <- pheno[,c(1:nbvariables)]
all.hm <- matrix(data=c('Trait', VCterms), nrow=1, byrow=TRUE)
all.herit <- matrix(data=c('Trait', 'h2plot',  'h2plot.pin', 'seplot.pin'), ncol=4, nrow=1, byrow=TRUE)
# boxcox.transformed <- matrix(data=c('Trait', 'Lambda', 'Transformed'), ncol=3, nrow=1, byrow=TRUE)
alltraitnames <- vector()

for (i in c(1:nbtraits)) {
  
  #### Current trait ####
  curr.trait <- colnames(pheno[ (nbvariables + i) ])
  alltraitnames <- c(alltraitnames, curr.trait)
  shell(paste('mkdir ', path.to.res, curr.trait, sep=''))
  setwd(paste(path.to.res, curr.trait, sep=''))
  cat( paste ('\n\n########################### Examining ', curr.trait, ' ##############################\n\n',sep='') )
  
  ## STEP 1 ## Remove outliers
  initial.outliers(pheno, curr.trait, random, fixed, nbvariables, BLUx) -> out
  as.matrix(out) -> out
  colnames(out) <- curr.trait
  if(i == 1) {
    cleanedpheno <- cbind(clean.dataset, out)
  } else {
    cleanedpheno <- cbind(cleanedpheno, out)
  }
  
  ## STEP 2 ## Find the best model and generate BLUEs or BLUPs
  fitmodel( cleanedpheno,  curr.trait, random, fixed, nested, PREDICTION, BLUx) -> bm
  assign(paste("bm_",curr.trait, sep=''), bm)
  ## For BLUEs, bm is a list with four elements, where:
  ## bm$model is the summary of the model,
  ## bm$best.model is the full list of asreml output objects,
  ## bm$blues are the BLUES
  ## bm$anova are the P-values for fixed effects if there are some
  ## For BLUPs, bm is a list with four elements, where:
  ## bm$model is the summary of the model,
  ## bm$best.model is the full list of asreml output objects,
  ## bm$blups are the BLUPS
  ## bm$anova are the P-values for fixed effects if there are some
  
  ### In addition to writing one BLUXs file in each directory, I also write all BLUXs in one single file for GAPIT or other GWAS
  if (BLUx == 'BLUES') {
    if(i == 1) {
      idx <- which(colnames(bm$blues) == 'predicted.value')
      allbluXs <- bm$blues[,1:idx]
      bluXscolnames <- colnames(bm$blues)[1:idx]
      bluXscolnames[idx] <- curr.trait
    } else {
      allbluXs <- cbind(allbluXs, as.character(bm$blues$predicted.value) )
      bluXscolnames[i+idx-1] <- curr.trait
    }
  }
  if (BLUx == 'BLUPS') {
    if(i == 1) {
      allbluXs <- bm$blups
      bluXscolnames <- c('predict', curr.trait)
    } else {
      allbluXs <- cbind(allbluXs, bm$blups)
      bluXscolnames <- c(bluXscolnames, curr.trait)
    }
  }
  
  ### And all the fixed effects into a single table
  if ( length(bm$anova) > 1 ) {
    Fstat <- bm$anova$F.con[2:length(bm$anova$F.con)]
    Pval <- bm$anova$Pr[2:length(bm$anova$Pr)]
    if (i == 1) {
      fix.header <- vector()
      for (j in c(2: length(rownames(bm$anova) ) ) ) {
        fix.header <- c(fix.header, paste('F.cond', rownames(bm$anova)[j], sep='_'), paste('Pvalue', rownames(bm$anova)[j], sep='_') )
      }
      allanova <- matrix(data=fix.header, ncol=length(fix.header), nrow=1)
    }
    fiXe <- vector()
    for (j in c(1:length(Fstat)) ) {
      fiXe <- c(fiXe, c( Fstat[j], Pval[j]))
    }
    allanova <- rbind(allanova, fiXe)
  }
  
  ## STEP 3 ## Calculate the heritability
  harmonic.mean( cleanedpheno, curr.trait, VCterms, genoname) -> harm.means
  hm <- matrix( data=c(curr.trait, harm.means), nrow=1, ncol=length(harm.means)+1, byrow=TRUE )
  all.hm <- rbind(all.hm, hm)
  fullmodel <- fit.full.model(cleanedpheno, curr.trait, random, fixed)
  ################################################################################
  #### JC: added the 'augmented' parameter
  ################################################################################
  herit.res <- calc.herit(fullmodel, curr.trait, harm.means, augmented)
  ################################################################################
  all.herit <- rbind(all.herit, c(curr.trait, herit.res$h2plot, herit.res$h2plot.pin, herit.res$seplot.pin))
  
}

setwd(path.to.infile)

write.table( cleanedpheno, file=paste('cleaned_', infile, sep=''),  quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
write.table(all.hm, file=paste("HM_", infile, sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
write.table(all.herit, file=paste("Heritability_", infile, sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


if (length(fixed > 0)) {
  alltraitnames <- c ('Trait', alltraitnames)
  allanova <- cbind(alltraitnames, allanova)
  write.table( allanova, file=paste('all_fixed_effects_', infile, sep=''),  quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
}

if (BLUx == 'BLUES') {
  colnames(allbluXs) <- bluXscolnames
  write.table( allbluXs, file=paste('BLUES_all_traits_', infile, sep=''),  quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
}
if (BLUx == 'BLUPS') {
  allbluXs <- cbind(row.names(allbluXs), allbluXs)
  colnames(allbluXs) <- bluXscolnames
  write.table( allbluXs, file=paste('BLUPS_all_traits_', infile, sep=''),  quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
}

#################     g.T not converged using wald.asreml(), so using update.asreml() to fix this bug
FM.asr <- asreml(fixed =  g.T  ~  CHECK  , random = ~ CHECK:YEAR + CHECK:GENOTYPE+ CHECK:YEAR:GENOTYPE:IS_EXPERIMENTAL
                 + YEAR + YEAR:TIER + YEAR:TIER:New_Block_Pass + YEAR:TIER:New_Block_Range + YEAR:PLATE , 
                 na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)
FM.asr2 <- update.asreml(FM.asr, fixed =  g.T  ~  CHECK + CHECK:GENOTYPE , random = ~ CHECK:YEAR + CHECK:YEAR:GENOTYPE:IS_EXPERIMENTAL 
                         + YEAR + YEAR:TIER + YEAR:TIER:New_Block_Pass + YEAR:TIER:New_Block_Range + YEAR:PLATE , 
                         na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)

### done
