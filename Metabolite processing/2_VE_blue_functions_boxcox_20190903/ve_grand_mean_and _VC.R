##############################
######### USER INPUT #########
##############################

dat <- read.table('ve_both_years_new_block.txt', header=TRUE, sep="\t", as.is=TRUE)
colnames(dat)
colnames(dat)[1:8] <- c("GENOTYPE","TIER","New_Block_Range","New_Block_Pass","YEAR","CHECK","IS_EXPERIMENTAL","PLATE")
write.table(dat,'input_test_with_checks.txt',sep='\t',row.names = F,quote=F)

# need to check whether the route changed

path.to.res <-'2_ve_blue_functions_boxcox_20190903\\'  
path.to.infile <- '2_ve_blue_functions_boxcox_20190903\\'
infile <- 'input_test_with_checks.txt'
modeltermfile <- 'model_terms_with_checks.txt'

library(asreml) #[64-bit] need to change the R version to C:\Program Files\R\R-3.3.2, 

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

### get variance components and grand_mean
shell(paste('mkdir ', path.to.res, "variance_and_grand_mean", sep=''))
setwd(paste(path.to.res, "variance_and_grand_mean", sep=''))

grand <- NULL
variance_conponent <- NULL

for (i in c(1:nbtraits)) {
  
  #### Current trait ####
  curr.trait <- colnames(pheno[ (nbvariables + i) ])
  alltraitnames <- c(alltraitnames, curr.trait)
  
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
  
  ###### get the variance component
  
  random_model <- paste( "asreml(fixed = ", curr.trait, "~ 1",
                         ", random = ~", paste(fixed,collapse=" + "), '+', paste(random,collapse=" + "),
                         ", na.method.X = 'omit', maxiter=50, data = cleanedpheno)" )
  RM.asr <- eval(parse( text = random_model))
  var.comps <- summary(RM.asr)$varcomp

  # put variance components into a vector
  VC <- var.comps[,2]
  VC <- matrix(VC, nrow=1)
  colnames(VC) <- row.names(var.comps)
  rownames(VC) <- curr.trait
  variance_conponent <- as.data.frame(rbind(variance_conponent,VC))

  ###### get the grand mean
  fullmodel <- paste( "asreml(fixed = ", curr.trait, " ~ ", paste(fixed,collapse=" + "),
                      ", random = ~", paste(random,collapse=" + "),
                      ", na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" ) 
  FM.asr <- eval(parse( text = fullmodel))
  fixed.ef <- FM.asr$coefficients$fixed
  # Get grand mean to add back to BLUPs
  intercept <- fixed.ef[grep('Intercept', names(fixed.ef))]
  itpt <- as.matrix(intercept)
  rownames(itpt) <- curr.trait
  grand <- rbind(grand,itpt)

}

### write variance components
variance_conponent$trait <- rownames(variance_conponent)
write.csv(variance_conponent[,c(11,1:10)],"ve_variance_conponent.csv",row.names = F)

### write grand mean
colnames(grand) <- "grand_mean"
grand <- as.data.frame(grand)
grand$trait <- row.names(grand)
write.csv(grand[,c(2,1)],"ve_grand_mean.csv",row.names = F)

### finished!
