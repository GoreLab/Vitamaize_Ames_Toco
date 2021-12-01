#boxcox_outlier_removal_function
# Created by Elodie Gazave; last modified Feb 12, 2016

############################################################
################# aREML outlier removal ####################
############################################################

# Input Data structure:
# The order of the columns does not matter as long as the first columns for the variables and the following ones for the traits
# E.g. the 5 first columns of the input file are Genotype, Year, Env, Row, Col
# Then, columns 6 to n are phenotypes (Yield, Height etc). The trait name in the header will be used to create output
# folder, so better not to use long weird names with spaces or special characters.
# Input should be sorted by Env, then Rep, then Row than Col
# That is, the higher level of nesting is sorted first (within Env we have Reps and cols and rows)
# then the second highest level of nesting (within Reps we have cols and rows)
# etc.

library(asreml)

initial.outliers <- function(pheno, curr.trait, random, fixed, nbvariables, BLUx) {
    originalpheno <- pheno ## keep an intact copy
  # nb.outliers <- 1
  # hold <- vector()
    significance <- vector()
  # iteration <- 0
    ### Open a pdf to show the outlier removal plots
    pdf(paste('outlier_',curr.trait,'_iterations.pdf', sep=''), height=7, width=7) ### length(which(is.na(studentizedRes))) is just to have unique values. Could be 1, 2 3 etc.
    Mycol <- rep('black', dim(pheno)[1] ) # needs to be defined here in case there are no outliers
    #
    ## if (nb.outliers > 0) {  # just to debug

        # Bonferroni corrected cutoff
        nbind <- eval(parse(text = paste('length(pheno[which(!is.na(pheno$',curr.trait,')),]$',curr.trait,')', sep='')))
        threshold = qt(p= 1-0.05/(2*nbind), df=(nbind - (length(fixed) + length(random) + 1) - 1) ) # the +1 is the parameter estimate for the mean (which is implicit in asREML)
        # n - p - 1 where p is the sum of all parameters estimates, mean included; and n-1 is because we remove 1 ind each time to estimate the rest
        # p 373 in the bible
        ### Here, check size of random and fixed vectors
        if (length(random) == 0){
            stop('No random terms in the model. We are not testing fixed effects alone, so there is nothing to test', call. = FALSE)
        }
        ### No fixed effect -> BLUx should be BLUPs
        if (length(fixed) == 0) {
            if ( BLUx != 'BLUPS' ) {
                stop('No fixed terms in the model. You cannot calculate BLUEs', call. = FALSE)
            }
            ### Full model with random effects alone
            fullmodel <- paste( "asreml(fixed = ", curr.trait, " ~ 1 ",
                               ", random = ~", paste(random,collapse=" + "),
                               ", na.method.X = 'include', maxiter=100, data = pheno, aom=TRUE)" )
        }
        if (length(fixed) > 0) {
            ### Full model with both fixed and random effects -> can be either BLUEs or BLUPs, depending whether Genotype is fixed or random
            fullmodel <- paste( "asreml(fixed = ", curr.trait, " ~ ", paste(fixed,collapse=" + "),
                               ", random = ~", paste(random,collapse=" + "),
                               ", na.method.X = 'include', maxiter=100, data = pheno, aom=TRUE, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
        }
        FM.asr <- eval(parse( text = fullmodel))
        dfFM <- FM.asr$nedf

        ### Outlier removal (gives closest results to SAS, except for mega outliers which are a bit more extrem here than in SAS)
        ### It does not matter as ANY method would have picked them up anyway
        ### The crucial part is to have good accuracy and reproducibility for marginal outliers

        stdRes <- resid(FM.asr, type="stdCond") #stdCond is the default if aom=TRUE; values are c("stdCond", "working", "deviance", "pearson", "response")
        studentizedRes <- stdRes / sqrt( (dfFM - stdRes^2) / (dfFM - 1) )

        ### we update nb.outliers
        nb.outliers <- length(which(abs(studentizedRes) > threshold))

        if (nb.outliers > 0) {
          # iteration = iteration + 1
            # Curr_trait columns only
            significance <- studentizedRes # so it will keep NA where there are no estimates
            significance[which(abs(studentizedRes) < threshold)] <- 'pass'
            significance[which(abs(studentizedRes) > threshold)] <- 'outlier_removed'
            ## this is to add those detected in a previous round
            # if (length(hold) > 0){
            #     significance[hold] <- 'outlier_removed'
            # }
            ### Update pheno, putting NA where we detected outliers
            #
            eval(parse(text= paste('pheno[which(abs(studentizedRes) > ', threshold, '),]$', curr.trait, '<- NA', sep='')))
          # hold <- c(hold, which(abs(studentizedRes) > threshold))
            #
            ### Write to pdf
            Mycol[ which(abs(studentizedRes) > threshold) ] <- 'red'
            plot(studentizedRes, col=Mycol)
            abline(h=threshold, col='red')
            abline(h=-threshold, col='red')
        } else {  # once we have removed everything, or if there are no outliers, we stop
            ## Generate a file with residuals
            #
            if (length(significance) == 0){
                significance <- rep('NA', length = length(studentizedRes) ) # we need null vector to write in out file in case we already have no outlier in the first round
            }
        }
            original.trait <- eval(parse(text = paste('originalpheno$', curr.trait, sep='')))
            formated <- cbind(originalpheno[,c(1:nbvariables)], original.trait, studentizedRes, significance)  ### Add the trait
            colnames(formated) <- c(colnames(originalpheno[,c(1:nbvariables)]), curr.trait, 'Studentized_Residuals','Outlier_Test')
            write.table(formated, file=paste('residuals_', curr.trait, sep=''), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
            #
            # Write the residuals without outliers in last page of pdf
        
        dev.off()
    ## ### we also want a final "clean" input file for records
    ## # it's not essential as the filtered object will be used as an input file for model testing
    ## # but we may want to keep record to the "clean", "final" dataset used for analysis
    ## ### The best way to proceed is to keep the trait vector with NA instead of outliers, and at the end of the loop over
    ## ### all traits, concatenate all the traits back together
    ## ### all traits with outliers removed are called "cleaned_<Traitname>" (e.g. cleaned_Height, cleaned_Yield)
    tmp.trait <- eval(parse(text = paste('pheno$', curr.trait, sep='')))
    return(tmp.trait)
}
