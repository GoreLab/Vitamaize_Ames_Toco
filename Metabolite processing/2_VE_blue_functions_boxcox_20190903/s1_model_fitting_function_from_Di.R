# Created by Elodie Gazave; last modified Feb 12, 2016

############################################################
################### aREML model fitting ####################
############################################################

# Input Data structure:
# The order of the columns does not matter as long as the first columns for the variables and the following ones for the traits
# E.g. the 5 first columns of the input file are Genotype, Year, Env, Row, Col
# Then, columns 6 to n are phenotypes (Yield, Height etc). The trait name in the header will be used to create output
# folder, so better not to use long weird names with spaces or special characters.
# Input should be sorted by Env, then Rep, then Row than Col
# That is, the higher level of organization is noted first (within Env we have Reps and cols and rows)
# then the second highest level of organization (within Reps we have cols and rows)
# etc.

library(asreml)


fitmodel <- function(cleanedpheno, curr.trait, random, fixed, nestedterms, PREDICTION, BLUx) {

    ### Here, check size of random and fixed vector
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
                           ", na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    }
    if (length(fixed) > 0) {
        ### Full model with both fixed and random effects -> can be either BLUEs or BLUPs, depending whether Genotype is fixed or random
        fullmodel <- paste( "asreml(fixed = ", curr.trait, " ~ ", paste(fixed,collapse=" + "),
                           ", random = ~", paste(random,collapse=" + "),
                           ", na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    }
    FM.asr <- eval(parse( text = fullmodel))
    FM <- FM.asr$loglik

    # #### Alternative model testing, dropping one term at a time
    # randomInBestModel <- vector()
    # #
    # for (i in c(1:length(random))) {
    #     termToKeep <- seq(1:length(random))
    #     termToKeep <- termToKeep[-i]
    #     tTK <- random[termToKeep]
    #     if (length(fixed) != 0) {
    #         alternatemodel <- paste( "asreml(fixed = ", curr.trait, " ~ ", paste(fixed,collapse=" + "),
    #                        ", random = ~", paste(tTK,collapse=" + "),
    #                        ", na.method.X = 'omit', maxiter=50, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    #     }
    #     if (length(fixed) == 0) {
    #         alternatemodel <- paste( "asreml(fixed = ", curr.trait, " ~ 1 ",
    #                        ", random = ~", paste(tTK,collapse=" + "),
    #                        ", na.method.X = 'omit', maxiter=50, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    #     }
    #     AM.asr <- eval(parse( text = alternatemodel))
    #     AM <- AM.asr$loglik
    #     #
    #     cat('\n\n')
    #     cat('####### Testing significance of terms ',random[i],' for',curr.trait,'#######')
    #     cat('\n\n')
    #     #
    #     effect = 2*(FM-AM)
    #     degreeOfFreedom <- length(random) - length(tTK)  #should be 1
    #     p.effect = pchisq(effect, df=degreeOfFreedom, lower.tail=FALSE)
    #     cat('Effect:', random[i], if(p.effect<0.05){'Significant'}else{'Not Significant'},'\n')
    #     cat('2*(log(full model)-log(simplified model)) =', effect, 'with',degreeOfFreedom,'degree(s) of freedom\n')
    #     cat('p-value = ', p.effect,'\n')
    #     cat('#####################################\n\n')
    #     #
    #     ### Keep the significant terms
    #     if (p.effect < 0.05) {
    #     ### Check for nested effects
    #         if (random[i] %in% nestedterms) {
    #             eachterm <- strsplit(random[i], ':')[[1]]
    #             all_terms <- eachterm[1]
    #             for (j in c(2: (length(eachterm) ) ) ) {
    #                 all_terms <- c(all_terms, paste(tail(all_terms, n=1), eachterm[j], sep=':') )
    #             }
    #             randomInBestModel <- c(randomInBestModel, all_terms)
    #         } else {
    #         randomInBestModel <- c(randomInBestModel, random[i])
    #         }
    #     }
    # }
    # 
    # ### Filter terms that may be in duplicate
    # randomInBestModel <- unique(randomInBestModel)
    # 
    # ### And make sure that all the terms were initially in the model
    # ### Just a sanity check, but in the loop above, we just decompose the interaction and re-create all possible nested terms
    # ### In theory, we want to keep "Year, Year:Loc, Year:Loc:Trt, Year:Loc:Trt:Rep" if Year:Loc:Trt:Rep:Col is significant
    # ### But we want to check that Year:Loc alone (for example) was in the model specified by the user
    # ### Which in this case may not be the case because we are not interested in the variance in Year:Loc, the variance
    # ### in Trt (treatment, wet or dry) being the interesting and most important source of differences
    # randomInBestModel <- randomInBestModel[randomInBestModel %in% random]
    # 
    # ### Check if we kept any random effect
    # if (length(fixed) != 0) {
    #     if (length(randomInBestModel) >= 1) {
    #         ### Best model
    #         bestmodel <- paste( "asreml(fixed = ", curr.trait, " ~", paste(fixed,collapse=" + "),
    #                        ", random = ~", paste(randomInBestModel,collapse=" + "),
    #                        ", na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    #         } else { bestmodel <- paste( "asreml(fixed = ", curr.trait, " ~", paste(fixed,collapse=" + "),
    #                        ", na.method.X = 'omit', maxiter=100, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    #     }
    # }
    # if (length(fixed) == 0) {
    #     if (length(randomInBestModel) >= 1) {
    #         ### Best model
    #         bestmodel <- paste( "asreml(fixed = ", curr.trait, " ~ 1 ",
    #                        ", random = ~", paste(randomInBestModel,collapse=" + "),
    #                        ", na.method.X = 'omit', maxiter=50, data = cleanedpheno, workspace = 1e+09, pworkspace = 1e+09)" )  ### adding "workspace=" by Xiaowei LI
    #         } else {
    #             stop('All terms were dropped. No terms left in the model', call. = F)
    #         }
    # 
    # }
    # 
    # # Output the best model QC and BLUES or BLUPs
    # BM.asr <- eval(parse( text = bestmodel))

    if (length(fixed) != 0) {
        ### Output the fixed effect ANOVA table
        # NB: Using 'numeric' and 'conditional' get us the closest to SAS output
        fixed.effects <- wald.asreml(FM.asr, denDF = 'numeric', ssType = 'conditional', pworkspace=320e6)
        fixed.effect.anova <- data.frame(fixed.effects$Wald)
        write.table(fixed.effect.anova, paste('Fixed_effect_ANOVA_', curr.trait, sep=''), quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
    } else { fixed.effect.anova <- 'There are no fixed effects to be tested in the model' }
    #
    pdf(paste(curr.trait, ".pdf", sep=''))
    plot(FM.asr)
    dev.off()
    #
    predict(FM.asr, classify=PREDICTION, data=cleanedpheno, pworkspace=320e6)-> blux
    #
    pdf(paste('QQplot_', curr.trait, ".pdf", sep=''))
    qqnorm(blux$fitted.values)
    qqline(blux$fitted.values)
    dev.off()
    #
    pdf(paste('Boxplot_', curr.trait, ".pdf", sep=''))
    boxplot(blux$fitted.values)
    dev.off()
    #
    pdf(paste('Histplot_', curr.trait, ".pdf", sep=''))
    hist(blux$fitted.values, freq=FALSE)
    curve(dnorm(x=x, mean=mean(blux$fitted.values, na.rm=TRUE), sd=sd(blux$fitted.values, na.rm=TRUE)),
           from=min(blux$fitted.values, na.rm=TRUE), to=max(blux$fitted.values, na.rm=TRUE),
           add=TRUE, col="blue")
    dev.off()
    #
    # mod.out <- c(Reduce(paste, deparse(BM.asr$fixed.formula)), Reduce(paste, deparse(BM.asr$random.formula)))
    # write.table(paste(mod.out, collapse=' '), paste('Best_model_', curr.trait, sep=''), quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    if (BLUx == 'BLUES') {
        write.table(blux$predictions$pvals, paste('BLUES_', curr.trait, sep=''), quote=FALSE, sep='\t', row.names=FALSE)
        return(list(model = summary(FM.asr), full.model = FM.asr, blues = blux$predictions$pvals, anova = fixed.effect.anova) )
    }
    if (BLUx == 'BLUPS') {
        # Same as BLUEs, except that we need to add the grand mean to the estimates to get the "real" BLUPs
        # Extract the model coefficients
        random.ef <- data.frame(FM.asr$coefficients$random)
        if (nrow(random.ef)>0){
          random.ef[random.ef==0] <- NA
        } # treat missing phenotype as NA, not adding grand mean
        fixed.ef <- FM.asr$coefficients$fixed
        # Get grand mean to add back to BLUPs
        intercept <- fixed.ef[grep('Intercept', names(fixed.ef))]
        BLUPS <- random.ef + intercept
        colnames(BLUPS) <- curr.trait
        write.table(BLUPS, paste('BLUPS_', curr.trait, sep=''), quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
        # We return the list of ALL predictions even for terms removed from the model; these will have NA
        # This is necessary so we can concatenate the results from the different traits, even if the best model is
        # different for different traits.
        # First, get the row.names of all the possible predictions
        all.random.ef <- data.frame(FM.asr$coefficients$random)
        all.predict <- row.names(all.random.ef)
        all.na <- rep(NA, length(all.predict))
        all.na <- as.matrix(all.na, ncol=1)
        row.names(all.na) <- all.predict
        n = 1
        for ( i in row.names(BLUPS) ) {
            if (i %in% row.names(all.na) ) {
                j <- which( row.names(all.na) %in% i )
                all.na[j,] <- BLUPS[n,]
            }
        n = n + 1
        }
        BLUPS <- all.na
        return(list(model = summary(FM.asr), full.model = FM.asr, blups = BLUPS, anova = fixed.effect.anova) )
    }
}
