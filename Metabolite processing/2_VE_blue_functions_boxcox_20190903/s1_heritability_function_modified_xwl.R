#boxcox_heritability_function
# Created by Elodie Gazave; last modified June 1, 2016

# Modified by James Chamness; last modification Oct 31, 2016
# Modification to allow use of augmented design; see in-code documentation
# and an addendum to the original readme file

############################################################
################# Heritability and s.e. ####################
############################################################

## This script is largely based on Duke Pauli's script
## Takes cleaned data as an input, as well as the Var Comp terms to extract


library(asreml)
library(reshape2)
library(reshape)
library(lattice)

fit.full.model <- function(cleanedpheno, curr.trait, random, fixed) {
  fullmodel <- paste( "asreml(fixed = ", curr.trait, "~ 1",
                      ", random = ~", paste(fixed,collapse=" + "), '+', paste(random,collapse=" + "),
                      ", na.method.X = 'omit', maxiter=50, data = cleanedpheno)" )
  FM.asr <- eval(parse( text = fullmodel))
  return(FM.asr)
}

#################################################################################
# calculate harmonic mean
harm.mean <- function(x){
  harmean <- 1/mean(1/x)
  return(harmean)
}

### Loop for curr.trait
harmonic.mean <- function(cleanedpheno, curr.trait, VCterms, genoname) {
  # To count the number of non-missing data per trait, we use the table() function
  # but to have count without missing, we temporarily remove the missing
  # so table() will output the actual valid data
  # We calculate the harmonic mean number of genotypes evaluated for different
  # variance component.
  # For the residual variance, we want the harmonic mean of all the genocode evaluated
  # across all the years, fields, conditions etc
  # For GXE terms, we want the harmonic mean number of environment per genocodes
  # These levels correspond to the way we want to partition the genetic variance component
  
  all.harm.mean <- vector()
  all.harm.mean.names <- vector()
  for (vct in VCterms) {  ## VCterms is a vector with (Genocode, Env)
    vctcol = unlist(strsplit(vct, ':'))
    if (length(vctcol) == 1 ) {
      #################################################################
      ### For ALL
      #################################################################
      if (!genoname %in% vctcol) {
        stop('One of the variance component terms you defined has no genotype. This wrapper only calculate heritability, not variance component analysis', call. = FALSE)
      }
      # Calculate harmonic mean for GENOCODE (the average number of times each genotype was evaluated in total)
      curr.trait.tmp <- cleanedpheno[ c(vctcol, curr.trait)]
      if ( length( which( is.na( curr.trait.tmp[curr.trait]) ) ) > 0 ) {
        curr.trait.tmp <- curr.trait.tmp [-which (is.na (curr.trait.tmp[curr.trait]) ), ]
      }
      count.geno.all <- table(curr.trait.tmp[,1])
      #If some genotypes have counts of '0' as it will return an error in 1/x in harmonic mean function, so we need to remove them
      if ( length (which(count.geno.all == 0 ) ) ) {
        count.geno.all <- count.geno.all[-which(count.geno.all == 0)]
      }
      #Take the harmonic mean for "plots"
      harm.mean.line <- harm.mean(count.geno.all)
      #assign(vct, harm.mean.line)
      ## Put the harmonic mean for plots and ENV onto one vector, with the first column being the trait.
      ## Not really necessary, just as a control to see the values
      all.harm.mean <- c(all.harm.mean, harm.mean.line)
      ## Here the harm mean of all lines (all genotypes) is called Residuals because we use this value to divide the Var Comp of Residuals
      all.harm.mean.names <- c(all.harm.mean.names, 'Residuals')
      
    } else if ( length(vctcol) == 2 ){
      ########################################################################
      ### For TWO-WAY
      ########################################################################
      # Calculate harmonic mean for GENOCODE in ENV (the average number of environments each genotype was evaluated in)
      if (!genoname %in% vctcol) {
        stop('One of the variance component terms you defined has no genotype. This wrapper only calculate heritability, not variance component analysis', call. = FALSE)
      }
      curr.trait.tmp <- cleanedpheno[ c(vctcol, curr.trait)]
      # Combine GENOCODE and the genetic variance partitioning variable (e.g. here, ENV)
      # to form unique entries for each value of trait
      # We don't want the number of unique genocode or the number of unique env
      # we want the number of unique genocode in env
      # e.g.   Genocode       20133 20141 20142
      #        USDA-SBN-001     1     2     2
      #        USDA-SBN-002     1     2     1
      #        USDA-SBN-003     1     2     2
      #        USDA-SBN-004     1     2     1
      melted.dat <- melt(curr.trait.tmp, id=c(vctcol[1], vctcol[2]), measure=curr.trait)
      # Remove NA because length(c(1,2,NA)) returns 3, so if an environment has no data, it will still counts it
      melted.dat <- melted.dat[which(!is.na(melted.dat$value)),]
      ### To cast, we want to have the genotype variable before the ' ~ ' because we always want to count the environment variables for each genotype
      ### Genotype is genoname
      tmpvctcol <- vctcol[-which(vctcol == genoname)]
      casted.dat <- cast(melted.dat, eval(parse(text = paste(genoname, '~', tmpvctcol) )), fun.aggregate =  length )
      # replace zeros with NAs. It's just a trick because it's easier to distinguish NAs than 0's from other numbers and therefore count them
      casted.dat[casted.dat == 0] <- NA
      # new vector with counts of years lines were observed
      count.not.NA <- apply (casted.dat, 1, function(x) length(which( !is.na( x[ c(1:dim(casted.dat)[2]) ] ) )) )
      # get harmonic mean for ENV
      harm.mean.GbyE <- harm.mean(count.not.NA) #aka in Duke's code treatment.harm.mean
      tmpname <- unlist(strsplit(vct, ':')) ### index in vct by length, and not "vct"
      tmpname <- paste(tmpname[1], tmpname[2], sep='_')
      # assign(tmpname, harm.mean.GbyE)
      all.harm.mean <- c(all.harm.mean, harm.mean.GbyE)
      all.harm.mean.names <- c(all.harm.mean.names, gsub("[:]","_", vct))
      
    } else if ( length(vctcol) == 3 ){
      ########################################################################
      ### For THREE-WAY
      ########################################################################
      if (!genoname %in% vctcol) {
        stop('One of the variance component terms you defined has no genotype. This wrapper only calculate heritability, not variance component analysis', call. = FALSE)
      }
      # Calculate harmonic mean for GENOCODE in ENV (the average number of environments each genotype was evaluated in)
      curr.trait.tmp <- cleanedpheno[ c(vctcol, curr.trait)]
      # Combine GENOCODE and the genetic variance partitioning variable (e.g. here, ENV)
      # to form unique entries for each value of trait
      # We don't want the number of unique genocode or the number of unique env
      # we want the number of unique genocode in env
      # e.g.   Genocode       20133 20141 20142
      #        USDA-SBN-001     1     2     2
      #        USDA-SBN-002     1     2     1
      #        USDA-SBN-003     1     2     2
      #        USDA-SBN-004     1     2     1
      melted.dat <- melt(curr.trait.tmp, id=c(vctcol[1], vctcol[2], vctcol[3]), measure=curr.trait)  ### Now do for 2 interactions, but can generalize for 3
      # Remove NA because length(c(1,2,NA)) returns 3, so if an environment has no data, it will still counts it
      melted.dat <- melted.dat[which(!is.na(melted.dat$value)),]
      ### To cast, we want to have the genotype variable before the ' ~ ' because we always want to count the environment variables for each genotype
      ### Genotype is genoname
      tmpvctcol <- vctcol[-which(vctcol == genoname)]
      casted.dat <- cast ( melted.dat, eval( parse(text = paste(genoname, '~', tmpvctcol[1], '+', tmpvctcol[2] ))), fun.aggregate =  length )
      # here is where I need to work to see how we would cast for 3 way
      # replace zeros with NAs. It's just a trick because it's easier to distinguish NAs than 0's from other numbers and therefore count them
      casted.dat[casted.dat == 0] <- NA
      # new vector with counts of years lines were observed
      count.not.NA <- apply (casted.dat, 1, function(x) length(which( !is.na( x[ c(1:dim(casted.dat)[2]) ] ) )) )
      # get harmonic mean for ENV
      harm.mean.GbyEbyE <- harm.mean(count.not.NA) #aka in Duke's code treatment.harm.mean
      tmpname <- unlist(strsplit(vct, ':')) ### index in vct by length, and not "vct"
      tmpname <- paste(tmpname[1], tmpname[2], tmpname[3], sep='_')
      # assign(tmpname, harm.mean.GbyEbyE)
      all.harm.mean <- c(all.harm.mean, harm.mean.GbyEbyE)
      all.harm.mean.names <- c(all.harm.mean.names, gsub("[:]","_", vct))
    } else if ( length(vctcol) == 4 ){  ### adding "length(vctcol) == 4" by Xiaowei Li
      ########################################################################
      ### For FOUR-WAY
      ########################################################################
      if (!genoname %in% vctcol) {
        stop('One of the variance component terms you defined has no genotype. This wrapper only calculate heritability, not variance component analysis', call. = FALSE)
      }
      # Calculate harmonic mean for GENOCODE in ENV (the average number of environments each genotype was evaluated in)
      curr.trait.tmp <- cleanedpheno[ c(vctcol, curr.trait)]
      # Combine GENOCODE and the genetic variance partitioning variable (e.g. here, ENV)
      # to form unique entries for each value of trait
      # We don't want the number of unique genocode or the number of unique env
      # we want the number of unique genocode in env
      # e.g.   Genocode       20133 20141 20142
      #        USDA-SBN-001     1     2     2
      #        USDA-SBN-002     1     2     1
      #        USDA-SBN-003     1     2     2
      #        USDA-SBN-004     1     2     1
      melted.dat <- melt(curr.trait.tmp, id=c(vctcol[1], vctcol[2], vctcol[3],vctcol[4]), measure=curr.trait)  ### Now do for 2 interactions, but can generalize for 3
      # Remove NA because length(c(1,2,NA)) returns 3, so if an environment has no data, it will still counts it
      melted.dat <- melted.dat[which(!is.na(melted.dat$value)),]
      ### To cast, we want to have the genotype variable before the ' ~ ' because we always want to count the environment variables for each genotype
      ### Genotype is genoname
      tmpvctcol <- vctcol[-which(vctcol == genoname)]
      casted.dat <- cast ( melted.dat, eval( parse(text = paste(genoname, '~', tmpvctcol[1], '+', tmpvctcol[2], '+', tmpvctcol[3] ))), fun.aggregate =  length )
      # here is where I need to work to see how we would cast for 3 way
      # replace zeros with NAs. It's just a trick because it's easier to distinguish NAs than 0's from other numbers and therefore count them
      casted.dat[casted.dat == 0] <- NA
      # new vector with counts of years lines were observed
      count.not.NA <- apply (casted.dat, 1, function(x) length(which( !is.na( x[ c(1:dim(casted.dat)[2]) ] ) )) )
      # get harmonic mean for ENV
      harm.mean.GbyEbyE <- harm.mean(count.not.NA) #aka in Duke's code treatment.harm.mean
      tmpname <- unlist(strsplit(vct, ':')) ### index in vct by length, and not "vct"
      tmpname <- paste(tmpname[1], tmpname[2], tmpname[3],tmpname[4], sep='_')
      # assign(tmpname, harm.mean.GbyEbyE)
      all.harm.mean <- c(all.harm.mean, harm.mean.GbyEbyE)
      all.harm.mean.names <- c(all.harm.mean.names, gsub("[:]","_", vct))
    }
    matrix(all.harm.mean, nrow=1) -> all.harm.mean
    colnames(all.harm.mean) <- all.harm.mean.names
  }
  return( all.harm.mean )
}

################################################################################
#### JC: added the 'augmented' boolean argument
################################################################################
pin <- function (fullmodel, transform, augmented=F) {
  # Usage example: FM.asr, herit.and.se ~ V2/(V1+V10),
  # with herit.and.se whatever row name we want to give to the return result table
  # fullmodel is an asreml object
  # Transform is a model of the type: herit.and.se ~ V2/(V1+V10)
  ### From http://www.homepages.ed.ac.uk/iwhite//asreml/uop
  # Need this function to calculate the standard errors of heritability;
  # It also produces heritability estimates.
  # This estimate and the one below (based on Holland et al.) should agree but the heritability
  pframe <- as.list(fullmodel$gammas)
  ################################################################################
  #### JC: if using augmented design, need to correct pattern-matching to identify
  #### the correct variance components
  ################################################################################
  if (augmented) {
  # gammaterms.with.G <- grep( paste('^', VCterms[1],':IS_EXPERIMENTAL!', VCterms[1], sep=''), names(pframe)) # MODIFIED FOR AUGMENTED DESIGN
    gammaterms.with.G <- grep(VCterms[1], names(pframe)) ### modified by Xiaowei LI
  } else {
    gammaterms.with.G <- grep( paste('^', VCterms[1],'!', VCterms[1], sep=''), names(pframe)) # ORIGINAL CODE: USE FOR STANDARD DESIGN
  }
  ################################################################################
  ### JC: End modification
  ################################################################################
  names(pframe)[gammaterms.with.G] <- gsub("[:]","_", VCterms[1])
  gammaterms.residuals <- grep('variance', names(pframe))
  names(pframe)[gammaterms.residuals] <- 'Residuals'
  ### loop over the 2-ways and three-ways to assign names
  
  for (vct in VCterms) {  ## VCterms is a vector with (Genocode, Env)
    vctcol = unlist(strsplit(vct, ':'))
    if (length(vctcol) > 2 ) { ### >2 modified by Xiaowei Li 
      gammaterms.with.GE <- grep( paste('^',vct,'!',sep='' ), names(pframe))
      vct <- gsub("[:]","_", vct)
      names(pframe)[gammaterms.with.GE] <- vct
    }
  }
  
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[fullmodel$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3)
    transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- fullmodel$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}


################################################################################
#### JC: added the 'augmented' boolean argument
################################################################################
calc.herit <- function (fullmodel, curr.trait, all.harm.mean, augmented=F){  # all_harm.mean is a vector with all the harmonic means calculated above
  # e.g. all_harm.mean is c(Entry_term, Entry_Trt_term, Entry_Year_term, Entry_Trt_Year_term)
  #
  # http://www4.ncsu.edu/~jholland/heritability/MultiEnvironLatticeHeritability.sas
  # harmonic.mean (cleanedpheno, curr.trait, VCterms) -> arm
  # arm$line -> harm.mean.line
  # arm$env -> harm.mean.GbyE
  # get variance components estimates and the variance-covariance matrix
  # fullmodel <- FM.asr
  covmat <- asreml:::asreml.ltri2mat(fullmodel$ai)
  var.comps <- summary(fullmodel)$varcomp
  # put variance components into a vector
  VC <- var.comps[,2]
  VC <- matrix(VC, nrow=1)
  colnames(VC) <- row.names(var.comps)
  
  # Rename the var comp terms in the asREML model output (the VC matrix) so they match the terms specified by the user (e.g. 'Genocode!Genocode.var' becomes 'Genocode')
  # Also subset VC to keep only the terms we will use to calculate heritability
  terms.to.keep <- vector()
  colnames.to.keep <- vector()
  # genetic variance term
  ################################################################################
  #### JC: same correction as previous function
  ################################################################################
  if (augmented) {
  # varcompterms.with.G <- grep( paste('^', VCterms[1],':IS_EXPERIMENTAL!', VCterms[1], sep=''), colnames(VC)) # MODIFIED FOR AUGMENTED DESIGN
    varcompterms.with.G <- grep(VCterms[1], colnames(VC)) ### modified by Xiaowei LI
  } else {
    varcompterms.with.G <- grep( paste('^', VCterms[1],'!', VCterms[1], sep=''), colnames(VC)) # ORIGINAL CODE: USE FOR STANDARD DESIGN
  }
  ################################################################################
  ### JC: End modification
  ################################################################################
  terms.to.keep <- c(terms.to.keep, varcompterms.with.G)
  colnames.to.keep <- c(colnames.to.keep, gsub("[:]","_", VCterms[1]))    ### modified by Xiaowei Li
  # residual term
  varcompterms.residuals <- grep('variance', colnames(VC))
  terms.to.keep <- c(terms.to.keep, varcompterms.residuals)
  colnames.to.keep <- c(colnames.to.keep, 'Residuals')
  # loop over the 2-ways and three-ways to assign names
  for (vct in VCterms) {  ## VCterms is a vector with (Genocode, Env)
    vctcol = unlist(strsplit(vct, ':'))
    if (length(vctcol) > 2 ) {
      varcompterms.with.GE <- grep( paste('^',vct,'!',sep='' ), colnames(VC))
      terms.to.keep <- c(terms.to.keep, varcompterms.with.GE)
      vct <- gsub("[:]","_", vct)
      colnames.to.keep <- c(colnames.to.keep, vct)
    }
  }
  
  VC <- VC[terms.to.keep]
  VC <- matrix(VC, nrow=1)
  colnames(VC) <- colnames.to.keep
  
  
  ###########################################
  ######### Calculate plot level heritability
  ###########################################
  # VC terms to extract
  
  Terms.with.G <- rep( 0, length(colnames(VC)) )
  Terms.with.G[which(colnames(VC) == gsub("[:]","_", VCterms[1]))] <- 1   ### modified by Xiaowei Li
  
  VG.p <- Terms.with.G * VC
  VP.p <- VC
  
  # Estimate plot level heritability
  Heritability.plot.basis <- sum(VG.p)/sum(VC)
  
  # Calculate SE of heritability on plot level basis with the pin function
  collapsedterms <- paste(colnames(VC), collapse=' + ')
  transformModel <- paste(gsub("[:]","_", VCterms[1]), ' / (' , collapsedterms , ')')   ### modified by Xiaowei Li
  
  ################################################################################
  ### JC: Modification to pass the argument for 'augmented' parameter
  ################################################################################
  se.heritability.pin.plot <- eval(parse(text = paste('pin (fullmodel, h2.plot ~ ', transformModel, ',',augmented,')', sep=''))) # modification
  ################################################################################
  ### JC: End modification
  ################################################################################
  
  # Extract value
  se.pin.plot <- se.heritability.pin.plot[1,2]
  h2.pin.plot <- se.heritability.pin.plot[1,1]
  
  # ############################################
  # ######### Calculate family mean heritability; here the var comp are divided by the harmonic mean
  # ############################################
  # # Use same TermsWithG as above, genetic variance is not devided by any harmonic mean
  # VG.f <- VG.p
  # ### All the other terms are divided by their respective harmonic mean
  # ### And the residuals are divided by the harmonic mean of all genotypes (harm.mean.line)
  # all.VC.over.HM <- vector()  # all.harm.mean <- harm.means
  # for (i in colnames(all.harm.mean)) {
  #   # i <- "CHECK_GENOTYPE" 
  #   # i <-  "CHECK_YEAR_GENOTYPE_IS_EXPERIMENTAL"
  #   
  #   VC.over.HM <- VC[which(colnames(VC) == i) ]/ all.harm.mean[which(colnames(all.harm.mean) == i)]
  #   all.VC.over.HM <- c(all.VC.over.HM, VC.over.HM)
  # }
  # all.VC.over.HM <- matrix(all.VC.over.HM, nrow=1)
  # colnames(all.VC.over.HM) <- colnames(all.harm.mean)
  # 
  # ### Estimate family level heritability
  # Heritability.family.basis <- sum(VG.f)/ ( sum(all.VC.over.HM) + sum(VG.f) )
  # 
  # # Calculate SE of heritability on family level basis with the pin function
  # collapsedterms <- c(gsub("[:]","_", VCterms[1]))   ### modified by Xiaowei Li
  # for (i in colnames(all.harm.mean)) {
  #   # i <- "CHECK_GENOTYPE" 
  #   # i <-  "CHECK_YEAR_GENOTYPE_IS_EXPERIMENTAL"
  #   
  #   n <- which(colnames(all.harm.mean) == i)
  #   tmp <- paste (i, '/', all.harm.mean[n], sep='')
  #   collapsedterms <- c(collapsedterms, tmp)
  # }
  # allcollapsedterms <- paste(collapsedterms, collapse=" + ")
  # transformModel <- paste(gsub("[:]","_", VCterms[1]), ' / (' , allcollapsedterms , ')')   ### modified by Xiaowei Li
  # ################################################################################
  # ### JC: Modification to pass the argument for 'augmented' parameter
  # ################################################################################
  # se.heritability.pin.family <- eval(parse(text = paste('pin (fullmodel, h2.family ~ ', transformModel, ',',augmented,')', sep=''))) # modification
  # ################################################################################
  # ### JC: End modification
  # ################################################################################
  # # Get value
  # se.pin.family <- se.heritability.pin.family[1,2]
  # h2.pin.family <- se.heritability.pin.family[1,1]
  return(list(h2plot=Heritability.plot.basis, h2plot.pin=h2.pin.plot, seplot.pin=se.pin.plot))
}
