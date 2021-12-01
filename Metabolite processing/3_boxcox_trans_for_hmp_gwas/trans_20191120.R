
############     remove sweet&popcon and perform BOX-COX transformation

library(data.table)

### readin line list (no sweet and popcorn)
d1462 <- read.delim('ready_for_imputation_samples.txt',header = F)

### ve 9 traits
blue_ve <- read.delim('BLUES_all_traits_input_test_with_checks_update.txt')
blue <- blue_ve
blue <- blue[which(blue$CHECK==99),]
blue <- blue[-which(blue$GENOTYPE=="14SJMY:Breed:0001"),-1]
### change Source ID for DE3
blue$GENOTYPE <- as.character(blue$GENOTYPE)
blue$GENOTYPE[which(blue$GENOTYPE=="13A089451A")] <- "PI 638551"
blue$GENOTYPE2 <- gsub(" ","",blue$GENOTYPE)
### merge
mer <- merge(blue,d1462,by.x="GENOTYPE2",by.y="V1")
mer <- mer[,-2]
colnames(mer)[1] <- "GENOTYPE"
write.csv(mer,file="VE_v6_locked_blue_raw_9_traits.csv", row.names = F)

### generate file for transformation use
ve_9 <- read.csv("VE_v6_locked_blue_raw_9_traits.csv")
mer_t <- ve_9

### write out constant
ct <- NULL

for (i in 2:10){
  m <- min(mer_t[,i],na.rm = T)
  ma <- c(colnames(mer_t)[i],m)
  ct <- rbind(ct,ma)
}

ct <- as.data.frame(ct)
colnames(ct) <- c("trait","minimum")
ct$minimum <- as.numeric(as.character(ct$minimum))
ct$constant <- ifelse(ct$minimum <= 0, abs(ct$minimum)+ 10^-9, 0)

fwrite(ct,"constant_added_for_toco_traits_before_box_cox.csv",row.names = F)

###### transformation
### adding small constant to traits with values <= 0
for (i in 2:10){
  
  if (length(which(mer_t[,i]<=0)) > 0) {
    mer_t[,i] <-  mer_t[,i]- min(mer_t[,i],na.rm = T) + 10^-9
  }
}

source('s1_boxcox_function.R')

traits <- c( "d.T3", "g.T3", "a.T3", "d.T", "g.T", "a.T", "Total.Tocotrienols", "Total.Tocopherols", "Total.Tocochromanols" )

infile <- 'input_test_with_checks.txt'

transf.dataset <- as.data.frame(mer_t$GENOTYPE)
colnames(transf.dataset) <- "GENOTYPE"
boxcox.transformed <- matrix(data=c('Trait', 'Lambda', 'Transformed'), ncol=3, nrow=1, byrow=TRUE)

shell(paste('mkdir ', "boxcox", sep=''))
setwd("boxcox")

for (i in c(1:9)) {
  
#### Current trait ####
curr.trait <- traits[i]

## STEP 4 ## Check if there is need for a BoxCox transformation
boxcox.fitting(mer_t, curr.trait) -> transformed.trait
as.matrix(transformed.trait$TRT) -> transf.trait
colnames(transf.trait) <- curr.trait

if (i==1){
  transfpheno <- cbind(transf.dataset, transf.trait)
}else{
  transfpheno <- cbind(transfpheno, transf.trait)
}

if (transformed.trait$LAMBDA == 1){
  boxcox.transformed <- rbind(boxcox.transformed, c(curr.trait, transformed.trait$LAMBDA, 'NO'))
}else{
  boxcox.transformed <- rbind(boxcox.transformed, c(curr.trait, transformed.trait$LAMBDA, 'YES'))
}

}

setwd("..\\")

write.table(transfpheno, file=paste("TOCO_v6_locked_blue_", 'transformed_tocotraits_9.txt',  sep=''),  quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
write.table(boxcox.transformed, file=paste('boxcox_transformation_tocotraits_9_applied_', infile, sep=''),  quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

