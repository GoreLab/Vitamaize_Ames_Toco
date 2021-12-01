### check the transformation mannually (2020-6-26)

### generate file for transformation use
ve_9 <- read.csv("VE_v6_locked_blue_raw_9_traits.csv")
mer_t <- ve_9

###### transformation
### adding small constant to traits with values <= 0
for (i in 2:10){
  
  if (length(which(mer_t[,i]<=0)) > 0) {
    mer_t[,i] <-  mer_t[,i]- min(mer_t[,i],na.rm = T) + 10^-9
  }
}

r_ra <- merge(ve_9,mer_t,by="GENOTYPE")
plot(r_ra$Total.Tocotrienols.x ~ r_ra$Total.Tocotrienols.y)
fwrite(r_ra,"ve_raw_blue_and_adding_samll_value_check.csv",row.names = F)

# if(lambda.to.use == -2) {
#   newvalue <- 1/(pheno[curr.trait]^2)
# } else if (lambda.to.use == -1) {
#   newvalue <- 1/(pheno[curr.trait])
# } else if( lambda.to.use == -0.5) {
#   newvalue <- 1/(sqrt(pheno[curr.trait]))
# } else if (lambda.to.use == 0) {
#   newvalue <- log(pheno[curr.trait])
# } else if (lambda.to.use == 0.5) {
#   newvalue <- sqrt(pheno[curr.trait])
#   # this is for the case where 1 is not in the CI but it is the closest convenient value
# } else if (lambda.to.use == 1) {
#   newvalue <- pheno[curr.trait]
# } else { # (lambda.to.use == 2)
#   newvalue <- pheno[curr.trait]^2
# }
# 
# 
# Trait	Lambda	Transformed
# d.T3	0	YES
# g.T3	0.5	YES
# a.T3	0.5	YES
# d.T	0.5	YES
# g.T	0.5	YES
# a.T	0.5	YES
# Total.Tocotrienols	0	YES
# Total.Tocopherols	0.5	YES
# Total.Tocochromanols	0.5	YES

mer_t$d.T3_m <- log(mer_t$d.T3)
plot(mer_t$d.T3_m~mer_t$d.T3)

mer_t$g.T3_m <-  sqrt(mer_t$g.T3)
plot(mer_t$g.T3_m~mer_t$g.T3)

mer_t$a.T3_m <-  sqrt(mer_t$a.T3)
mer_t$d.T_m <- sqrt(mer_t$d.T)
mer_t$g.T_m <-  sqrt(mer_t$g.T)
mer_t$a.T_m <-  sqrt(mer_t$a.T)
mer_t$Total.Tocotrienols_m <- log(mer_t$Total.Tocotrienols)
mer_t$Total.Tocopherols_m <- sqrt(mer_t$Total.Tocopherols)
mer_t$Total.Tocochromanols_m <- sqrt(mer_t$Total.Tocochromanols)
mer_t <- mer_t[,-c(2:10)]

library(data.table)

tb <- fread("TOCO_v6_locked_blue_transformed_tocotraits_9.txt",data.table = F)
mer <- merge(tb,mer_t,by="GENOTYPE")
c <- mer[mer$d.T3==mer$d.T3_m,]
c <- mer[mer$g.T3==mer$g.T3_m,]

plot(mer$d.T3~mer$d.T3_m)
plot(mer$g.T3~mer$g.T3_m)
plot(mer$a.T3~mer$a.T3_m)
plot(mer$d.T~mer$d.T_m)
plot(mer$g.T~mer$g.T_m)
plot(mer$a.T~mer$a.T_m)
plot(mer$Total.Tocotrienols~mer$Total.Tocotrienols_m)
plot(mer$Total.Tocopherols~mer$Total.Tocopherols_m)
plot(mer$Total.Tocochromanols~mer$Total.Tocochromanols_m)

fwrite(mer,"ve_blue_check.csv",row.names = F)

