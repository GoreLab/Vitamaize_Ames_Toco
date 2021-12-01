# Args <- commandArgs(TRUE)
# ch <- Args[1]
# t <- Args[2]

# ch=1
# sink()
setwd("/workdir/dw524/vitamaize/eQTL/GAPIT/")
dir <- "/workdir/dw524/vitamaize/eQTL/GAPIT/"
geno_dir <- "/workdir/dw524/vitamaize/genotype/10.kinship_for_TWAS_v1_1/"


#                                             #GWAS#
##########################################################################################################
#Step 0: Set directory to load GAPIT source files (this step is omitted for using package)
library('MASS')
# source("http://bioconductor.org/biocLite.R")
# biocLite("multtest")
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R 
library("scatterplot3d")
#source('http://www.zzlab.net/GAPIT/previous/gapit_functions20190714.txt')
source("GAPIT_function_20190716_version.txt")
source("http://zzlab.net/GAPIT/emma.txt")
library(data.table)
library(rrBLUP)
library(matrixcalc)

# source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")

myY <- fread(paste(dir,"PeerResiduals_RmOut_v1.1_B73.txt",sep=''), header = TRUE,data.table = F)
myKI<- read.csv(paste(dir, 'GAPIT.Kin.VanRaden.csv', sep = ''), header = F)
names(myKI)[2:ncol(myKI)]=as.character(myKI[,1])

finished=list.files(path ='ch10/', pattern='*GWAS.Results.csv')
finished=gsub('.GWAS.Results.csv','',finished)
finished=gsub('GAPIT.MLM.','',finished)

genes=names(myY)[-1]
genes=genes[-which(genes %in% finished)]

myY$ID=gsub('_','',myY$ID)
myY=myY[,c('ID',genes)]

ind=which(myKI[,1] %in% myY$ID)
myKI_2=myKI[ind,c(1,ind+1)]
myY=myY[myY$ID %in% myKI_2[,1],]

# run eQTL -------------------------------------------------------------------
try(system(paste("mkdir -p ", dir,"ch",ch, sep='')))
setwd(paste(dir,"ch",ch, sep=''))
curr_chr=fread(paste(geno_dir, 'AGPv4_Ames_vitamaize', ch, '_exp_subset.012', sep = ''),head=F,data.table = F)
indv=read.delim(paste(geno_dir, 'AGPv4_Ames_vitamaize', ch, '_exp_subset.012.indv', sep = ''),head=F)
map=read.delim(paste(geno_dir, 'AGPv4_Ames_vitamaize', ch, '_exp_subset.012.pos', sep = ''),head=F)
map$ID=paste(map$V1,map$V2,sep='_')

myGD=curr_chr
names(myGD)=c('taxa',map$ID)
myGD$taxa=indv$V1
myGD$taxa=as.character(myGD$taxa)

myGM=map[,c(3,1,2)]
names(myGM)=c('Names','Chromosome','Position')

myGD[1:5,1:5]
myGM[1:5,1:3]

myGM_sub=myGM
myGD_sub=myGD[myGD$taxa %in% myY$ID,]


log_file <- paste(dir,"ch",ch,'/SAM_log_GAPIT_',ch,'_all.txt',  sep = '')
sink(log_file)

myGAPIT <- GAPIT(Y=myY,#This is phenotype data
                 GD=myGD_sub,
                 GM=myGM_sub,
                 KI=myKI_2,
                 group.from=nrow(myY),		#Lower bound for number of group
                 group.to=nrow(myY),			#Upper bound for number of group
                 group.by=1,				#range between 1 and number of individuals, smaller the finer optimization
                 SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
                 SNP.impute = "Middle",
                 cutOff = 0.05,
                 Major.allele.zero = TRUE,
                 file.output = T
)


rm(list=ls())
