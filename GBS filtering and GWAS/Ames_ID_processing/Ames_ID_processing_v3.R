### This Ames ID processing script (by Xiaowei Li) contains following analyses:
# (1) compare Pedigree info from Maria end and that from Entry_coder (Laura)
# (2) check "location_to_Source_15-17_key.csv" file, this file was from Laura 2019-6-20
# (3) compare Pedigree info from Maria end and that from key file
# (4) connect lines from Maria end to source ID from Laura end
# (5) generate Ames ID full list for GBS extraction
# (6) retrive full list of GBS records from Panzea_ZeaGBSv2.7
# (7) remove sweet corn and popcorn
# (8) compare sweet corn and popcorn list with Laura's
# (9) remove sweet corn and popcorn from full list
# (10) retrive common line list of GBS records from Panzea_ZeaGBSv2.7
# (11) remove "ae" and seven accessions to generate locked list
# (12) retrive common line list of GBS records from Panzea_ZeaGBSv2.7

#################################     (1) compare Pedigree info from Maria end and that from Entry_coder (Laura)    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

### check entry_coder file, this file was compiled by Xiaowei Li using files from Laura
coder <- read.csv('entry_coder_isu.csv',header = TRUE)
coder_ex <- coder[which(coder$Pedigree!= "B73"),]
coder_p_dup <- coder_ex[duplicated2(coder_ex$Pedigree),] # 12 pedigree_dup
coder_s_dup <- coder_ex[duplicated2(coder_ex$Source),] # zero source_dup
write.csv(coder_p_dup ,file="coder_p_dup.csv", row.names = F)

### d2015_vb
d2015_vb <- read.csv('Ames 2015 - BVits Final_021518.csv',header = TRUE)  # 1802 lines
d2015_vb <- d2015_vb[,c(6,8,9)]
d2015_vb$Range_Pass_2015_vb <- paste(d2015_vb$Range,d2015_vb$Pass,sep="_")
d2015_vb_c <- merge(d2015_vb,coder,by.x="Range_Pass_2015_vb",by.y="Range_Pass_2015")
d2015_vb_c_dif <- d2015_vb_c[which(as.character(d2015_vb_c$Pedigree.x)!=as.character(d2015_vb_c$Pedigree.y)),]

### d2015_ve
d2015_ve <- read.csv('AMES15 - TOCOCHROMANOLS - REVISED_022219.csv',header = TRUE) # 1801 lines
d2015_ve <- d2015_ve[,c(6,8,9)]
d2015_ve$Range_Pass_2015_ve <- paste(d2015_ve$Range,d2015_ve$Pass,sep="_")
d2015_ve_c <- merge(d2015_ve,coder,by.x="Range_Pass_2015_ve",by.y="Range_Pass_2015")
d2015_ve_c_dif <- d2015_ve_c[which(as.character(d2015_ve_c$Pedigree.x)!=as.character(d2015_ve_c$Pedigree.y)),]

d2015_dif <-merge(d2015_vb_c_dif,d2015_ve_c_dif,by.x="Range_Pass_2015_vb",by.y="Range_Pass_2015_ve",all=T)
write.csv(d2015_dif ,file="d2015_dif.csv", row.names = F)

### d2017_vb
d2017_vb <- read.csv('Ames 2017_BVits_FINAL_d3Thiamine corrected and plate info corrected_022219.csv', header = TRUE) # 1748 lines
d2017_vb <- d2017_vb[,c(6,8,9)]
d2017_vb$Range_Pass_2017_vb <- paste(d2017_vb$Range,d2017_vb$Pass,sep="_")
d2017_vb_c <- merge(d2017_vb,coder,by.x="Range_Pass_2017_vb",by.y="Range_Pass_2017")
d2017_vb_c_dif <- d2017_vb_c[which(as.character(d2017_vb_c$Pedigree.x)!=as.character(d2017_vb_c$Pedigree.y)),]

### d2017_ve
d2017_ve <- read.csv('AMES 2017_TOCOCHROMANOLS_FINAL_.csv', header = TRUE) # 1738 lines
d2017_ve <- d2017_ve[,c(6,8,9)]
d2017_ve$Range_Pass_2017_ve <- paste(d2017_ve$Range,d2017_ve$Pass,sep="_")
d2017_ve_c <- merge(d2017_ve,coder,by.x="Range_Pass_2017_ve",by.y="Range_Pass_2017")
d2017_ve_c_dif <- d2017_ve_c[which(as.character(d2017_ve_c$Pedigree.x)!=as.character(d2017_ve_c$Pedigree.y)),]

d2017_dif <-merge(d2017_vb_c_dif,d2017_ve_c_dif,by.x="Range_Pass_2017_vb",by.y="Range_Pass_2017_ve",all=T)
write.csv(d2017_dif ,file="d2017_dif.csv", row.names = F)
### compare Pedigree info, finished!

#################################     (2) check "location_to_Source_15-17_key.csv" file, this file was from Laura 2019-6-20    ########################################

### read in raw file
key <- read.csv('location_to_Source_15-17_key.csv',header = TRUE) # 3790 rows
### add "id" column
key$id <- paste(key$Year,key$Range,key$Pass,sep="_")
key_dup <- key[duplicated2(key$id),] # 10 rows with duplications
write.csv(key_dup ,file="key_dup.csv", row.names = F)
### remove 10 rows dups
key_u <- distinct(key, id, .keep_all = TRUE) # 3790 - 10 = 3780 rows
write.csv(key_u ,file="location_to_Source_15-17_key_update.csv", row.names = F)

### read in updated key file
key <- read.csv('location_to_Source_15-17_key_update.csv',header = TRUE) # 3780 rows

### check 2015 year data
key_15 <- key[which(key$Year== "2015"),]
key_15_ex <- key_15[which(key_15$Pedigree!= "B73"),]
key_15_ex_p_dup <- key_15_ex[duplicated2(key_15_ex$Pedigree),] # 12 pedigree_dup
key_15_ex_s_dup <- key_15_ex[duplicated2(key_15_ex$Source),] # zero source_dup

### check 2017 year data
key_17 <- key[which(key$Year== "2017"),]
key_17_ex <- key_17[which(key_17$Pedigree!= "B73"),]
key_17_ex_p_dup <- key_17_ex[duplicated2(key_17_ex$Pedigree),] # 12 pedigree_dup
key_17_ex_s_dup <- key_15_ex[duplicated2(key_17_ex$Source),] # zero source_dup
### the key file is fine! finished!

#################################     (3) compare Pedigree info from Maria end and that from key file    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

### read in updated key file
key <- read.csv('location_to_Source_15-17_key_update.csv',header = TRUE) # 3780 rows

### d2015_vb
d2015_vb <- read.csv('Ames 2015 - BVits Final_021518.csv',header = TRUE)  # 1802 lines
d2015_vb <- d2015_vb[,c(6,8,9)]
d2015_vb$id_2015_vb <- paste("2015",d2015_vb$Range,d2015_vb$Pass,sep="_")
d2015_vb_k <- merge(d2015_vb,key,by.x="id_2015_vb",by.y="id")
d2015_vb_k_dif <- d2015_vb_k[which(as.character(d2015_vb_k$Pedigree.x)!=as.character(d2015_vb_k$Pedigree.y)),]
### the d2015_vb_k_dif discrepancy from format and spelling errors and all of them can be accepted!

### d2015_ve
d2015_ve <- read.csv('AMES15 - TOCOCHROMANOLS - REVISED_022219.csv',header = TRUE) # 1801 lines
d2015_ve <- d2015_ve[,c(6,8,9)]
d2015_ve$id_2015_ve <- paste("2015",d2015_ve$Range,d2015_ve$Pass,sep="_")
d2015_ve_k <- merge(d2015_ve,key,by.x="id_2015_ve",by.y="id")
d2015_ve_k_dif <- d2015_ve_k[which(as.character(d2015_ve_k$Pedigree.x)!=as.character(d2015_ve_k$Pedigree.y)),]
### the d2015_ve_k_dif discrepancy from format and spelling errors and all them can be accepted!

### d2017_vb
d2017_vb <- read.csv('Ames 2017_BVits_FINAL_d3Thiamine corrected and plate info corrected_022219.csv', header = TRUE) # 1748 lines
d2017_vb <- d2017_vb[,c(6,8,9)]
d2017_vb$id_2017_vb <- paste("2017",d2017_vb$Range,d2017_vb$Pass,sep="_")
d2017_vb_k <- merge(d2017_vb,key,by.x="id_2017_vb",by.y="id")
d2017_vb_k_dif <- d2017_vb_k[which(as.character(d2017_vb_k$Pedigree.x)!=as.character(d2017_vb_k$Pedigree.y)),]

### This planting error means that B73 was really planted in Range 5, Pass 74, Tier 4 in 2017; 
### H22w was really planted in Range 5 Pass 75 Tier 4 that year. (comment from Laura)
### all d2017_vb_k_dif discrepancies can be accepted!

### d2017_ve
d2017_ve <- read.csv('AMES 2017_TOCOCHROMANOLS_FINAL_.csv', header = TRUE) # 1738 lines
d2017_ve <- d2017_ve[,c(6,8,9)]
d2017_ve$id_2017_ve <- paste("2017",d2017_ve$Range,d2017_ve$Pass,sep="_")
d2017_ve_k <- merge(d2017_ve,key,by.x="id_2017_ve",by.y="id")
d2017_ve_k_dif <- d2017_ve_k[which(as.character(d2017_ve_k$Pedigree.x)!=as.character(d2017_ve_k$Pedigree.y)),]

### This planting error means that B73 was really planted in Range 5, Pass 74, Tier 4 in 2017; 
### H22w was really planted in Range 5 Pass 75 Tier 4 that year. (comment from Laura)
### all d2017_ve_k_dif discrepancies can be accepted!
### the pedigree column in key file is fine and can be used for analysis!

#################################     (4) connect lines from Maria end to source ID from Laura end    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

### read in updated key file
key <- read.csv('location_to_Source_15-17_key_update.csv',header = TRUE) # 3780 rows

### d2015_vb
d2015_vb <- read.csv('Ames 2015 - BVits Final_021518.csv',header = TRUE)  # 1802 lines
d2015_vb$id_2015_vb <- paste("2015",d2015_vb$Range,d2015_vb$Pass,sep="_")
d2015_vb_k <- merge(d2015_vb,key,by.x="id_2015_vb",by.y="id")
d2015_vb_k2 <- d2015_vb_k[,c(2:6,31:33,27:28,8:26)]
colnames(d2015_vb_k2)[5] <- "Source"
colnames(d2015_vb_k2)[7] <- "Pedigree"
colnames(d2015_vb_k2)[8] <- "ID"
colnames(d2015_vb_k2)[12] <- "Range"
colnames(d2015_vb_k2)[13] <- "Pass"
### check whether each line was measured one time
d2015_vb_k2_dup <- d2015_vb_k2[duplicated2(d2015_vb_k2$PedId),] # yes
write.csv(d2015_vb_k2 ,file="Ames_2015_vb_adding_ID.csv", row.names = F)
### Ames_2015_vb_adding_ID file finished!

### d2015_ve
d2015_ve <- read.csv('AMES15 - TOCOCHROMANOLS - REVISED_022219.csv',header = TRUE) # 1801 lines
d2015_ve$id_2015_ve <- paste("2015",d2015_ve$Range,d2015_ve$Pass,sep="_")
d2015_ve_k <- merge(d2015_ve,key,by.x="id_2015_ve",by.y="id")
d2015_ve_k2 <- d2015_ve_k[,c(2:6,31:33,27:28,8:26)]
colnames(d2015_ve_k2)[5] <- "Source"
colnames(d2015_ve_k2)[7] <- "Pedigree"
colnames(d2015_ve_k2)[8] <- "ID"
colnames(d2015_ve_k2)[12] <- "Range"
colnames(d2015_ve_k2)[13] <- "Pass"
### check whether each line was measured one time
d2015_ve_k2_dup <- d2015_ve_k2[duplicated2(d2015_ve_k2$PedId),] # yes
write.csv(d2015_ve_k2 ,file="Ames_2015_ve_adding_ID.csv", row.names = F)
### Ames_2015_ve_adding_ID file finished!

### d2017_vb
d2017_vb <- read.csv('Ames 2017_BVits_FINAL_d3Thiamine corrected and plate info corrected_022219.csv', header = TRUE) # 1748 lines
d2017_vb$id_2017_vb <- paste("2017",d2017_vb$Range,d2017_vb$Pass,sep="_")
d2017_vb_k <- merge(d2017_vb,key,by.x="id_2017_vb",by.y="id")
d2017_vb_k2 <- d2017_vb_k[,c(2:6,31:33,27:28,8:26)]
colnames(d2017_vb_k2)[5] <- "Source"
colnames(d2017_vb_k2)[7] <- "Pedigree"
colnames(d2017_vb_k2)[8] <- "ID"
colnames(d2017_vb_k2)[12] <- "Range"
colnames(d2017_vb_k2)[13] <- "Pass"
### check whether each line was measured one time
d2017_vb_k2_dup <- d2017_vb_k2[duplicated2(d2017_vb_k2$PedId),] # 10 samples were measured two times
### to be consistent, remove those 10 measurements from Plate 19
remove_id <- d2017_vb_k2_dup[which(d2017_vb_k2_dup$LCMS.Plate=="Ames17-BV-19L"),c(1:3)]
all_row <- NULL # get row number
for (i in 1:10){
  tmp_row <- which(d2017_vb_k2$LCMS.Plate=="Ames17-BV-19L" & d2017_vb_k2$PedId==remove_id$PedId[i] )
  all_row <- c(all_row,tmp_row)
}
d2017_vb_k3 <- d2017_vb_k2[-all_row,]
d2017_vb_k3_dup <- d2017_vb_k3[duplicated2(d2017_vb_k3$PedId),]

### This planting error means that B73 was really planted in Range 5, Pass 74, Tier 4 in 2017; 
### H22w was really planted in Range 5 Pass 75 Tier 4 that year. (comment from Laura)

### correct the mistakes for Event.Name column
r5p74_75 <- d2017_vb_k3[which((d2017_vb_k3$Range=="5" & d2017_vb_k3$Pass=="74")|(d2017_vb_k3$Range=="5" &d2017_vb_k3$Pass=="75")),]
d2017_vb_k3$Event.Name[which(d2017_vb_k3$Range=="5" & d2017_vb_k3$Pass=="74")] <- "Check"
d2017_vb_k3$Event.Name[which(d2017_vb_k3$Range=="5" & d2017_vb_k3$Pass=="75")] <- "0"

write.csv(d2017_vb_k3 ,file="Ames_2017_vb_adding_ID.csv", row.names = F)
### Ames_2017_vb_adding_ID file finished!

### d2017_ve
d2017_ve <- read.csv('AMES 2017_TOCOCHROMANOLS_FINAL_.csv', header = TRUE) # 1738 lines
d2017_ve$id_2017_ve <- paste("2017",d2017_ve$Range,d2017_ve$Pass,sep="_")
d2017_ve_k <- merge(d2017_ve,key,by.x="id_2017_ve",by.y="id")
d2017_ve_k2 <- d2017_ve_k[,c(2:6,30:32,26:27,8:25)]
colnames(d2017_ve_k2)[5] <- "Source"
colnames(d2017_ve_k2)[7] <- "Pedigree"
colnames(d2017_ve_k2)[8] <- "ID"
colnames(d2017_ve_k2)[12] <- "Range"
colnames(d2017_ve_k2)[13] <- "Pass"
### check whether each line was measured one time
d2017_ve_k2_dup <- d2017_ve_k2[duplicated2(d2017_ve_k2$PedId),] # all samples were measured one time

### This planting error means that B73 was really planted in Range 5, Pass 74, Tier 4 in 2017; 
### H22w was really planted in Range 5 Pass 75 Tier 4 that year. (comment from Laura)

### correct the mistakes for Event.Name column
r5p74_75 <- d2017_ve_k2[which((d2017_ve_k2$Range=="5" & d2017_ve_k2$Pass=="74")|(d2017_ve_k2$Range=="5" &d2017_ve_k2$Pass=="75")),]
d2017_ve_k2$Event.Name[which(d2017_ve_k2$Range=="5" & d2017_ve_k2$Pass=="74")] <- "Check"
d2017_ve_k2$Event.Name[which(d2017_ve_k2$Range=="5" & d2017_ve_k2$Pass=="75")] <- "0"

write.csv(d2017_ve_k2 ,file="Ames_2017_ve_adding_ID.csv", row.names = F)
### Ames_2017_ve_adding_ID file finished!
### connecting lines from Maria end to source ID from Laura end finished!

#################################     (5) generate Ames ID full list for GBS extraction    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

###### read in 2015 raw table
d2015_vb <- read.csv('Ames_2015_vb_adding_ID.csv',header = TRUE)[,c(7:9,14)]  # 1802 lines
d2015_ve <- read.csv('Ames_2015_ve_adding_ID.csv',header = TRUE)[,c(7:9,14)] # 1801 lines
### exclude control B73
d2015_vb_ex <- d2015_vb[-which(d2015_vb$Pedigree=="B73"),] # 1802 - 103 = 1699 lines
d2015_ve_ex <- d2015_ve[-which(d2015_ve$Pedigree=="B73"),] # 1801 - 103 = 1698 lines
### identify local control
d2015_vb_ex_dup <- d2015_vb_ex[duplicated2(d2015_vb_ex$ID),] # zero local control
d2015_ve_ex_dup <- d2015_ve_ex[duplicated2(d2015_ve_ex$ID),] # zero local control
d2015 <- merge(d2015_vb_ex, d2015_ve_ex,by='ID',all = T) # 1699 unique IDs
d2015[which(is.na(d2015$Comments.x)|is.na(d2015$Comments.y)),] # 1 (CML 154Q)
### test whether the Comments are identical
comment_diff <- d2015[which(d2015$Comments.x!=d2015$Comments.y),] # comments are the same
d2015 <- d2015_vb_ex
colnames(d2015)[4] <- "Comments_d2015"
write.csv(d2015,file="d2015_id.csv", row.names = F)

###### read in 2017 raw table
d2017_vb <- read.csv('Ames_2017_vb_adding_ID.csv', header = TRUE)[,c(7:9,14)] # 1738 lines
d2017_ve <- read.csv('Ames_2017_ve_adding_ID.csv', header = TRUE)[,c(7:9,14)] # 1738 lines
### exclude control B73
d2017_vb_ex <- d2017_vb[-which(d2017_vb$Pedigree=="B73"),] # 1738 - 44 = 1694 lines
d2017_ve_ex <- d2017_ve[-which(d2017_ve$Pedigree=="B73"),] # 1738 - 44 = 1694 lines
### identify local control
d2017_vb_ex_dup <- d2017_vb_ex[duplicated2(d2017_vb_ex$ID),] # zero (local control)
d2017_ve_ex_dup <- d2017_ve_ex[duplicated2(d2017_ve_ex$ID),] # zero (local control)
d2017 <- merge(d2017_vb_ex, d2017_ve_ex,by='ID',all = T) # 1694 unique IDs
### test whether the Comments are identical
comment_diff <- d2017[which(d2017$Comments.x!=d2017$Comments.y),] # comments are the same
d2017 <- d2017_vb_ex
colnames(d2017)[4] <- "Comments_d2017"
write.csv(d2017,file="d2017_id.csv", row.names = F)

###### merge 2015 and 2017;identify discrepancy; generate ames_isu
d2015 <- read.csv('d2015_id.csv',header = TRUE)
d2017 <- read.csv('d2017_id.csv',header = TRUE)
ames_isu <- merge(d2015, d2017, by='ID',all = T)
d2015_d2017_diff <- ames_isu[which(is.na(ames_isu$Pedigree.x)|is.na(ames_isu$Pedigree.y)),]
### using clean format to check the diff
d2015_d2017_diff$ID_cleanformat <- d2015_d2017_diff$ID
d2015_d2017_diff$ID_cleanformat <- gsub(' ','',d2015_d2017_diff$ID_cleanformat)
d2015_d2017_diff$ID_cleanformat <- tolower(d2015_d2017_diff$ID_cleanformat)
dup <- d2015_d2017_diff[duplicated2(d2015_d2017_diff$ID_cleanformat),] # the diff is true
### merge pedigree column
ames_isu$Pedigree <- as.character(ames_isu$Pedigree.x)

for (i in 1:length(ames_isu$ID)){
  if (is.na(ames_isu$Pedigree.x[i])){
    ames_isu$Pedigree[i] <- as.character(ames_isu$Pedigree.y[i])
  }
}

### merge comment column
ames_isu$d2015_d2017 <- paste(ames_isu$Comments_d2015,ames_isu$Comments_d2017,sep="_")
table(ames_isu$d2015_d2017)
ames_isu$Comments <- "other"
ames_isu$Comments[which(ames_isu$d2015_d2017=='-_popcorn')] <- "popcorn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='-_sweet corn')] <- "sweet corn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='NA_popcorn')] <- "popcorn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='NA_sweet corn')] <- "sweet corn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='popcorn_NA')] <- "popcorn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='popcorn_popcorn')] <- "popcorn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='popcorn_sweet corn')] <- "popcorn or sweet corn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='sweet corn_NA')] <- "sweet corn"
ames_isu$Comments[which(ames_isu$d2015_d2017=='sweet corn_sweet corn')] <- "sweet corn"
ames_isu <- ames_isu[,c(1,8,10)] # 1762 unique ID
dup_p <- ames_isu[duplicated2(ames_isu$Pedigree),]
write.csv(dup_p,file="ames_isu_pedigree_dup.csv", row.names = F)
write.csv(ames_isu,file="ames_isu_id.csv", row.names = F)
### ames_id full list finished!

#################################     (6) retrive full list of GBS records from Panzea_ZeaGBSv2.7    ########################################
### also note that DE3's source is known as both "13A089451A" and "PI638551" in different files (from Laura)

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

ames_isu <- read.csv(file="ames_isu_id.csv",stringsAsFactors=FALSE)
# change ID for DE3 ("13A089451A" was not available in v2.7)
ames_isu$ID[which(ames_isu$ID=="13A089451A")] <- as.character("PI 638551")
v2.7 <- read.csv(file="Panzea_ZeaGBSv2.7_id.csv")
ames_isu$ID_clean <- ames_isu$ID
ames_isu$ID_clean <- gsub(' ','',ames_isu$ID_clean)

ames_isu$ID_clean <- tolower(ames_isu$ID_clean)
v2.7$acc_v2.7_clean <- v2.7$acc_v2.7
v2.7$acc_v2.7_clean <- gsub(' ','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('-','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('_','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\(','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\)','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- tolower(v2.7$acc_v2.7_clean)

ames_isu_v2.7 <- merge(ames_isu, v2.7, by.x="ID_clean", by.y="acc_v2.7_clean",all = T)
ames_isu_v2.7 <- ames_isu_v2.7[which(!is.na(ames_isu_v2.7$Pedigree)),]
ames_isu_v2.7_na <- ames_isu_v2.7[which(is.na(ames_isu_v2.7$INDV)),] # 4 IDs without GBS records
ames_isu_v2.7_uniq <- distinct(ames_isu_v2.7, ID, .keep_all = TRUE)
ames_isu_v2.7_all <- ames_isu_v2.7[which(!is.na(ames_isu_v2.7$INDV)),] # 2087 GBS records for 1758 IDs
write.csv(ames_isu_v2.7_all[,-1],file="ames_isu_v2.7_all.csv", row.names = F)

ames_isu_v2.7_all_dup <- ames_isu_v2.7_all[duplicated2(ames_isu_v2.7_all$ID_clean),]
write.csv(ames_isu_v2.7_all_dup[,-1],file="ames_isu_v2.7_all_dup_599.csv", row.names = F)  # dup_599 list

### identify GBS records with pedigree dup
isu <- read.csv('ames_isu_id.csv',header = TRUE)
isu_dup <- isu[duplicated2(isu$Pedigree),]

all <- read.csv('ames_isu_v2.7_all.csv',header = TRUE)
isu_dup_all <- merge(isu_dup, all, by="ID")
isu_dup_all <- isu_dup_all[,-c(2:3)]
colnames(isu_dup_all)[2:3] <- c("Pedigree","Comments")
write.csv(isu_dup_all,file="ames_isu_v2.7_all_dup_p_26.csv", row.names = F)  # dup_26 list

#################################     (7) remove sweet corn and popcorn     ########################################

ames_2013gb <- read.csv('13059_2013_3103_MOESM1_ESM.csv', header = TRUE) # accession list from Cinta's paper
ames_isu <- read.csv('ames_isu_id.csv', header = TRUE)
### adding clean format ID for matching
ames_2013gb$ID_cleanformat <- ames_2013gb$Accesion.N
ames_2013gb$ID_cleanformat <- gsub(' ','',ames_2013gb$ID_cleanformat)
ames_2013gb$ID_cleanformat <- tolower(ames_2013gb$ID_cleanformat)

ames_isu$ID_cleanformat <- ames_isu$ID
ames_isu$ID_cleanformat <- gsub(' ','',ames_isu$ID_cleanformat)
ames_isu$ID_cleanformat <- tolower(ames_isu$ID_cleanformat)

ames_isu_gb <- merge(ames_isu,ames_2013gb,by="ID_cleanformat",all=T)
ames_isu_gb <- ames_isu_gb[which(!is.na(ames_isu_gb$ID)),]

### correct dup error from Cinta's list
ames_isu_gb_dup <- ames_isu_gb[duplicated2(ames_isu_gb$ID_cleanformat),]
ames_isu_gb <- ames_isu_gb[-which(ames_isu_gb$ID_cleanformat == "ames27101" & ames_isu_gb$Breeding.program == "Ontario"),]
ames_isu_gb <- ames_isu_gb[-which(ames_isu_gb$ID_cleanformat == "pi543850" & ames_isu_gb$Breeding.program == "Other"),]
ames_isu_gb <- ames_isu_gb[,-c(5:11)]
ames_isu_gb$Comments_merged <- paste(ames_isu_gb$Comments,ames_isu_gb$Pop.structure,sep="_")
table(ames_isu_gb$Comments_merged)
ames_isu_gb$Comments_isu_gb <- "other"
ames_isu_gb$Comments_note <- "ok"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "popcorn or sweet corn_unclassified")] <- "popcorn or sweet corn_unclassified"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "popcorn or sweet corn_unclassified")] <- "in_question"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "popcorn_popcorn")] <- "popcorn"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "popcorn_unclassified")] <- "popcorn_unclassified"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "popcorn_unclassified")] <- "in_question"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "sweet corn_NA")] <- "sweet corn_NA"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "sweet corn_NA")] <- "in_question"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "sweet corn_non-stiff stalk")] <- "sweet corn_non-stiff stalk"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "sweet corn_non-stiff stalk")] <- "in_question"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "sweet corn_stiff stalk")] <- "sweet corn_stiff stalk"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "sweet corn_stiff stalk")] <- "in_question"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "sweet corn_sweet corn")] <- "sweet corn"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "sweet corn_unclassified")] <- "sweet corn_unclassified"
ames_isu_gb$Comments_note[which(ames_isu_gb$Comments_merged == "sweet corn_unclassified")] <- "in_question"

ames_isu_gb <- ames_isu_gb[,-c(4:6)]
in_ques_38 <- ames_isu_gb[which(ames_isu_gb$Comments_note=="in_question"),] ### 38 accessions
write.csv(in_ques_38,file="ames_isu_in_question_38.csv", row.names = F)

### check these accessions in NPGS and fixed the discrepancies
npgs_38 <- read.csv('ames_isu_in_question_38_npgs.csv', header = TRUE,stringsAsFactors=FALSE) ### 36 of the 38 in_question were fixed
ames_isu_gb_ok <- ames_isu_gb[which(ames_isu_gb$Comments_note=="ok"),]

for (i in c(1:38)) {
  if (npgs_38$NPGS[i]=="sweet corn"|npgs_38$NPGS[i]=="popcorn") { 
    npgs_38$Comments_note[i] <- "ok" 
    npgs_38$Comments_isu_gb[i] <- npgs_38$NPGS[i]
  }
  
  if (npgs_38$Comments_note[i]=="in_question") { 
    npgs_38$Comments_note[i] <- paste(npgs_38$Comments_note[i],"_NPGS_", npgs_38$NPGS[i],sep="")
  }
}

ames_isu_gb_npgs38 <- rbind(ames_isu_gb_ok,npgs_38[,c(1:5)])
write.csv(ames_isu_gb_npgs38,file="ames_isu_gb_npgs38.csv", row.names = F)

### combine info from Di's list
ames_isu_gb_npgs38 <- read.csv('ames_isu_gb_npgs38.csv', header = TRUE)
check_83 <- read.csv('lines_to_check_sweet_pop_vitamaize_master.csv', header = TRUE)
check_83$ID_cleanformat <- gsub(' ','',check_83$Accession)
check_83$ID_cleanformat <- tolower(check_83$Accession)
check_83 <- check_83[,c(1,9,3,8)]
mer <- merge(ames_isu_gb_npgs38,check_83,by="ID_cleanformat",all=T)
mer <- mer[which(!is.na(mer$ID)),]
table(mer$GRIN)

for (i in c(1:1762)) {
  if (!is.na(mer$GRIN[i]) & mer$GRIN[i]=="popcorn") { 
    mer$Comments_isu_gb[i] <- mer$GRIN[i]
    mer$Comments_note[i] <- "ok"
  }
  if (!is.na(mer$GRIN[i]) & mer$GRIN[i]=="unknown") { 
    mer$Comments_note[i] <- "in_question_NPGS_unknown"
  }
}
write.csv(mer[,c(1:5)],file="ames_isu_gb_npgs38_npgs83.csv", row.names = F) ### with four accessions still in question
### finished!

###### generate ames_1762_population_structure file

ames_isu <- read.csv("ames_isu_gb_npgs38_npgs83.csv",stringsAsFactors=FALSE)
ames_isu$ID[1] <- "PI 638551" ### change acc ID for DE3
ames_isu$ID_cleanformat[1] <- "pi638551" ### change acc ID for DE3

ames_2013gb <- read.csv('13059_2013_3103_MOESM1_ESM.csv', header = TRUE) # accession list from Cinta's paper
### adding clean format ID for matching
ames_2013gb$ID_cleanformat <- ames_2013gb$Accesion.N
ames_2013gb$ID_cleanformat <- gsub(' ','',ames_2013gb$ID_cleanformat)
ames_2013gb$ID_cleanformat <- tolower(ames_2013gb$ID_cleanformat)

ames_isu_gb <- merge(ames_isu,ames_2013gb,by="ID_cleanformat",all=T)
ames_isu_gb <- ames_isu_gb[which(!is.na(ames_isu_gb$ID)),]

### correct dup error from Cinta's list
ames_isu_gb_dup <- ames_isu_gb[duplicated2(ames_isu_gb$ID_cleanformat),]
ames_isu_gb <- ames_isu_gb[-which(ames_isu_gb$ID_cleanformat == "ames27101" & ames_isu_gb$Breeding.program == "Ontario"),]
ames_isu_gb <- ames_isu_gb[-which(ames_isu_gb$ID_cleanformat == "pi543850" & ames_isu_gb$Breeding.program == "Other"),]
ames_isu_gb <- ames_isu_gb[,-c(7:11)]
ames_isu_gb$Comments_isu_gb_s <- ames_isu_gb$Comments_isu_gb

ames_isu_gb$Comments_merged <- paste(ames_isu_gb$Comments_isu_gb,ames_isu_gb$Pop.structure,sep="_")
table(ames_isu_gb$Comments_merged)

ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_landraces")] <- "landraces"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_NA")] <- "unknown"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_non-stiff stalk")] <- "non-stiff stalk"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_stiff stalk")] <- "stiff stalk"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_tropical")] <- "tropical"
ames_isu_gb$Comments_isu_gb[which(ames_isu_gb$Comments_merged == "other_unclassified")] <- "unclassified"

write.csv(ames_isu_gb[,c(1:3,10,4:8)],file="ames_1762_pop_str_20191024.csv", row.names = F)

#################################     (8) compare sweet corn and popcorn list with Laura's    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

sp <- read.csv('ames_isu_gb_npgs38_npgs83.csv', header = TRUE)
sp_remove <- sp[which(sp$Comments_isu_gb!="other"|sp$Comments_note!="ok"),] ### 246 IDs will be removed from our side
lau <- read.csv('sweet.and.pop.for.Di.Source_Laura.csv', header = TRUE) ### 246 IDs will be removed from Laura's side

sp[duplicated2(sp$ID),] ### unique
lau[duplicated2(lau$Source),] ### unique

mer <- merge(sp,lau,by.x="ID",by.y="Source",all=T)
### check whether sp contain additional sweet or pop
mer_oth <- mer[is.na(mer$Gen),]
mer_oth_sp <- mer_oth[which(mer_oth$Comments_isu_gb!="other"),] ### sp did not contain additional sweet or pop
mer_oth_sp2 <- mer_oth[which(mer_oth$Comments_note!="ok"),] ### sp did not contain additional sweet or pop

### check consistency between sp and lau
mer_con <- mer[!is.na(mer$Gen),]
mer_con$comment_merge <- paste(mer_con$Comments_isu_gb,"#",mer_con$Comments,sep="")
table(mer_con$comment_merge)
ques <- mer_con[which(mer_con$Comments_isu_gb=="other"),] ### these 2 lines are in-question from the Comments_note
### finished: our sweet&pop list is the same as Laura's

#################################     (9) remove sweet corn and popcorn from full list    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

### read in d1762 full list with finalized comments for sweet corn and popcorn
sp <- read.csv('ames_isu_gb_npgs38_npgs83.csv', header = TRUE)
### read in d1726 list from dup_599 directory
d1726 <- read.csv('ames_1726_geno_pheno_ID.csv', header = TRUE)

mer <- merge(d1726,sp,by = "ID", all =T) ### For line DE3, PI 638551 was use as its source ID to retrive GBS entry,
                                         ### while 13A089451A was used as its source ID in the field record. PI 638551 was retained in subsequent analysis;
                                         ### for B73, B73(PI550473):250024957 was used to retrive GBS entry and retained
mer2 <- mer[which(!is.na(mer$Comments)),]
mer2 <- mer2[,c(1,2,4,5,8,9)]
colnames(mer2)[2] <-  "Pedigree"
### fill the "NA" cells introduced by merge
mer2$Comments_isu_gb[which(mer2$Pedigree == "B73"|mer2$Pedigree == "DE3")] <- "other"
mer2$Comments_note[which(mer2$Pedigree == "B73"|mer2$Pedigree == "DE3")] <- "ok"
table(mer2$Comments_note)
### remove lines with comment in question
mer3 <- mer2[which(mer2$Comments_note=="ok"),] ### 4 IDs in-question were removed
table(mer3$Comments_isu_gb)
### remove sweetcorn and popcorn
mer4 <- mer3[-which(mer3$Comments_isu_gb=="sweet corn"|mer3$Comments_isu_gb=="popcorn"),] ### 238 sweet or pop IDs were removed

write.csv(mer4,file="ames_1484_geno_pheno_ID.csv", row.names = F) ### 1484 IDs after removing sweet corn and popcorn
### finished!

#################################     (10) retrive common line list of GBS records from Panzea_ZeaGBSv2.7    ########################################
### also note that DE3's source is known as both "13A089451A" and "PI638551" in different files (from Laura)

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

ames_1762 <- read.csv(file="ames_isu_gb_npgs38_npgs83.csv",stringsAsFactors=FALSE)
ames_1762$ID[which(ames_1762$ID=="13A089451A")] <- as.character("PI 638551")
ames_1516 <- ames_1762[which(ames_1762$Comments_isu_gb=="other"&ames_1762$Comments_note=="ok"),]

ames_1516$ID_clean <- ames_1516$ID
ames_1516$ID_clean <- gsub(' ','',ames_1516$ID_clean)
ames_1516$ID_clean <- tolower(ames_1516$ID_clean)

v2.7 <- read.csv(file="Panzea_ZeaGBSv2.7_id.csv")
v2.7$acc_v2.7_clean <- v2.7$acc_v2.7
v2.7$acc_v2.7_clean <- gsub(' ','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('-','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('_','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\(','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\)','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- tolower(v2.7$acc_v2.7_clean)

ames_1516_v2.7 <- merge(ames_1516, v2.7, by.x="ID_clean", by.y="acc_v2.7_clean",all = T)
ames_1516_v2.7 <- ames_1516_v2.7[which(!is.na(ames_1516_v2.7$Pedigree)),]
ames_1516_v2.7_na <- ames_1516_v2.7[which(is.na(ames_1516_v2.7$INDV)),] # 4 IDs without GBS records

ames_1516_v2.7_all <- ames_1516_v2.7[which(!is.na(ames_1516_v2.7$INDV)),] 
ames_1516_v2.7_uniq <- distinct(ames_1516_v2.7_all, ID, .keep_all = TRUE) # 1801 GBS records for 1512 IDs

write.csv(ames_1516_v2.7_all[,-1],file="ames_1512_v2.7_all.csv", row.names = F)
ames_1512 <- read.csv("ames_1512_v2.7_all.csv",h=T)
ames_1512_dup <- ames_1512[duplicated2(ames_1512$ID_clean),]
write.csv(ames_1512_dup[,-1],file="ames_1512_v2.7_all_dup_520.csv", row.names = F)  # dup_520 list

###
dup_520 <- read.csv(file="ames_1512_v2.7_all_dup_520.csv",h=T)
dup_f <- dup_520[which(dup_520$F_MISS < 0.6),]
dup_f <- dup_f[duplicated2(dup_f$ID),]
dup_f_u <- distinct(dup_f, ID, .keep_all = TRUE)
write.csv(dup_f,file="ames_1512_v2.7_all_dup_514.csv", row.names = F)  # dup_514 list

#################################     (11) remove "ae" and seven accessions to generate locked list    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

ames_1762 <- read.csv(file="ames_1762_pop_str_20191024.csv",stringsAsFactors=FALSE)[,1:6]
ames_1762$ID <- gsub(" ","",ames_1762$ID)
ames_1516 <- ames_1762[which(ames_1762$Comments_isu_gb!="popcorn"&ames_1762$Comments_isu_gb!="sweet corn"&ames_1762$Comments_note=="ok"),]
### remove "ae" lines

### adding ae comment
ames_1516$ae <- grepl('ae', ames_1516$Pedigree)
ames_1516$ae <- gsub('TRUE','ae',ames_1516$ae)
ames_1516$ae <- gsub('FALSE','not_ae',ames_1516$ae)
ames_1504 <- ames_1516[-which(ames_1516$ae=="ae"),]

### remove seven lines
seven <- read.csv("seven_accessions_removed_by_kernel_type.csv")
ames_1497 <- ames_1504[-which(ames_1504$ID%in%seven$Accession),]
write.csv(ames_1497[,2:6],"ames_vita_1497_pheno_locked_list_20191118.csv",row.names = F)

#################################     (12) retrive common line list of GBS records from Panzea_ZeaGBSv2.7    ########################################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)

ames_1497 <- read.csv("ames_vita_1497_pheno_locked_list_20191118.csv")
ames_1497$ID_clean <- ames_1497$ID
ames_1497$ID_clean <- tolower(ames_1497$ID_clean)

v2.7 <- read.csv(file="Panzea_ZeaGBSv2.7_id.csv")
v2.7$acc_v2.7_clean <- v2.7$acc_v2.7
v2.7$acc_v2.7_clean <- gsub(' ','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('-','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('_','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\(','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- gsub('\\)','',v2.7$acc_v2.7_clean)
v2.7$acc_v2.7_clean <- tolower(v2.7$acc_v2.7_clean)

ames_1497_v2.7 <- merge(ames_1497, v2.7, by.x="ID_clean", by.y="acc_v2.7_clean",all = T)
ames_1497_v2.7 <- ames_1497_v2.7[which(!is.na(ames_1497_v2.7$Pedigree)),]
ames_1497_v2.7_na <- ames_1497_v2.7[which(is.na(ames_1497_v2.7$INDV)),] # 4 IDs without GBS records

ames_1497_v2.7_all <- ames_1497_v2.7[which(!is.na(ames_1497_v2.7$INDV)),] 
ames_1497_v2.7_uniq <- distinct(ames_1497_v2.7_all, ID, .keep_all = TRUE) # 1779 GBS records for 1493 IDs
write.csv(ames_1497_v2.7_all[,-1],file="ames_1493_v2.7_all.csv", row.names = F)

ames_1493 <- read.csv("ames_1493_v2.7_all.csv",h=T)
ames_1493_dup <- ames_1493[duplicated2(ames_1493$ID),]
write.csv(ames_1493_dup,file="ames_1493_v2.7_all_dup_514.csv", row.names = F)  # dup_514 list
### finished!
