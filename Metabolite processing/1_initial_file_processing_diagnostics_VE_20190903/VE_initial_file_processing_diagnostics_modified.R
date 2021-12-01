### this script file (modified from Di Wu's by Xiaowei Li) includes following analyses:
# (1) format raw data
# (2) calculate correlation between years
# (3) add and format other columns for wrapper
# (4) check plate variation

########################################################################################
##############################     (1) format raw data     #############################

###load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)

### read in files
d2015 <- read.csv('Ames_2015_ve_adding_ID.csv',header = TRUE,stringsAsFactors=FALSE)[,-26] # remove PC8 trait, 1801 rows
d2015$uni <- paste(d2015$Year,d2015$Pass,d2015$Range,sep="_")
d2017 <- read.csv('Ames_2017_ve_adding_ID.csv', header = TRUE) # 1738 lines
d2017$uni <- paste(d2017$Year,d2017$Pass,d2017$Range,sep="_")

### read in updated key file (adding two lines)
# By that logic, it looks like the N215 is in New_Block_Range 5 and New_Block_Pass 1, and the B73 is in New_Block_Range 7 and New_Block_Pass 3. (comment from Laura)
newkey <- read.csv('Source_to_Block_key_15-17_adding_2_lines.csv',header = TRUE,stringsAsFactors=FALSE)
newkey$uni <- paste(newkey$Year,newkey$Pass,newkey$Range,sep="_")
newkey <- newkey[,c(2,6,9,10,12,13)]
### remove dups in key file
newkey_dup <- newkey[duplicated2(newkey$uni),] # 10 samples were measured twice
### to keep consistecy, remove those 10 measurements from Plate 19
remove_id <- newkey_dup[which(newkey_dup$plate==19),]
all_row <- NULL
for (i in 1:10){
  tmp_row <- which(newkey$plate==19 & newkey$uni==remove_id$uni[i] )
  all_row <- c(all_row,tmp_row)
}
newkey <- newkey[-all_row,]

### adding new_pass and new_range for d2015
d2015_n <- merge(d2015,newkey,by="uni")
emp <- d2015_n[which(d2015_n$ID!=d2015_n$Source.y),]
d2015_n <- d2015_n[,c(2:14,32:33,15:29)]
colnames(d2015_n) <- gsub("\\.x","",colnames(d2015_n))
write.csv(d2015_n,'Ames_2015_ve_adding_ID_new_block.csv',row.names = F) ### 1801 rows

### adding new_pass and new_range for d2017
d2017_n <- merge(d2017,newkey,by="uni")
emp <- d2017_n[which(d2017_n$ID!=d2017_n$Source.y),]
d2017_n <- d2017_n[,c(2:14,32:33,15:29)]
colnames(d2017_n) <- gsub("\\.x","",colnames(d2017_n))
write.csv(d2017_n,'Ames_2017_ve_adding_ID_new_block.csv', row.names = F)

### read in new files
d2015 <- read.csv('Ames_2015_ve_adding_ID_new_block.csv', header = TRUE) # 1801 rows
d2017 <- read.csv('Ames_2017_ve_adding_ID_new_block.csv', header = TRUE) # 1738 rows

### check colnames before renaming
n2015 <- colnames(d2015)
n2017 <- colnames(d2017)
nnr <- as.data.frame(rbind(n2015,n2017))
(nnrdif <- nnr[ ,nnr[1,]!=nnr[2,]])

### rename
names(d2017) <- names(d2015)

### 2015 check-line
check <- which(d2015$Event.Name=='Check')
B73 <- which(d2015$Pedigree=='B73')
which(check!=B73) # should be 0
### 2017 check-line
check <- which(d2017$Event.Name=='Check')
B73 <- which(d2017$Pedigree=='B73')
which(check!=B73) # should be 0

###### plot field map
d2015$check <- 99
d2015$check[which(d2015$Event.Name=='Check')] <- 1
d2017$check <- 99
d2017$check[which(d2017$Event.Name=='Check')] <- 1
d2015$block <- paste(d2015$New_Block_Pass,d2015$New_Block_Range,sep="_")
d2017$block <- paste(d2017$New_Block_Pass,d2017$New_Block_Range,sep="_")

pdf('field_map_new_range_row.pdf', width=12, height=10)
ggplot(d2015, aes(Range, Pass)) + 
  geom_tile(aes(fill = factor(block)),colour = "white")+
  labs(title="2015")+
  scale_fill_manual(values = rep(c("cornflowerblue","yellow","forestgreen","orange","cyan"),500))+
  geom_text(aes(label = block),size=2)+
  scale_colour_manual(values = c("red","black"))+
  geom_point(data = d2015[which(d2015$check==1),],mapping = aes(Range, Pass),color="red",size=3)+
  theme(legend.position = "none")

ggplot(d2017, aes(Range, Pass)) + 
  geom_tile(aes(fill = factor(block)),colour = "white")+
  labs(title="2017")+
  scale_fill_manual(values = rep(c("cornflowerblue","yellow","forestgreen","orange","cyan"),500))+
  geom_text(aes(label = block),size=2)+
  scale_colour_manual(values = c("red","black"))+
  geom_point(data = d2017[which(d2017$check==1),],mapping = aes(Range, Pass),color="red",size=3)+
  theme(legend.position = "none")

dev.off()

############################################################################################################
##############################     (2) calculate correlation between years     #############################

### merge file to get correlation
traits <- c( "d.T3", "g.T3", "a.T3", "d.T", "g.T", "a.T", "Total.Tocotrienols", "Total.Tocopherols", "Total.Tocochromanols" )

### remove B73 first, then merge using ID
curr_raw <- merge(d2015[-which(d2015$check==1),c('ID',"HPLC.Plate.Number",traits)],
               d2017[-which(d2017$check==1),c('ID',"HPLC.Plate.Number",traits)],by='ID') # intersection set contains 1630 lines
curr_raw_dup <- curr_raw[duplicated2(curr_raw$ID),] # no dup

### correlation and Scatterplot
 correlation <- matrix(NA,nrow=0,ncol=2)
 
 pdf('correlation_and_comparison_between_years_raw_ve.pdf',width = 10, height = 8)
 for (i in traits){
   index <- grep(paste(i,".",sep=''),names(curr_raw),fixed = TRUE) 
   curr <- as.data.frame(cbind(curr_raw[index[1]],curr_raw[index[2]]))
   ### correlation
   r <- cor(curr[,1],curr[,2],use='complete.obs')
   p <- cor.test(curr[,1],curr[,2])$p.value
   correlation <- rbind(correlation,c(r[1],p))
   ### overall correlation Scatterplot
   plot(curr[,2]~curr[,1],main=paste(i),xlab='2015',ylab='2017',pch=16)
   abline(0,1,col='red')
   text(x=max(curr[,1])*0.8,y=max(curr[,2])*0.8,paste('r =',round(r,2),sep=''))
  
   ### check 2015_plate effect
   g1 <- ggplot(curr, aes(curr[,1], curr[,2]))+
     geom_point(aes(col=curr_raw$HPLC.Plate.Number.x)) +
     geom_smooth(aes(col=curr_raw$HPLC.Plate.Number.x), method="lm", se=F)+
     scale_color_manual(values =rainbow(length(unique(curr_raw$HPLC.Plate.Number.x))))+
     theme_bw()+
     labs( title=paste(i," grouped by 2015_plate",sep=''))+
     ylab("2017") +
     xlab("2015") +
     theme(legend.title = element_text(size = 8,face = "bold",colour = "black"),
           legend.text = element_text(size = 5),
           plot.title = element_text(size=16,face = "bold",color="forestgreen"),
           axis.title=element_text(size=14,face="bold"),
           axis.text.x = element_text(size=10,face="bold",colour = "black",angle = 0),
           axis.text.y = element_text(size=10,face="bold",colour = "black",angle = 0))+
     guides(col = guide_legend(title = "2015_plate"))+
     annotate("text", x=max(curr[,1])*0.8,y=max(curr[,2])*0.8,label=paste('r = ',round(r,2),sep=''),size =5, colour = "black")
   print(g1)

   ### check 2017_plate effect
   g2 <- ggplot(curr, aes(curr[,1], curr[,2]))+
     geom_point(aes(col=curr_raw$HPLC.Plate.Number.y)) +
     geom_smooth(aes(col=curr_raw$HPLC.Plate.Number.y), method="lm", se=F)+
     scale_color_manual(values =rainbow(length(unique(curr_raw$HPLC.Plate.Number.y))))+
     theme_bw()+
     labs( title=paste(i," grouped by 2017_plate",sep=''))+
     ylab("2017") +
     xlab("2015") +
     theme(legend.title = element_text(size = 8,face = "bold",colour = "black"),
           legend.text = element_text(size = 5),
           plot.title = element_text(size=16,face = "bold",color="blue"),
           axis.title=element_text(size=14,face="bold"),
           axis.text.x = element_text(size=10,face="bold",colour = "black",angle = 0),
           axis.text.y = element_text(size=10,face="bold",colour = "black",angle = 0))+
     guides(col = guide_legend(title = "2017_plate"))+
     annotate("text", x=max(curr[,1])*0.8,y=max(curr[,2])*0.8,label=paste('r = ',round(r,2),sep=''),size =5, colour = "black")
   print(g2)

   ### histogram
   hist(curr[,1],col=rgb(255, 0, 0, max = 255, alpha = 100),breaks=50,xlab='',main=paste(i))
   hist(curr[,2], add=T,col=rgb(0, 0, 255, max = 255, alpha = 100),breaks=50)
   legend('topright',fill=c(rgb(255, 0, 0, max = 255, alpha = 100),rgb(0, 0, 255, max = 255, alpha = 100)),legend=c('2015','2017'))
 }
 dev.off()

### write correlation table
correlation <- as.data.frame(correlation)
row.names(correlation) <- traits
names(correlation) <- c('r','p')
write.csv(correlation,'correlation_ve_between_years_raw.csv')

##############################################################################################################
##############################     (3) add and format other columns for wrapper     ##########################

### adding is.experimental column
d2015$is.experimental <- 1
d2015$is.experimental[which(d2015$Event.Name=='Check')] <- 0
d2017$is.experimental <- 1
d2017$is.experimental[which(d2017$Event.Name=='Check')] <- 0

### adding plate column
table(d2015$HPLC.Plate.Number)
d2015$plate <- gsub('Ames15-CT','',d2015$HPLC.Plate.Number)
d2015$plate <- gsub('H','',d2015$plate)
d2015$plate <- as.numeric(d2015$plate)
table(d2015$plate) # should be the same

table(d2017$HPLC.Plate.Number)
d2017$plate <- gsub('Ames17-Toc-','',d2017$HPLC.Plate.Number)
d2017$plate <- gsub('Ames17_Toc-','',d2017$plate)
d2017$plate <- gsub('H','',d2017$plate)
d2017$plate <- gsub('L','',d2017$plate)
d2017$plate <- as.numeric(d2017$plate)
table(d2017$plate)

### merge two years together
### check colnames before merge
n2015 <- colnames(d2015)
n2017 <- colnames(d2017)
nnr <- as.data.frame(rbind(n2015,n2017))
(nnrdif <- nnr[ ,nnr[1,]!=nnr[2,]])
all <- rbind.data.frame(d2015,d2017)
all_raw <- all[,c('ID','Tier','Range','Pass','Year','check','is.experimental','plate',traits)]
all_new <- all[,c('ID','Tier','New_Block_Range','New_Block_Pass','Year','check','is.experimental','plate',traits)]

write.table(all_raw,'ve_both_years_raw.txt',sep='\t',row.names = F,quote=F)
write.table(all_new,'ve_both_years_new_block.txt',sep='\t',row.names = F,quote=F)

###########################################################################################
#############################     (4) check plate variation     ###########################

### load libraries
library(gdata)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(dplyr)

raw <- read.delim('ve_both_years_raw.txt')
raw_ck <- raw[which(raw$is.experimental=='0'),]

### whether every plate with check?
table(raw_ck[which(raw_ck$Year=='2015'),]$plate)
table(raw_ck[which(raw_ck$Year=='2017'),]$plate) # 11th plate no checks

### whether check source unique?
raw_ck$ckid <- paste(raw_ck$Year,raw_ck$Range,raw_ck$Pass, sep="_")
raw_ck_dup <- raw_ck[duplicated2(raw_ck$ckid),]

### box plot to show trait range for each plate
pdf('2015_2017_Ames_ve_plate_box_raw.pdf', width=12, height=8)

for (i in 9:17){
b <- ggplot(raw, aes(x=factor(plate), y=raw[,i],fill=plate))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.9),
             aes(group=plate,colour = factor(is.experimental)))+
  scale_color_manual(values =c('red','black'),
                     name="Lines type",
                     breaks=c("0", "1"),
                     labels=c("Control_B73", "Ames_line"))+
  facet_grid(.~Year)+
  theme_bw()+
  labs('')+
  ylab(colnames(raw)[i])+
  xlab("LCMS plate") +
  theme(legend.title = element_text(size = 8,face = "bold",colour = "black"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size=10,face = "bold",color="forestgreen"),
        axis.title=element_text(size=10,face="bold"),
        axis.text.x = element_text(size=9,face="bold",colour = "black",angle = 0),
        axis.text.y = element_text(size=9,face="bold",colour = "black",angle = 0))
print(b)

}

dev.off()


