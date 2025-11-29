install.packages("openxlsx")

info = read.table('C:/Users/ohuyh/OneDrive/바탕 화면/DATA/Sample_Information.csv', header = TRUE, row.names = 1, sep='\t')
ex = read.table('C:/Users/ohuyh/OneDrive/바탕 화면/DATA/Expression_Table.csv', header = TRUE, row.names = 1, sep='\t')
annot = read.table('C:/Users/ohuyh/OneDrive/바탕 화면/DATA/Annotations.csv', header = TRUE, row.names = 1, sep='\t')
sh = read.table('C:/Users/ohuyh/OneDrive/바탕 화면/DATA/DE_SA_vs_HC.csv', header = TRUE, row.names = 1, sep='\t')
gh = read.table('C:/Users/ohuyh/OneDrive/바탕 화면/DATA/DE_GOUT_vs_HC.csv', header = TRUE, row.names = 1, sep='\t')

## summary statistics for the clinical information 

#Neutrophil
Neu_GO = info[10:18,3] 
MeanNeu_GO = round(mean(Neu_GO),2)
SDNeu_GO = round(sd(Neu_GO),2)

Neu_SA = info[19:27,3] 
MeanNeu_SA = round(mean(Neu_SA),2)
SDNeu_SA = round(sd(Neu_SA),2)

pGS = t.test(Neu_GO,Neu_SA)
pGS = pGS$p.value
pGS

#SEX 
info_GO = info[10:18,]
sex_GO_female = subset(info_GO,SEX=="F")
num_GO_female = nrow(sex_GO_female)
sex_GO_male = subset(info_GO,SEX=="M")
num_GO_male = nrow(sex_GO_male)

info_SA = info[19:27,]
sex_SA_female = subset(info_SA,SEX=="F")
num_SA_female = nrow(sex_SA_female)
sex_SA_male = subset(info_SA,SEX=="M")
num_SA_male = nrow(sex_SA_male)

#---------------------------------------------------------------------------------------------------------

##Extract significant genes from HC vs GOUT and HC VS SA  
#gh
gh_sig = subset(gh,p.adj < 0.05)
gh_sig = subset(gh_sig,log2Fold >= 1)
summary(gh_sig)
#sh
sh_sig = subset(sh,p.adj < 0.05)
sh_sig = subset(sh_sig,log2Fold >= 1)
summary_sh = summary(sh_sig)



##Create Expression level Tables with significant genes 
# significant genes from two dataset 
geneID_gh = row.names(gh_sig)
geneID_sh = row.names(sh_sig)

#samples 
HC = names(ex[1:9])
GO = names(ex[10:18])
SA= names(ex[19:27])

allsamples = names(ex[1:27])
GH = names(ex[1:18])

# HC vs GOUT
GH_ID_sig = ex[geneID_gh,GH] 
GH_sig_all=ex[geneID_gh,allsamples]


# HC vs SA 
HA_ID_sig = ex[geneID_sh,HC]
SA_ID_sig = ex[geneID_sh,SA]
SH_ID_sig = merge(HA_ID_sig,SA_ID_sig, by.x=0,by.y=0)
row.names(SH_ID_sig)=SH_ID_sig[,1]
SH_ID_sig = SH_ID_sig[-1]

SH_sig_all=ex[geneID_sh,allsamples]

## In HC vs GOUT(GH_ID_sig), Make a plot of the expression levels of all significant genes in the each group (HC, GOUT)

install.packages("ggplot2")
library(ggplot2)

# calculate the mean  from HC and GOUT 
Mean_Healthy = round(rowMeans(GH_ID_sig[,HC]),2)
Mean_GOUT  = round(rowMeans(GH_ID_sig[,GO]),2)

Mean_GH_sig = data.frame(Mean_Healthy,Mean_GOUT)

# Make the histogram and density 

HC_G = ggplot(Mean_GH_sig,aes(x=log10(Mean_Healthy))) + 
  geom_histogram(aes(y=..density..),binwidth = 0.2, color = 'black',fill="white",size =0.5)+
  geom_density(fill='lightgreen',alpha=0.4) + 
  labs(x="log10(Mean value of Healthy)", y="Density", title ="HC vs GOUT")+
  scale_x_continuous(limits = c(-0.6 ,4.0))+
  scale_y_continuous(limits = c(0,1.1))

HC_G

GO_H = ggplot(Mean_GH_sig,aes(x=log10(Mean_GOUT))) + 
  geom_histogram(aes(y=..density..),binwidth = 0.2, color = 'black',fill="white",size =0.5)+
  geom_density(fill='skyblue',alpha=0.4) + 
  labs(x="log10(Mean value of GOUT )", y="Density", title ="HC vs GOUT")+
  scale_x_continuous(limits = c(-0.6,4.0))+
  scale_y_continuous(limits = c(0,1.1))
GO_H


## In HC vs SA(SH_ID_sig), Make a plot of the expression levels of all significant genes in the each group (HC, SA)
# calculate the mean  from HC and GOUT 
Mean_Healthy = round(rowMeans(SH_ID_sig[,HC]),2)
Mean_SEPSIS  = round(rowMeans(SH_ID_sig[,SA]),2)

Mean_SH_sig = data.frame(Mean_Healthy,Mean_SEPSIS)

# Make the histogram and density 

HC_S = ggplot(Mean_SH_sig,aes(x=log10(Mean_Healthy))) + 
  geom_histogram(aes(y=..density..),binwidth = 0.2, color = 'black',fill="white",size =0.5)+
  geom_density(fill='lightgreen',alpha=0.4) + 
  labs(x="log10(Mean value of Healthy)", y="Density", title ="HC vs SEPSIS")+
  scale_x_continuous(limits = c(-2.5 ,6))+
  scale_y_continuous(limits = c(0,0.4))
HC_S

SA_H = ggplot(Mean_SH_sig,aes(x=log10(Mean_SEPSIS))) + 
  geom_histogram(aes(y=..density..),binwidth = 0.2, color = 'black',fill="white",size =0.5)+
  geom_density(fill='salmon',alpha=0.3) + 
  labs(x="log10(Mean value of SEPSIS )", y="Density", title ="HC vs SEPSIS")+
  scale_x_continuous(limits = c(-2.5, 6))+
  scale_y_continuous(limits = c(0,0.4))
SA_H

#---------------------------------------------------------------------------------------------------

##Are these significant genes affected by any of the clinical measurements?

# ENSG00000179023 : have the lowest p.adj value in HC vs GOUT 

gene179 = ex["ENSG00000179023",]
gene179 = t(gene179)
gene179 = data.frame(gene179)

gene179$sex = info$SEX
gene179$neutropils = info$NEUTROPHILS

##Effect on SEX## 
# density plot
p179_sex_d = ggplot(gene179,aes(x=ENSG00000179023,fill=sex))+
            geom_density(alpha=0.4) + labs(x="ENSG00000179023", y="Density", title ="179_sex_density")
p179_sex_d

# box plot
p179_sex_b = ggplot(gene179,aes(x=sex,y=ENSG00000179023, fill=sex))+
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) + 
  stat_summary(fun=mean, colour="red") + 
  labs(y="Expression level", title ="ENSG00000179023")
p179_sex_b




## Neutrophil ## 
round(cor(gene179$neutropils,gene179$ENSG00000179023),2) 
model179 = lm(gene179$ENSG00000179023~gene179$neutropils)
anova(model179)


# ENSG00000198074 : have the lowest p.adj value in HC vs SEPSIS
gene198 = ex["ENSG00000198074",]
gene198 = t(gene198)
gene198= data.frame(gene198)

gene198$sex = info$SEX
gene198$neutropils = info$NEUTROPHILS

##Effect on SEX## 
# density plot
p198_sex_d = ggplot(gene198,aes(x=ENSG00000198074,fill=sex))+
  geom_density(alpha=0.4) + labs(x="ENSG00000198074", y="Density", title ="198_sex_density")
p198_sex_d

# box plot
p198_sex_b = ggplot(gene198,aes(x=sex,y=ENSG00000198074, fill=sex))+
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) + 
  stat_summary(fun=mean, colour="red") + 
  labs(y="Expression level", title ="ENSG00000198074")
p198_sex_b



## Neutrophil ## 
round(cor(gene198$neutropils,gene198$ENSG00000198074),2) 
model198 = lm(gene198$ENSG00000198074~gene198$neutropils)
anova(model198)
summary(model198)

g198_neu = ggplot(gene198, aes(x=neutropils, y=ENSG00000198074)) +  
  geom_point() +  
  labs(x="Neutrophils", y="Expression", title ='ENSG00000198074' ) + 
  geom_smooth(method = "lm", se = FALSE) 
g198_neu


## Make the plots for observing the distribution of significant genes 

#change the columns and rows 
GH_ID_t = t(GH_sig_all)
GH_ID_t = data.frame(GH_ID_t)
SH_ID_t = t(SH_sig_all)
SH_ID_t = data.frame(SH_ID_t)

##merging the sample group column from  info data to specify which group each sample belongs to 
samplegroup = info[,1]

GH_ID_t = cbind(samplegroup,GH_ID_t)
colnames(GH_ID_t)[1] = "Sample Group" 

GS_ID_sig = GH_ID_t[10:27,] # Include both GOUT and SA group based on significant genes of HC vs GOUT
GOUT_GH = GH_ID_t[GH_ID_t$`Sample Group`=="GOUT",-1] #only GOUT group on significant genes of HC vs GOUT
SA_GH = GH_ID_t[GH_ID_t$`Sample Group`=="SEPSIS",-1 ] #only SA group on significant genes of HC vs GOUT

SH_ID_t = cbind(samplegroup,SH_ID_t)
colnames(SH_ID_t)[1] = "Sample Group" 

SG_ID_sig = SH_ID_t[10:27,]# Include both GOUT and SA group based on significant genes of HC vs SA
GOUT_SH = SH_ID_t[SH_ID_t$`Sample Group`=="GOUT",-1 ] #only GOUT group on significant genes of HC vs SA
SA_SH = SH_ID_t[SH_ID_t$`Sample Group`=="SEPSIS",-1 ] #only SA group on significant genes of HC vs SA

#----------------------------------------------------------------------------------------------------------
## Exclude the significant genes between GOUT vs SA among significant genes HC VS Gout

#for confirmation of the loop. 
GOUT_179 = GOUT_GH[,"ENSG00000179023"]
SA_179 = SA_GH[,"ENSG00000179023"]

meanG = mean(GOUT_179)
meanS = mean(SA_179)
log2= log2(meanG) - log2(meanS)
log2


# Make the dataset(de_GS), which has log2fold, p value from GOUT and SA samples about HC VS GOUT significant genes 
# Make the loop 
de_GS = as.data.frame(matrix(0,ncol = 2, nrow=nrow(GH_ID_sig)))
de_GS
names(de_GS) = c("log2fold",'p')
row.names(de_GS) = row.names(GH_ID_sig)

for (num in 1:nrow(GH_ID_sig))
{
  mean_GOUT = mean(GOUT_GH[,num])
  mean_SA = mean(SA_GH[,num])
  log2fold = log2(mean_GOUT) - log2(mean_SA)
  p = t.test(GOUT_GH[,num],SA_GH[,num])
  p = p$p.value
  de_GS[num,'log2fold'] = log2fold
  de_GS[num,'p'] = p
}

#sort the significant genes, p>0.05 and log2fold >=2

de_GS = subset(de_GS,p<0.05)
de_GS = subset(de_GS,log2fold >=2)

#significant genes 
de_GS_sig = row.names(de_GS)
de_GS_sig 


## Exclude the significant genes between GOUT vs SA among significant genes HC VS SA


# Make the dataset(de_SG), which has log2fold, p value from GOUT and SA samples about HC VS SA significant genes 
# Make the loop 
de_SG = as.data.frame(matrix(0,ncol = 2, nrow=nrow(SH_ID_sig)))
de_SG
names(de_SG) = c("log2fold",'p')
row.names(de_SG) = row.names(SH_ID_sig)

for (num in 1:nrow(SH_ID_sig))
{
  mean_GOUT = mean(GOUT_SH[,num])
  mean_SA = mean(SA_SH[,num])
  log2fold = log2(mean_SA) - log2(mean_GOUT)
  p = t.test(GOUT_SH[,num],SA_SH[,num])
  p = p$p.value
  de_SG[num,'log2fold'] = log2fold
  de_SG[num,'p'] = p
}

#change INF to 0 of on log2fold columns

de_SG$log2fold[is.infinite(de_SG$log2fold)] = 0 

#sort the significant genes, p>0.05 and log2fold >=2
de_SG = subset(de_SG,p<0.05)
de_SG = subset(de_SG,log2fold >=2)

#significant genes 
de_SG_sig = row.names(de_SG)
de_SG_sig 

## find the common gene -> result = Empty. 
common_gene = intersect(de_GS_sig,de_SG_sig)
## -> There is no common gene between GH and SH 


## Distribution of significant genes 


# significant GOUT genes ("ENSG00000071909" "ENSG00000267934")
G7190= ggplot(GS_ID_sig,aes(x=`Sample Group`, y=ENSG00000071909, colour = `Sample Group`))+
      geom_boxplot() + 
      geom_point() +
      coord_cartesian(ylim=c(0,230))
      
G7190

G2679= ggplot(GS_ID_sig,aes(x=`Sample Group`, y=ENSG00000267934, colour = `Sample Group`))+
  geom_boxplot() + 
  geom_point() +
  coord_cartesian(ylim=c(0,230))
G2679

# significant SA genes ("ENSG00000065328" "ENSG00000012048")->Two genes picked by the order of the highest log2fold values
G6532= ggplot(SG_ID_sig,aes(x=`Sample Group`, y=ENSG00000065328, colour = `Sample Group`))+
  geom_boxplot() + 
  geom_point()
G6532

G1204= ggplot(SG_ID_sig,aes(x=`Sample Group`, y=ENSG00000012048, colour = `Sample Group`))+
  geom_boxplot() + 
  geom_point()
G1204

#--------------------------------------------------------------------------------------------------------
# Download the dataframes 

library(openxlsx)
write.xlsx(sh_sig,"sh_sig.xlsx", row.names =TRUE)
write.xlsx(gh_sig,"gh_sig.xlsx", row.names =TRUE)
write.xlsx(de_GS,"de_GS.xlsx", row.names =TRUE)
write.xlsx(de_SG,"de_SG.xlsx", row.names =TRUE)
