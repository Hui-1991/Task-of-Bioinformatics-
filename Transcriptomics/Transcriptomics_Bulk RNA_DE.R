library(ashr)
library(DESeq2) 
library(devEMF)
library(extrafont)
library(ggplot2)
library(ggrepel)
library(limma)
library(vsn)
library(hexbin)


###################################################################
####                        Functions                          ####
###################################################################

# For creating PCA plot using ggplot 
pca_plot = function(rlog, intgroup, color, ntop=500)
{
  pca_data = plotPCA(rlog, intgroup=intgroup, returnData=TRUE)
  
  # Calculate percentage variance explained
  pca_variance = round(100 * attr(pca_data, "percentVar"), 1)
  
  # Insert each sample name 
  pca_data$sample = rownames(pca_data)
  
  # Create PCA plot with variance in axis labels
  ggplot(pca_data, aes(PC1, PC2, color=color)) +
    geom_point(size=3) +  # Set dot size
    geom_text_repel(aes(label = sample), size=5, show.legend = FALSE)+
    xlab(paste0("PC1 ", "(",pca_variance[1],"%)")) +  # Add PC1 variance
    ylab(paste0("PC2 ","(",pca_variance[2],"%)"))+    # Add PC2 variance
    xlim(-5,7.5)+ # adjust x axis 
    ylim(-3,4)+   # adjust y axis
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(size = 14),  # Adjust axis label size
      legend.title = element_blank(),  # Removes legend title
      legend.position = "bottom",      # change lengend position
      legend.text = element_text(size =12), 
      legend.key = element_rect(fill = 'white', color = 'white'),
      panel.background = element_rect(fill = 'white', colour = 'black') # adjust paenl background
      )
}  

##################################################################
####         Read gene and transcript matrices               ##### 
##################################################################

geneTable= read.csv('gene_count_matrix.csv', row.names = 1)
transcriptTable = read.csv('transcript_count_matrix.csv', row.names = 1)


##################################################################
####               Read experiment design table               ####
##################################################################
# For one factor : Group 
colTable.g = read.csv("Design.Group.csv", row.names = 1 )
# For batch effect (two factor): Batch and Group 
colTable.gb = read.csv("Design.Group.Batch.csv", row.names = 1 )



##################################################################
####        Create dds object for each matrix                #####
##################################################################

#---------------------- gene matrix -----------------------------#   
# dds with gene matrics and batch effect 
dds.gene = DESeqDataSetFromMatrix(countData = geneTable,
                                  colData = colTable.gb,
                                  design = ~ batch + group) 

# check gene's number
dim(dds.gene)[1]

#--------------------- transcript matrix -----------------------# 
# dds with transcript matrics and batch effect 
dds.transcript = DESeqDataSetFromMatrix(countData = transcriptTable,
                                        colData = colTable.gb,
                                        design = ~ batch + group)
# check gene's number
dim(dds.transcript)[1]



##################################################################
#####                 Pre-filtering                          #####
##################################################################

smallestSamplesize = 3 # assign least sample size

#---------------------- gene matrix -----------------------------# 
# Filter the genes haing a min 10 read in at least 3 samples 
filter = rowSums(counts(dds.gene) >=10) >= smallestSamplesize
dds.gene = dds.gene[filter,]

# check gene's number after pre-filtering
dim(dds.gene)[1]

#--------------------- transcript matrix ------------------------# 
# Filter the genes haing a min 10 read in at least 3 samples 
filter = rowSums(counts(dds.transcript) >=10) >= smallestSamplesize
dds.transcript = dds.transcript[filter,]

# check gene's number after pre-filtering
dim(dds.transcript)[1]



##################################################################
#####                      DE analysis                       #####
##################################################################

#---------------------- gene matrix -----------------------------#

dds.gene = DESeq(dds.gene)
dds.gene

#--------------------- transcript matrix ------------------------# 

dds.transcript = DESeq(dds.transcript)
dds.transcript


##################################################################
####                Create Dispersion plots                   ####
##################################################################

# Generate plot and save as png with dds.gene
png("dispersion_gene.png", width = 3000, height = 1500, res=300) 
par(family = "Arial") # assign font 
plotDispEsts(dds.gene,
             main ="Gene count matrix")
            
dev.off()  

# Generate plot and save as png with dds.transcript
png("dispersion_transcirpt.png", width = 3000, height = 1500, res=300) 
par(family = "Arial") # assign font 
plotDispEsts(dds.transcript,
             main ="Transcript count matrix")
        
dev.off()  



#################################################################
####                Create rlog-based PCA                   #####
#################################################################

#---------------------- gene matrix -----------------------------#
# perform rlog
rlog.gene = rlog(dds.gene, blind = TRUE)

#create PCA plot 
pca_gene.g = pca_plot(rlog.gene,intgroup = "group", color = rlog.gene$group)
  pca_gene.g
pca_gene.b = pca_plot(rlog.gene,intgroup = "batch", color = rlog.gene$batch)
pca_gene.b

# Save plots as png 
ggsave("pca_gene_g.png", plot = pca_gene.g, width = 10, height = 6, dpi = 300)
ggsave("pca_gene_b.png", plot = pca_gene.b, width = 10, height = 6, dpi = 300)


#--------------------- transcript matrix ------------------------#

# perform rlog
rlog.transcript = rlog(dds.transcript, blind = TRUE)

#create PCA plot 
pca_transcript.g = pca_plot(rlog.transcript,intgroup = "group", color = rlog.gene$group)
pca_transcript.g
pca_transcript.b = pca_plot(rlog.transcript,intgroup = "batch", color = rlog.gene$batch)
pca_transcript.b

# Save plots as png 
ggsave("pca_transcirpt_g.png", plot = pca_transcript.g, width = 10, height = 8, dpi = 300)
ggsave("pca_transcirpt_b.png", plot = pca_transcript.b, width = 10, height = 8, dpi = 300)



##################################################################
####                     Create meanSDplot                   #####
##################################################################

# Normalization -> Log2(n+1)
log2norm = normTransform(dds.gene)

# Generate meansdplot with log2 normalized data 
log2.meanplot = meanSdPlot(assay(log2norm))
log2.meanplot = log2.meanplot$gg + 
  ggtitle("Log2 Normalization") +          # Add title
  xlab("Rank(mean)")+
  ylab("Standatd Deviation")+
  ylim(c(0,4.5))+
  theme(
    text = element_text(family = "Arial"), # Change font
    legend.position = "bottom")            # Adjust legend location

log2.meanplot

# Generate meandsdplot with rlog data 
rlog.meanplot  = meanSdPlot(assay(rlog.gene))
rlog.meanplot = rlog.meanplot$gg + 
  ggtitle("Regularized Log Transformation")+  # Add title
  xlab("Rank(mean)")+
  ylab("Standatd Deviation")+
  ylim(c(0,4.5))+
  theme(
    text = element_text(family = "Arial"), # Change font
    legend.position = "bottom")            # Adjust legend location

rlog.meanplot

# Save plots as png 
ggsave("meanSdPlot_log2norm.png", plot = log2.meanplot, width = 8, height = 10, dpi = 300)
ggsave("meanSdPlot_rlog.png", plot = rlog.meanplot, width = 8, height = 10, dpi = 300)



###################################################################
#####         Extract DE result: Unshrunken LFCs               ####
###################################################################

# NULL hypothesis of LogfoldChange(LFC)=0 
res0.BA = results(dds.gene, contrast = c("group","B","A")) 
res0.CA = results(dds.gene, contrast = c("group","C","A"))

summary(res0.BA, alpha=0.05)
summary(res0.CA, alpha=0.05)
# NULL hypothesis of LogfoldChange(LFC)<1
res1.BA = results(dds.gene, contrast = c("group","B","A"), lfcThreshold = 1)
res1.CA = results(dds.gene, contrast = c("group","C","A"), lfcThreshold = 1)

summary(res1.BA, alpha=0.05)
summary(res1.CA, alpha=0.05)

# Generate MAplots and save as png 
## Compare results for B vs A between NULL hypotheses LFC=0 and LFC < 1 
png("MAplot_unshrunken_BA.png", width = 2000, height = 1000, res = 300)
par(mfrow = c(1,2))
DESeq2::plotMA(res0.BA, alpha=0.05,
               main = "B vs A (LFC = 0)", # title
               ylim = c(-10,10),          # y axis
               colSig = "darkgreen")      # dot color 
abline(h=c(-1,1), col="black")            
DESeq2::plotMA(res1.BA,alpha=0.05, 
               main = "B vs A (LFC < 1)",
               ylim = c(-10,10),
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")

dev.off()

## Compare results for C vs A between LFC=0 and LFC < 1 
png("MAplot_unshrunken_CA.png", width = 2000, height = 1000, res = 300)
par(mfrow = c(1,2))
DESeq2::plotMA(res0.CA, alpha=0.05, 
               main = "C vs A (LFC = 0)", 
               ylim = c(-10,10),
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")
DESeq2::plotMA(res1.CA,alpha=0.05,
               main = "C vs A (LFC < 1)",
               ylim = c(-10,10),
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")

dev.off()


###################################################################
####          Extract DE result: shrunken LFCs                 ####
###################################################################
results(dds.gene)
resultsNames(dds.gene)

# NULL hypothesis of LogfoldChange(LFC)=0 
## compare between group: B vs A and  C vs A 
res0.sBA = lfcShrink(dds.gene, coef = "group_B_vs_A", type = "ashr")
res0.sCA = lfcShrink(dds.gene, coef = "group_C_vs_A", type = "ashr")

summary(res0.sBA, alpha=0.05)
summary(res0.sCA, alpha=0.05)

# NULL hypothesis of LogfoldChange(LFC)=1
## compare between group: B vs A, and  C vs A 
res1.sBA = lfcShrink(dds.gene, coef = "group_B_vs_A", type = "ashr", lfcThreshold = 1)
res1.sCA = lfcShrink(dds.gene, coef = "group_C_vs_A", type = "ashr", lfcThreshold = 1)

summary(res1.sBA, alpha=0.05)
summary(res1.sCA, alpha=0.05)

## extract the result with adj.p < 0.05 & sort LFC 
## write csv file 
### B vs A 
res1.sBA.sig = subset(res1.sBA, res1.sBA$padj<0.05) 
res1.sBA.sig.sort = res1.sBA.sig[order(res1.sBA.sig$log2FoldChange),]
head(res1.sBA.sig.sort)     # check sorted table 
dim(res1.sBA.sig.sort)[1]   # check the number of sorted table 
write.csv(res1.sBA.sig.sort,file="Significant.Genes.BvsA.csv") # write csv file 

### C vs A 
res1.sCA.sig = subset(res1.sCA, res1.sCA$padj<0.05) 
res1.sCA.sig.sort = res1.sCA.sig[order(res1.sCA.sig$log2FoldChange),]
head(res1.sCA.sig.sort)     # check sorted table 
dim(res1.sCA.sig.sort)[1]   # check the number of sorted table 
write.csv(res1.sCA.sig.sort,file="Significant.Genes.CvsA.csv") # write csv file 

#-----------------------------------------------------------------#

# Create plotMA and save as png
## Compare results for B vs A between LFC=0 and LFC < 1 
png("MAplot_shrunken_BA.png", width = 2000, height = 1000, res = 300)
par(mfrow = c(1,2))
DESeq2::plotMA(res0.sBA, alpha=0.05, 
               ylim = c(-10,10),
               main = "B vs A (LFC = 0)",
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")
DESeq2::plotMA(res1.sBA,alpha=0.05, 
               ylim = c(-10,10),
               main = "B vs A (LFC < 1)",
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")
dev.off()

## Compare results for C vs A between LFC=0 and LFC < 1 
png("MAplot_shrunken_CA.png", width = 2000, height = 1000, res = 300)
par(mfrow = c(1,2))
DESeq2::plotMA(res0.sCA, alpha=0.05, 
               ylim = c(-10,10),
               main = "C vs A (LFC = 0)",
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")
DESeq2::plotMA(res1.sCA,alpha=0.05, 
               ylim = c(-10,10),
               main = "C vs A (LFC < 1)",
               colSig = "darkgreen")
abline(h=c(-1,1), col="black")
dev.off()


#################################################################
####                Remove batch effect                      ####
#################################################################

# Create dds.g object (without batch) and rlog
dds.g = DESeqDataSetFromMatrix(countData = geneTable,
                               colData = colTable.g,
                               design = ~ group)
filter = rowSums(counts(dds.g) >=10) >= smallestSamplesize
dds.g = dds.g[filter,] # pre-filtering 
dds.g = DESeq(dds.g) # Run DE analysis
rlog.g = rlog(dds.g, blind =  TRUE) # Perform rlog 


# Perform batch correction 
mydesign = model.matrix(design(dds.g), colData(dds.g))
batch.corrected = limma::removeBatchEffect(assay(rlog.g),
                                       batch = colData(dds.gene)$batch,
                                       design = mydesign)
# Export correct rlog-transformed matrix 
write.csv(batch.corrected, file = "BatchCorrected.Rlog.csv", quote = FALSE)

# Generate PCA plot with batch correctd 
## Convert to a DESeqTransfrom object befor creating PCA plot 
## -> reference  code : https://www.biostars.org/p/403053/#473727
batch.rlog = rlog.g 
assay(batch.rlog) = batch.corrected
pca_b.correctd = pca_plot(batch.rlog, 
                          intgroup = "group", 
                          color = colData(dds.gene)$group)
pca_b.correctd
# save as png 
ggsave("pca_b_corrected.png", plot = pca_b.correctd, width = 10, height = 6, dpi = 300)
