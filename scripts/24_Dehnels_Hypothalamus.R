###Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon
#Reads were preprocessed on Noctillio server, where they were trimmed using fastp and aligned/quantified using Kalisto
#First I want to get the quantification files onto R and convert transcript abundance to gene abundance
#Then I will create a handful of quant files
#First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(readr)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)

#Set up your working directories
indir <- "./"
outputdir <- "../analysis/DESeq2/"
##Create a transcript to gene key using the gff file for the shrew
sortxdb <- makeTxDbFromGFF("../data/0_refs/GCF_000181275.1_SorAra2.0_genomic.gff.gz")
k <- keys(sortxdb, keytype="TXNAME")       
sortx2gene <- select(sortxdb, k, "GENEID", "TXNAME")
sortx2gene <- sortx2gene[!duplicated(sortx2gene[,1]),]
sortx2gene <- na.omit(sortx2gene)
#gets rid of the XRs, which are some misc_rnas, and do not apppear to be associated with genes

cs_hyp_samples <- read.table("../1_ids/completecycle_hypothalamus.txt", header = T)
cs_hyp_files <-file.path("../data/4_kallisto", cs_hyp_samples$Sample_name, "abundance.tsv")
names(cs_hyp_files) <- paste0("sample_", cs_hyp_samples$Sample_name)
all(file.exists(cs_hyp_files))
#
cs_hyp.count.tsv <- tximport(cs_hyp_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_hyp.tpm.tsv <- tximport(cs_hyp_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
write.table(cs_hyp.tpm.tsv$abundance, "cs_hyp.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

###HYPOTHALAMUS
hyp_stages <- factor(c(cs_hyp_samples$Run))
hyp_organs <- factor(c(cs_hyp_samples$Condition))
hyp_full1 <- factor(c(cs_hyp_samples$Run))
hyp_sex <- factor(c(cs_hyp_samples$Sex))
hyp_stages_organ_frame <-cbind(as.data.frame(hyp_stages),as.data.frame(hyp_organs),as.data.frame(hyp_full1),as.data.frame(hyp_sex))
ggplot(hyp_stages_organ_frame, aes(x=hyp_full1))+
  geom_bar(stat = 'count')+
  theme_bw()
#3 RNA Quality Control
vis_hyp_samples <- read.table("~/CompleteShrew/data/samples/completecycle_hypothalamus_4vis.txt", header = T)
ggplot(vis_hyp_samples,aes(x=RIN,y=Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
#3B removing outlier with high reads
hyp_stages <- factor(c(cs_hyp_samples$Run))
hyp_stages <- hyp_stages[-6]
hyp_organs <- factor(c(cs_hyp_samples$Condition))
hyp_organs <- hyp_organs[-6]
hyp_full1 <- factor(c(cs_hyp_samples$Run))
hyp_full1 <- hyp_full1[-6]
hyp_sex <- factor(c(cs_hyp_samples$Sex))
hyp_sex <- hyp_sex[-6]
hyp_stages_organ_frame <-cbind(as.data.frame(hyp_stages),as.data.frame(hyp_organs),as.data.frame(hyp_full1),as.data.frame(hyp_sex))
#everything looks pretty good
#Make DESeq onject
colnames(cs_hyp.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
dds_hyp_all <- DESeqDataSetFromMatrix(round(subset(cs_hyp.count.tsv$counts, select=-Stg2_4)), DataFrame(hyp_stages_organ_frame), ~hyp_sex + hyp_full1)
mcols(dds_hyp_all) <- cbind(mcols(dds_hyp_all), row.names(cs_hyp.count.tsv$counts))
rownames(dds_hyp_all) <- row.names(cs_hyp.count.tsv$counts)
#Run DESeq
dds_hyp_all <- DESeq(dds_hyp_all)
vst_dds_hyp_all <- vst(dds_hyp_all)
#Generate PCA
pcaData_hyp_all<- plotPCA(vst_dds_hyp_all,intgroup=c("hyp_stages","hyp_organs"), ntop=300, returnData=TRUE)
ggplot(pcaData_hyp_all, aes(x = PC1, y = PC2, color = factor(hyp_stages))) +
  geom_point(size=2)+
  theme_bw()r
#Look at distances if desired
hypsampleDists <- dist(t(assay(vst_dds_hyp_all)))
hypsampleDistMatrix <- as.matrix(hypsampleDists)
colnames(hypsampleDistMatrix) <- NULL
##make the heatmap
pheatmap(hypsampleDistMatrix, clustering_distance_rows=hypsampleDists,
         clustering_distance_cols = hypsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
hyptopVarGenes <- head( order( rowVars( assay(vst_dds_hyp_all) ), decreasing=TRUE ), 40 )
heatmap.2( assay(vst_dds_hyp_all)[ hyptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

###Quick brain mass phenotype
brainmass <- c(0.2990,	0.2590,	0.2650,	0.3150,	0.2620, 0.2433,0.2305,0.2316,0.2130,0.2002,0.2482,0.2410,0.3786,0.2731,0.2727,0.2627,0.2710,0.2347,0.2614,0.2527, 0.2488)
brainstage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",4),rep("Stage4",4),rep("Stage5",5))
brainmass_df <- data.frame(brainmass,brainstage)
ggplot(brainmass_df,aes(x=brainstage,y=brainmass))+
  geom_boxplot()+
  scale_y_continuous(name="brainMass", limits=c(0.1,.4))+
  theme_bw()

#hypothalamus DESeq analysis
#hypothalamus 2-4
lfc <- 0
hyp24res <- results(dds_hyp_all, contrast = c("hyp_full1","Stage4","Stage2"))
hyp24resSig <- subset(hyp24res,hyp24res$padj<.05)
hyp24resSigLog <- subset(hyp24resSig,abs(hyp24resSig$log2FoldChange)>=lfc)
hyp24up <- subset(hyp24resSigLog,(hyp24resSigLog$log2FoldChange)>=0)
hyp24down <- subset(hyp24resSigLog,(hyp24resSigLog$log2FoldChange)<=0)
DESeq2::plotMA(hyp24res, ylim = c (-4,4)); #drawLines()
as.data.frame(rownames(hyp24down))
rownames(hyp24down)
as.data.frame(rownames(hyp24up))
Reduce(intersect, list(row.names(hyp24resSig),row.names(hip24resSig),row.names(cor13resSig)))
length(hyp24up$log2FoldChange)
length(hyp24down$log2FoldChange)
#write.table(as.data.frame(rownames(hyp24down)), file='../analysis/DESeq2/hyp24down', quote=FALSE, sep='\t')
#write.table(as.data.frame(rownames(hyp24up)), file='../analysis/DESeq2/hyp24up', quote=FALSE, sep='\t')
subset(hyp24resSigLog,(rownames(hyp24resSig))=="GADD45B")
plotCounts(dds_hyp_all, gene="FOXO4", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="BCL2L1", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="FOS", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="NFKBIA", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="CTSK", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="FANCD2", intgroup="hyp_full1")

###RUN LISTS THROUGH KEGG online, can then make figures below
hyp24_kegg <- read_delim("../analysis/DESeq2/KEGG_workable.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
hyp24_kegg
hyp24_kegg <- subset(hyp24_kegg, PValue <= .05)
library(dplyr)
library(ggplot2)

hyp24_kegg$Term <- factor(hyp24_kegg$Term, levels = hyp24_kegg$Term[order(-hyp24_kegg$PValue)])

ggplot(hyp24_kegg , aes(Term,-log(PValue),color = Count))+
  geom_point(size=6.5)+
  scale_color_viridis_c(option="A", name="continuous\nvalue")+
  coord_flip()+
  theme_bw()

ggplot(hyp24_kegg, aes(x=Term,y=Count))+
  geom_point(stat = 'identity', linewidth=1)+
  coord_flip()+
  scale_fill_manual(values=c("#FE4F49","#5059FF","#5059FF","#5059FF", "#5059FF"))+
  scale_y_continuous(breaks = seq(-5, 5, len = 6))+
  theme_bw()

#Apoptosis gene graph
a_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
#try it with z-score on same graph in red
a_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
GADD45B_hyp <-a_hyp[rownames(a_hyp) %in% "GADD45B", ]
BCL2L1_hyp <-a_hyp[rownames(a_hyp) %in% "BCL2L1", ]
NFKBIA_hyp <-a_hyp[rownames(a_hyp) %in% "NFKBIA", ]
FOS_hyp <-a_hyp[rownames(a_hyp) %in% "FOS", ]
CALM3_hyp <-a_hyp[rownames(a_hyp) %in% "CALM3", ]
BDKRB1_hyp <-a_hyp[rownames(a_hyp) %in% "BDKRB1", ]
HMOX1_hyp <-a_hyp[rownames(a_hyp) %in% "HMOX1", ]
ZBTB16_hyp <-a_hyp[rownames(a_hyp) %in% "ZBTB16", ]
MAX_hyp <-a_hyp[rownames(a_hyp) %in% "MAX", ]
c_hyp  <- rep(c(c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))),9)
antiapop <- rbind(GADD45B_hyp,NFKBIA_hyp,BCL2L1_hyp,FOS_hyp,CALM3_hyp,BDKRB1_hyp,HMOX1_hyp,ZBTB16_hyp,MAX_hyp)
antiapop
z <- matrix(0,nrow(antiapop),ncol(antiapop))
z
for (i in 1:nrow(antiapop)) {
  for (j in 1:ncol(antiapop)) {
    z[i,j] <- (antiapop[i,j]-mean(antiapop[i,]))/sd(antiapop[i,])
  }
}
z
zz  <- matrix(0,nrow(z),5)
for (i in 1:nrow(z)) {
  zz[i,1] <- sum(z[i,1],z[i,2],z[i,3],z[i,4],z[i,5])/5
  zz[i,2] <- sum(z[i,6],z[i,7],z[i,8])/3
  zz[i,3] <- sum(z[i,9],z[i,10],z[i,11],z[i,12],z[i,13])/5
  zz[i,4] <- sum(z[i,14],z[i,15],z[i,16],z[i,17],z[i,18])/5
  zz[i,5] <- sum(z[i,19],z[i,20],z[i,21],z[i,22],z[i,23])/5
}
zz
real_antiapop <- as.matrix(rbind(as.matrix(zz[1,]),as.matrix(zz[2,]),as.matrix(zz[3,]),as.matrix(zz[4,]),as.matrix(zz[5,]),as.matrix(zz[6,]),as.matrix(zz[7,]),as.matrix(zz[8,]),as.matrix(zz[9,])))
close <-as.matrix(c(rep("GADD45B",5),rep("NFKBIA",5),rep("BCL2L1",5),rep("FOS",5),rep("CALM3",5),rep("BDKRB1",5),rep("HMOX1",5),rep("ZBTB16",5),rep("MAX",5)))
vclose<-as.matrix(c(rep("NoTC",20),rep("C",25)))
soclose<- rep(c("Stg1","Stg2","Stg3","Stg4","Stg5"),9)
bcl<-as.matrix(c(rep("NoTB",10),rep("YesB",5),rep("NoTB",30)))
realreal <-  cbind(as.data.frame(real_antiapop),as.data.frame(close),as.data.frame(soclose),as.data.frame(vclose),as.data.frame(bcl))
colnames(realreal) <- c("Expr","Gene","Stage","Cancer","BCL")
realreal
realreal <- realreal %>% 
  arrange(BCL == "YesB")
ggplot(realreal) +
  geom_line(data = subset(realreal, BCL != "YesB"),aes(x = Stage, colour = Cancer, y =Expr,group=Gene,linetype = BCL),size=1.5)+
  scale_color_manual(values=c("black","#FE4F49"))+
  geom_line(
    data = subset(realreal, BCL == "YesB"),
    aes(x = Stage, y = Expr, colour = Cancer, group = Gene, linetype = BCL),
    size = 1.5
  ) +
  theme_bw()


####MAKE VOLCANO PLOT FOR HYPTHALAMUS
hyp24resX <-  results(dds_hyp_all, contrast = c("hyp_full1","Stage4","Stage2"))
for (i in 1:length(hyp24resX$padj)) {
  if  (hyp24resX$padj[i]<1e-4 & !is.na (hyp24resX$padj[i])) {
    hyp24resX$padj[i] <- 1e-4
  }
  if (hyp24resX$log2FoldChange[i]>6 & !is.na (hyp24resX$log2FoldChange[i])) {
    hyp24resX$log2FoldChange[i] <- 6
  }
  if (hyp24resX$log2FoldChange[i]< -6 & !is.na (hyp24resX$log2FoldChange[i])) {
    hyp24resX$log2FoldChange[i] <- -6
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(hyp24resX$log2FoldChange) == 6, 17,
  ifelse(hyp24resX$padj==1e-4, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- '<log1e-15'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  hyp24resX$padj > 0.05, 'grey',
  ifelse(hyp24resX$log2FoldChange <= -1.58, 'red',
         ifelse(hyp24resX$log2FoldChange >= 1.58, 'blue',
                ifelse(hyp24resX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals
keyvals[is.na(keyvals)] <- 'grey'
str(keyvals)
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
length(keyvals)
length(keyvals.shape)
####
EnhancedVolcano(hyp24resX,
                lab = rownames(hyp24resX),
                xlim=c(-6 ,6),
                ylim=c(0,4.5),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = c('BCL2L1','GADD45B','NFKBIA','FOS','LMX1A', 'CCDC22', 'KCNS1', 'PAQR4'),
                pCutoff = .05,
                FCcutoff = 10,
                colCustom = keyvals,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5)

###Timeseq Hypothalamus with ABS LFC 0.5 change)
con_hyp3 <- c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
hyp3timeseq <- data.frame(sampleID = 1:23, group = c(1, 1, 1,1,1,2, 2, 2, 3, 3, 3,3,3,4, 4, 4, 4,4,5, 5, 5, 5,5),
                          timepoint = con_hyp3)
gf <- data.frame(chr = c(rep('chr1', 15296), rep('chr2', 2000), rep('chr4', 2000)),
                 start = rep(100, 19296),
                 end = rep(500, 19296),
                 id = row.names(cs_hyp.tpm.tsv$abundance))
library(TCseq)
hyp3tca <- TCA(design = hyp3timeseq, counts = round(DESeq2::counts(dds_hyp_all, normalized=TRUE)), gf)
hyp3tca <- DBanalysis(hyp3tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
hyp3tca <- timecourseTable(hyp3tca, value = "expression",  lib.norm = FALSE, filter = TRUE,abs.fold = 0.5)
hyp3_t <- tcTable(hyp3tca)
length(hyp3_t[,1])
#Gene counts: 19k to 14377 to 786
#Determine how many clusters using timeclust2 with clusGap
clusGap(hyp3_t,
        FUNcluster = timeclust2,
        algo = 'cm',
        K.max = 20,
        B = 20)
##12 for Hypothalamus .5LFC
chyp3tca <- timeclust(hyp3tca, algo = "cm", k = 12, standardize = TRUE)
XXXchyp3_px <-timeclustplot(chyp3tca, value = "z-score(PRKM)", cols = 2,cl.color = "gray50",membership.color= hcl.colors(30, "Viridis",rev = TRUE))
hyp3cxxx<-clustResults(chyp3tca)
hyp_clusters<-as.data.frame(hyp3cxxx@cluster)
write.table(hyp3cxxx@membership, "../analysis/TimeSeq/MemberShip", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#print individual cluster figures
rownames(subset(hyp_clusters,`hyp3cxxx@cluster` == 2))
hyp3cqqq <- 1
print(XXXchyp3_px[[hyp3cqqq]])

#Quantifying developmental vs seasonal
cluster <- read_table("../analysis/TimeSeq/cluster.tsv")
hyp3_t_combo <- merge(cluster, hyp3_t, by.x = "Gene", by.y = "row.names")
cluster1 <- subset(hyp3_t_combo,Cluster==1)
cluster2 <- subset(hyp3_t_combo,Cluster==2)
cluster3 <- subset(hyp3_t_combo,Cluster==3)
cluster4 <- subset(hyp3_t_combo,Cluster==4)
cluster5 <- subset(hyp3_t_combo,Cluster==5)
cluster6 <- subset(hyp3_t_combo,Cluster==6)
cluster7 <- subset(hyp3_t_combo,Cluster==7)
cluster8 <- subset(hyp3_t_combo,Cluster==8)
cluster9 <- subset(hyp3_t_combo,Cluster==9)
cluster10 <- subset(hyp3_t_combo,Cluster==10)
cluster11 <- subset(hyp3_t_combo,Cluster==11)
cluster12 <- subset(hyp3_t_combo,Cluster==12)
#pull zscore
zcluster1 <- t(apply(cluster1[,3:7], 1, function(x) scale(x)))
zcluster2 <- t(apply(cluster2[,3:7], 1, function(x) scale(x)))
zcluster3 <- t(apply(cluster3[,3:7], 1, function(x) scale(x)))
zcluster4 <- t(apply(cluster4[,3:7], 1, function(x) scale(x)))
zcluster5 <- t(apply(cluster5[,3:7], 1, function(x) scale(x)))
zcluster6 <- t(apply(cluster6[,3:7], 1, function(x) scale(x)))
zcluster7 <- t(apply(cluster7[,3:7], 1, function(x) scale(x)))
zcluster8 <- t(apply(cluster8[,3:7], 1, function(x) scale(x)))
zcluster9 <- t(apply(cluster9[,3:7], 1, function(x) scale(x)))
zcluster10 <- t(apply(cluster10[,3:7], 1, function(x) scale(x)))
zcluster11 <- t(apply(cluster11[,3:7], 1, function(x) scale(x)))
zcluster12 <- t(apply(cluster12[,3:7], 1, function(x) scale(x)))

#plot to check
library(tidyr)
# Convert to a data frame
zcluster12 <- as.data.frame(zcluster12)

# Restore column names (assuming they correspond to Stg1, Stg2, Stg3, etc.)
colnames(zcluster12) <- colnames(cluster12[, 3:7])

# Restore row names as a new column for gene names
zcluster12$Gene <- rownames(cluster12)
zcluster12<-as.data.frame(zcluster12)
zcluster1_long <- zcluster1 %>%
  pivot_longer(cols = starts_with("Stg"), 
             names_to = "Stage", 
             values_to = "Expression")

# Create the plot
ggplot(zcluster1_long, aes(x = Stage, y = Expression, group = Gene, color = Gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Gene Expression Across Stages",
    x = "Stage",
    y = "Expression (Z-score or raw)"
  ) +
  theme(legend.title = element_blank())

###mean zscores for all genes
zclus1_means <- colMeans(zcluster1[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus2_means <- colMeans(zcluster2[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus3_means <- colMeans(zcluster3[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus4_means <- colMeans(zcluster4[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus5_means <- colMeans(zcluster5[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus6_means <- colMeans(zcluster6[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus7_means <- colMeans(zcluster7[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus8_means <- colMeans(zcluster8[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus9_means <- colMeans(zcluster9[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus10_means <- colMeans(zcluster10[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus11_means <- colMeans(zcluster11[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])
zclus12_means <- colMeans(zcluster12[, c("Stg1", "Stg2", "Stg3", "Stg4", "Stg5")])

total_mean<-as.data.frame(rbind(zclus1_means,zclus2_means,zclus3_means,zclus4_means,zclus5_means,zclus6_means,zclus7_means,zclus8_means,zclus9_means,zclus10_means,zclus11_means,zclus12_means))
total_mean$Cluster <- c(1,2,3,4,5,6,7,8,9,10,11,12)
subset(total_mean, abs(Stg1-Stg5) > 1.25 & (abs(Stg3-Stg5) < .5 ))
subset(total_mean, abs(Stg1-Stg5) > 1.25 & (abs(Stg3-Stg5) < .5 | abs(Stg2-Stg3)<.5 ))
subset(total_mean, abs(Stg1-Stg5) > 1.25 & abs(Stg2-Stg4) < .75)
subset(total_mean, abs(Stg1-Stg5) > 1.25 & (abs(Stg3-Stg5) < .25 | abs(Stg3-Stg1) < .25))
#########################
#KEGG for TCSEQ figures
library(readr)
total_workaable <- read_delim("../analysis/TimeSeq/DAVID_workable2.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
total_workaable <- subset(total_workaable, PValue <= .05)
library(dplyr)
library(ggplot2)

total_workaable$Term <- factor(total_workaable$Term, levels = total_workaable$Term[order(-total_workaable$PValue)])
ggplot(total_workaable , aes(Term,-log(PValue),color = Count))+
  geom_point(size=9)+
  scale_color_viridis_c(option="A", name="continuous\nvalue")+
  coord_flip()+
  theme_bw()

#figure 4
interesting_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
purp_dot <- c("N","N","N","Y","N")
#CCDC22
plotCounts(dds_hyp_all, gene="CCDC22", intgroup="hyp_full1")
CCDC22_hyp <-as.data.frame(interesting_hyp[rownames(interesting_hyp) %in% "CCDC22", ])
CCDC22_hyp$Stage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(CCDC22_hyp) <- c("Expression", "Stage")
#plot
summary_CCDC22 <- CCDC22_hyp %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_CCDC22$purp <- purp_dot
# Step 3: Plot the data
CCDC22_plot<-ggplot(summary_CCDC22 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "black", size =2) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  geom_point(size = 9, aes(color = purp_dot)) +
  theme_classic()+
  scale_color_manual(values = c("black", "#8470FF"))
CCDC22_plot
ggsave("../analysis/DESeq2CCDC22_plot.png", CCDC22_plot,width = 2.5, height = 4, dpi =300,)

#KCNS1
plotCounts(dds_hyp_all, gene="KCNS1", intgroup="hyp_full1")
KCNS1_hyp <-as.data.frame(interesting_hyp[rownames(interesting_hyp) %in% "KCNS1", ])
KCNS1_hyp$Stage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(KCNS1_hyp) <- c("Expression", "Stage")
#plot
summary_KCNS1 <- KCNS1_hyp %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_KCNS1$purp <- purp_dot
# Step 3: Plot the data
KCNS1_plot<-ggplot(summary_KCNS1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "black", size =2) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  geom_point(size = 9, aes(color = purp_dot)) +
  theme_classic()+
  scale_color_manual(values = c("black", "#8470FF"))
KCNS1_plot
ggsave("../analysis/DESeq2KCNS1_plot.png", KCNS1_plot,width = 2.5, height = 4, dpi =300,)

#LMX1A
plotCounts(dds_hyp_all, gene="LMX1A", intgroup="hyp_full1")
LMX1A_hyp <-as.data.frame(interesting_hyp[rownames(interesting_hyp) %in% "LMX1A", ])
LMX1A_hyp$Stage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(LMX1A_hyp) <- c("Expression", "Stage")
#plot
summary_LMX1A <- LMX1A_hyp %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_LMX1A$purp <- purp_dot
# Step 3: Plot the data
LMX1A_plot<-ggplot(summary_LMX1A , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "black", size =2) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  geom_point(size = 9, aes(color = purp_dot)) +
  theme_classic()+
  scale_color_manual(values = c("black", "#8470FF"))
LMX1A_plot
ggsave("/Users/bill/ShrewProjects/Sorex_Hypothalamus_Transcriptomics2/analysis/DESeq2LMX1A_plot.png", LMX1A_plot,width = 2.5, height = 4, dpi =300,)

#PAQR4 #8470FF #3F3FFF #FE4F49
plotCounts(dds_hyp_all, gene="PAQR4", intgroup="hyp_full1")
PAQR4_hyp <-as.data.frame(interesting_hyp[rownames(interesting_hyp) %in% "PAQR4", ])
PAQR4_hyp$Stage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(PAQR4_hyp) <- c("Expression", "Stage")
#plot
summary_PAQR4 <- PAQR4_hyp %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_PAQR4$purp <- purp_dot
# Step 3: Plot the data
PAQR4_plot<-ggplot(summary_PAQR4 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "black", size =2) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  geom_point(size = 9, aes(color = purp_dot)) +
  theme_classic()+
  scale_color_manual(values = c("black", "#8470FF"))
PAQR4_plot
ggsave("../analysis/DESeq2PAQR4_plot.png", PAQR4_plot,width = 2.5, height = 4, dpi =300,)
