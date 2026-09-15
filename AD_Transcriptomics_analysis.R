#Downstream analysis of the Rush_AD dataset to find differentially expressed regular genes

#Loading the required libraries for differential expression analysis 
library(here)
library(DESeq2)
library(GEOquery)
library(pheatmap)
library(tidyverse)
library(tidyr)
library(stringr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#***************************************************
#STEP 1- Pre-processing
#***************************************************

#Loading the counts file 
counts <- read.csv(here("Data","cleaned_counts.csv"))
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
dim(counts)

#Extracting the metadata
data <- read.csv(here("Data","experiment_report.tsv"),header=TRUE,sep="\t",skip=1)
head(data)
dim(data)


#Mapping ENCSR to ENCFF
mapping <- data.frame(Accession=data$Accession, Files=data$Files, Condition=data$Biosample.summary)
mapping <- mapping %>%
  dplyr::select(Accession, Files, Condition)
mapping <- mapping %>%
  separate_rows(Files, sep=",")
mapping$Files <- gsub("","",mapping$Files)
mapping$Files <- gsub("/files/","",mapping$Files)
mapping$Files <- gsub("/","",mapping$Files)
head(mapping)


#Mapping ENCFF files belonging to fastq files
ENCFF <- readLines(here("Data","ENCFF_ids.txt"))


#checking if all required ENCFF files are present
all(ENCFF %in% mapping$Files)


#keeping only the ENCFF fastq files
mapping2 <- mapping[mapping$Files %in% ENCFF, ]


#Describing the conditions
mapping2 <- mapping2 %>%
  mutate(Condition=str_extract_all(mapping2$Condition,"(?<=with ).*(?=;)"))
mapping2$Condition[mapping2$Condition == "character(0)"] <- "No cognitive impairment"


#creating the column data
colData <- data.frame(Files=mapping2$Files, Condition=unlist(mapping2$Condition))


#Performing QC
all(colData$Files%in%colnames(counts))
colData <- colData[match(colnames(counts), colData$Files), ] #reordering
all(colData$Files == colnames(counts))
length(colnames(counts))
rownames(colData) <- colData$Files
identical(rownames(colData), colnames(counts))
is.na(colData$Condition)

#***************************************************
#STEP 2- Differential Expression analysis
#***************************************************

#Performing DESeq2
table(colData$Condition)
colData$Condition <- gsub("Cognitive impairment, Alzheimer's disease",
                          "Alzheimer's disease, Cognitive impairment",
                          colData$Condition)
colData$Condition <- factor(colData$Condition)
colData$Condition <- factor(colData$Condition,
                            levels = c("No cognitive impairment",
                                       "mild cognitive impairment",
                                       "Cognitive impairment",
                                       "Alzheimer's disease",
                                       "Alzheimer's disease, Cognitive impairment"))

args(DESeqDataSetFromMatrix)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design = ~Condition)
dim(dds)
head(dds)


#Filtering genes
genes_before <- nrow(dds)
genes_before
dds <- dds[rowSums(counts(dds)) >= 10, ]
genes_after <- nrow(dds)
genes_after


#Runing DESeq2
dds_deseq <- DESeq(dds)
sizeFactors(dds_deseq)
results(dds_deseq)
resultsNames(dds_deseq)
colData(dds_deseq)

#Results
res_mild <- results(dds_deseq, name="Condition_mild.cognitive.impairment_vs_No.cognitive.impairment")
png(here("Figures","MA_mild.png"), width = 1200, height = 900, res = 150)
plotMA(res_mild)
dev.off()
sig_res_mild <- subset(res_mild, padj<0.05 & abs(log2FoldChange)>1)
sig_res_mild # 2 DEGs are observed

res_cog <-results(dds_deseq, name="Condition_Cognitive.impairment_vs_No.cognitive.impairment")
png(here("Figures","MA_cog.png"), width = 1200, height = 900, res = 150)
plotMA(res_cog)
dev.off()
sig_res_cog <- subset(res_cog, padj<0.05 & abs(log2FoldChange)>1)
sig_res_cog #5 DEGs are observed


res_Alz <- results(dds_deseq, name="Condition_Alzheimer.s.disease_vs_No.cognitive.impairment")
png(here("Figures","MA_Alz.png"), width = 1200, height = 900, res = 150)
plotMA(res_Alz)
dev.off()
sig_res_alz <- subset(res_Alz, padj<0.05 & abs(log2FoldChange)>1)
sig_res_alz #0 DEGs are observed

res_Alz_cog <- results(dds_deseq,name="Condition_Alzheimer.s.disease..Cognitive.impairment_vs_No.cognitive.impairment")
png(here("Figures","MA_Alz_cog.png"), width = 1200, height = 900, res = 150)
plotMA(res_Alz_cog)
dev.off()
sig_res_Alz_cog <- subset(res_Alz_cog, padj<0.05 & abs(log2FoldChange)>1)
sig_res_Alz_cog #1 DEG is observed


#***************************************************
#STEP 3- Visualization of DEGs
#***************************************************

#Visualization of DEGs using volcano plot

#For No cognitive impairment vs mild cognitive impairment
png(here("Figures","Volcano_mild.png"), width = 1200, height = 900, res = 150)
plot(res_mild$log2FoldChange, -log10(res_mild$padj),
     pch=20,
     main="Volcano plot for Control VS mild cognitive impairment group",
     xlab="Log2 Fold Change",
     ylab="-log10 of p-adj value",
     col=ifelse(res_mild$padj <0.05 & res_mild$log2FoldChange>1, "brown", "grey"))
grid()
dev.off()

#For No cognitive impairment vs cognitive impairment
png(here("Figures","Volcano_cog.png"), width = 1200, height = 900, res = 150)
plot(res_cog$log2FoldChange, -log10(res_cog$padj),
     pch=20,
     main="Volcano plot for Control VS cognitive impairment group",
     xlab="Log2 Fold Change",
     ylab="-log10 P-adj value",
     col=ifelse(res_cog$padj<0.05 & res_cog$log2FoldChange>1, "blue",
                ifelse(res_cog$padj<0.05 & res_cog$log2FoldChange < -1, "brown","grey")))
grid()
dev.off()

#For No cognitive impairment vs Alzheimer's disease
png(here("Figures","Volcano_Alz.png"), width = 1200, height = 900, res = 150)
plot(res_Alz$log2FoldChange, -log10(res_Alz$padj),
     pch=20,
     main="Volcano plot for Control VS Alzheimer's disease group",
     xlab="Log2 Fold Change",
     ylab="-log10 P-adj value",
     col=ifelse(res_Alz$padj<0.05 & res_Alz$log2FoldChange>1, "blue",
                ifelse(res_Alz$padj<0.05 & res_Alz$log2FoldChange < -1, "brown","grey")))
grid()
dev.off()

#For No cognitive impairment vs Alzheimer's disease & cognitive impairment
png(here("Figures","Volcano_Alz_cog.png"), width = 1200, height = 900, res = 150)
plot(res_Alz_cog$log2FoldChange, -log10(res_Alz_cog$padj),
     pch=20,
     main="Volcano plot for Control VS Alzheimer's disease \n& cognitive impairment group",
     xlab="Log2 Fold Change",
     ylab="-log10 P-adj value",
     col=ifelse(res_Alz_cog$padj<0.05 & res_Alz_cog$log2FoldChange>1, "blue",
                ifelse(res_Alz_cog$padj<0.05 & res_Alz_cog$log2FoldChange < -1, "brown","grey")))
grid()
dev.off()

#PCA plot
vsd <- vst(dds_deseq, blind=FALSE)
colData(vsd)
png(here("Figures","PCA.png"), width = 1200, height = 900, res = 150)
plotPCA(vsd, intgroup="Condition")
dev.off()


#***************************************************
#STEP 4- Ranking data for Enrichment analysis
#***************************************************

#Ranking genes

#For mild cognitive impariment vs control group
nrow(res_mild)
sum(!is.na(res_mild$stat))
res_mild_rank <- res_mild[!is.na(res_mild$stat),]
ranked_df_mild <- data.frame(gene=rownames(res_mild_rank),
                             stat=as.numeric(res_mild_rank$stat))
ranked_df_mild <- ranked_df_mild[order(ranked_df_mild$stat, decreasing=TRUE),]


#For cognitive impariment vs control group
res_cog_rank <- res_cog[!is.na(res_cog$stat),]
ranked_df_cog <- data.frame(gene=rownames(res_cog_rank),
                            stat=as.numeric(res_cog_rank$stat))
ranked_df_cog <- ranked_df_cog[order(ranked_df_cog$stat, decreasing=TRUE),]


#For Alzheimer's vs control group
res_Alz_rank <- res_Alz[!is.na(res_Alz$stat),]
ranked_df_Alz <- data.frame(gene=rownames(res_Alz_rank),
                            stat=as.numeric(res_Alz_rank$stat))
ranked_df_Alz <- ranked_df_Alz[order(ranked_df_Alz$stat, decreasing=TRUE),]


#For Alz & cognitive impariment vs control group
res_Alz_cog_rank <- res_Alz_cog[!is.na(res_Alz_cog$stat),]
ranked_df_Alz_cog <- data.frame(gene=rownames(res_Alz_cog_rank),
                                stat=as.numeric(res_Alz_cog_rank$stat))
ranked_df_Alz_cog <- ranked_df_Alz_cog[order(ranked_df_Alz_cog$stat, decreasing=TRUE),]


#Setting as names
rank_mild   <- setNames(ranked_df_mild$stat, ranked_df_mild$gene)
rank_cog    <- setNames(ranked_df_cog$stat, ranked_df_cog$gene)
rank_Alz    <- setNames(ranked_df_Alz$stat, ranked_df_Alz$gene)
rank_Alz_cog <- setNames(ranked_df_Alz_cog$stat, ranked_df_Alz_cog$gene)


#Formatting ENSEMBL IDs
names(rank_mild) <- sub("\\..*","",names(rank_mild))
names(rank_mild)

names(rank_cog) <- sub("\\..*","",names(rank_cog))
names(rank_cog)

names(rank_Alz) <- sub("\\..*","",names(rank_Alz))
names(rank_Alz)

names(rank_Alz_cog) <- sub("\\..*","",names(rank_Alz_cog))
names(rank_Alz_cog)


#checking for duplicates
sum(duplicated(names(rank_mild)))
sum(duplicated(names(rank_cog)))
sum(duplicated(names(rank_Alz)))
sum(duplicated(names(rank_Alz_cog)))


#*******************************************************************
#STEP 5- Gene Set Enrichment Analysis using ClusterProfile package
#*******************************************************************

#For mild cognitive impairment vs control group

set.seed(42)
#GO analysis
gse_mild <- gseGO(geneList     = rank_mild,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENSEMBL",        
                  ont          = "ALL",
                  minGSSize    = 15,
                  maxGSSize    = 800,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH")
gse_mild_df <- as.data.frame(gse_mild)
head(gse_mild_df)
dotplot(gse_mild, showCategory=15, font.size =7)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_go_mild.png"), width=8, height=8, dpi=300)

cnetplot(gse_mild, showCategory = 5, foldChange = rank_mild,  node_label = "category")
ggsave(here("Figures","cnetplot_go_mild.png"), width=8, height=8, dpi=300)

#KEGG analysis
#Converting EMSEMBL to ENTREZ
mild_ids <- bitr(names(rank_mild), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(mild_ids)
rank_mild_entrez <- rank_mild[mild_ids$ENSEMBL]
names(rank_mild_entrez) <- mild_ids$ENTREZID
rank_mild_entrez <- rank_mild_entrez[!duplicated(names(rank_mild_entrez))]
rank_mild_entrez <- sort(rank_mild_entrez, decreasing = TRUE)
length(rank_mild_entrez)

gse_kegg_mild <- gseKEGG(
  geneList      = rank_mild_entrez,
  organism      = "hsa",       
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  seed          = TRUE)

kegg_mild_df <- as.data.frame(gse_kegg_mild)
nrow(kegg_mild_df)
dotplot(gse_kegg_mild, showCategory = 10)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_kegg_mild.png"), width=8, height=8, dpi=300)



#For cognitive impairment vs control group

set.seed(42)
#GO analysis
gse_cog <- gseGO(geneList     = rank_cog,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "ENSEMBL",     
                 ont          = "ALL",
                 minGSSize    = 15,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
gse_cog_df <- as.data.frame(gse_cog)
head(gse_cog_df)
dotplot(gse_cog, showCategory=15, font.size=7)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_go_cog.png"), width=8, height=8, dpi=300)

cnetplot(gse_cog, showCategory=5, foldChange=rank_cog, node_label = "category")
ggsave(here("Figures","cnetplot_go_cog.png"), width=8, height=8, dpi=300)

#KEGG analysis
#Converting EMSEMBL to ENTREZ
cog_ids <- bitr(names(rank_cog), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(cog_ids)
rank_cog_entrez <- rank_cog[cog_ids$ENSEMBL]
names(rank_cog_entrez) <- cog_ids$ENTREZID
rank_cog_entrez <- rank_cog_entrez[!duplicated(names(rank_cog_entrez))]
rank_cog_entrez <- sort(rank_cog_entrez, decreasing = TRUE)
length(rank_cog_entrez)

gse_kegg_cog <- gseKEGG(
  geneList      = rank_cog_entrez,
  organism      = "hsa",       
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  seed          = TRUE)

kegg_cog_df <- as.data.frame(gse_kegg_cog)
nrow(kegg_cog_df)
dotplot(gse_kegg_cog, showCategory = 10)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_kegg_cog.png"), width=8, height=8, dpi=300)



#For Alzheimer's disease vs control group

set.seed(42)
#GO analysis
gse_Alz <- gseGO(geneList     = rank_Alz,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "ENSEMBL",     
                 ont          = "ALL",
                 minGSSize    = 15,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
gse_Alz_df <- as.data.frame(gse_Alz)
head(gse_Alz_df)
dotplot(gse_Alz, showCategory=15, font.size=7)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_go_Alz.png"), width=8, height=8, dpi=300)

cnetplot(gse_Alz, showCategory=5, foldChange=rank_Alz, node_label = "category")
ggsave(here("Figures","cnetplot_go_Alz.png"), width=8, height=8, dpi=300)

#KEGG analysis
#Converting EMSEMBL to ENTREZ
Alz_ids <- bitr(names(rank_Alz), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(Alz_ids)
rank_Alz_entrez <- rank_Alz[Alz_ids$ENSEMBL]
names(rank_Alz_entrez) <- Alz_ids$ENTREZID
rank_Alz_entrez <- rank_Alz_entrez[!duplicated(names(rank_Alz_entrez))]
rank_Alz_entrez <- sort(rank_Alz_entrez, decreasing = TRUE)
length(rank_Alz_entrez)

gse_kegg_Alz <- gseKEGG(
  geneList      = rank_Alz_entrez,
  organism      = "hsa",       
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  seed          = TRUE)

kegg_Alz_df <- as.data.frame(gse_kegg_Alz)
nrow(kegg_Alz_df)
dotplot(gse_kegg_Alz, showCategory = 10)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_kegg_Alz.png"), width=8, height=8, dpi=300)


#For Alzheimer's disease & cognitive impairment vs control group

set.seed(42)
#GO analysis
gse_Alz_cog <- gseGO(geneList     = rank_Alz_cog,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENSEMBL",     
                     ont          = "ALL",
                     minGSSize    = 15,
                     maxGSSize    = 800,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
gse_Alz_cog_df <- as.data.frame(gse_Alz_cog)
head(gse_Alz_cog_df)
dotplot(gse_Alz_cog, showCategory=15, font.size=7)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_go_Alz_cog.png"), width=8, height=8, dpi=300)

cnetplot(gse_Alz_cog, showCategory=5, foldChange=rank_Alz_cog, node_label = "category")
ggsave(here("Figures","cnetplot_go_Alz_cog.png"), width=8, height=8, dpi=300)


#KEGG analysis
#Converting EMSEMBL to ENTREZ
Alz_cog_ids <- bitr(names(rank_Alz_cog), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(Alz_cog_ids)
rank_Alz_cog_entrez <- rank_Alz_cog[Alz_cog_ids$ENSEMBL]
names(rank_Alz_cog_entrez) <- Alz_cog_ids$ENTREZID
rank_Alz_cog_entrez <- rank_Alz_cog_entrez[!duplicated(names(rank_Alz_cog_entrez))]
rank_Alz_cog_entrez <- sort(rank_Alz_cog_entrez, decreasing = TRUE)
length(rank_Alz_cog_entrez)

gse_kegg_Alz_cog <- gseKEGG(
  geneList      = rank_Alz_cog_entrez,
  organism      = "hsa",       
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0,
  seed          = TRUE)

kegg_Alz_cog_df <- as.data.frame(gse_kegg_Alz_cog)
nrow(kegg_Alz_cog_df)
dotplot(gse_kegg_Alz_cog, showCategory = 10)+ facet_grid(.~.sign)
ggsave(here("Figures","dotplot_kegg_Alz_cog.png"), width=8, height=8, dpi=300)
