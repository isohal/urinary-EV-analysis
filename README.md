# urinary-EV-analysis
This page describes the code used to conduct bioinformatics analysis in this publication.

## A cross-week analysis of urinary extracellular vesicles after respiratory tract exposure intervention in a cohort of isocyanate-exposed workers

Ikjot S. Sohal<sup>1*</sup>, Anoop Pal<sup>2</sup>, Jack Lepine<sup>3</sup>, Pengyuan Liu<sup>4</sup>, Adam Wisnewski<sup>5</sup>, Carrie A. Redlich<sup>5</sup>, Dhimiter Bello<sup>6*</sup>

<sup>1</sup>Department of Biological Sciences, College of Science
University of North Texas
Denton, TX 76201, USA

<sup>2</sup>Izon Science Ltd.
196 Boston Ave, Medford, MA 02155

<sup>3</sup>Next Generation Sequencing Core Facility
University of Massachusetts Lowell
Lowell, MA 01854 

<sup>4</sup>Department of Chemistry
University of Massachusetts Lowell
Lowell, MA 01854, USA

<sup>5</sup>Yale University
Department of Internal Medicine, 
Yale University School of Medicine
New Haven, CT 06510, USA

<sup>6</sup>Department of Biomedical and Nutritional Sciences, Zuckerberg College of Health Sciences
University of Massachusetts Lowell
Lowell, MA 01854, USA

\*To whom correspondence should be addressed:

Ikjot S. Sohal – ikjot.sohal@unt.edu

Dhimiter Bello – dhimiter_bello@uml.edu

## R script for EV-mRNA, EV-ncRNA, and EV-proteomics analysis
The following script describes the code to extract read count matrix and generate expression dataset using *DESeq2* package. The code is shown using EV-mRNA dataset as an example. 
``` r
library(readr)

#Extracting the count matrix
CountMatrix <- read.csv(file = "<filepath>/ExpMatrix_mRNAs.csv")
CountMatrix.Final = CountMatrix[,-1]
rownames(CountMatrix.Final) = CountMatrix$gene_id
colnames(CountMatrix.Final) = colnames(CountMatrix_MA[,-1])
CountMatrix.Final[is.na(CountMatrix.Final)] = 0

#Creating condition labels for samples in the normalized matrix
deseqData = cbind.data.frame(
  sample = c("S1_Pre_RESP", "S2_Pre_RESP", "S3_Pre_RESP", "S4_Pre_RESP","S5_Pre_RESP","S1_Post_RESP", "S2_Post_RESP", "S3_Post_RESP", "S4_Post_RESP", "S5_Post_RESP"),
  worker = c("S1", "S2", "S3", "S4", "S5", "S1", "S2", "S3", "S4", "S5"),
  condition = c("Pre_RESP", "Pre_RESP", "Pre_RESP", "Pre_RESP", "Pre_RESP", "Post_RESP", "Post_RESP", "Post_RESP", "Post_RESP", "Post_RESP"))
rownames(deseqData) = colnames(CountMatrix.Final)

#Confirming that colnames of normalized matrix is equal to colnames of "condition" vector
all(colnames(CountMatrix.Final) %in% rownames(deseqData))


library(DESeq2)
#Generating DESeq expression dataset from normalized count matrix
dds = DESeqDataSetFromMatrix(
  countData = CountMatrix.Final,
  colData = deseqData,
  design = ~ worker + condition)
colSums(counts(dds[,1:10]))

#Run DESeq
dds_mRNAs = DESeq(dds)

#Principal Component Analysis
library(ggplot2)
dds_mRNAs$condition = as.factor(c("Pre_RESP", "Pre_RESP", "Pre_RESP", "Pre_RESP", "Pre_RESP", "Post_RESP", "Post_RESP", "Post_RESP", "Post_RESP", "Post_RESP"))
vsd_mRNAs = varianceStabilizingTransformation(dds_mRNAs)
pca_mRNAs <- prcomp(t(assay(vsd_mRNAs)))
p_mRNAs = plotPCA(vsd_mRNAs, intgroup = c("condition"), returnData = T)

t.test(PC1 ~ group, data = p_mRNAs)
t.test(PC2 ~ group, data = p_mRNAs)

#Extracting results for pre-RESP vs post-RESP comparison
DEmRNAs_results = results(dds_mRNAs, contrast = c("condition", "Pre_RESP", "Post_RESP"))
DEmRNAs_results = DEmRNAs_results[!is.na(DEmRNAs_results_MvsF[,"pvalue"]),] #remove pvalues as NA
DEmRNAs_results = DEmRNAs_results[DEmRNAs_results[,"pvalue"] < 0.05 & abs(DEmRNAs_results$log2FoldChange) > 0.58,] #keep pvalues < 0.05 and log2FC > 0.58 and < -0.58
DEmRNAs_results_clean = DEmRNAs_results[!duplicated(DEmRNAs_results$hgnc_name) & !duplicated(DEmRNAs_results$hgnc_name, fromLast = TRUE),]

#Generate volcano plots for pre-RESP vs post-RESP comparison
library(EnhancedVolcano)

#Assign color values to upregulated, downregulated and non-significant mRNAs
keyvals = rep('#BBBBBB', nrow(DEmRNAs_results_clean))
names(keyvals) = rep('Non-significant', nrow(DEmRNAs_results_clean))
keyvals[which(DEmRNAs_results_clean$log2FoldChange > 0.58 & DEmRNAs_results_clean$pvalue < 0.05)] = '#EE6677'
names(keyvals)[which(DEmRNAs_results_clean$log2FoldChange > 0.58 & DEmRNAs_results_clean$pvalue < 0.05)] = 'Upregulated'
keyvals[which(DEmRNAs_results_clean$log2FoldChange < -0.58 & DEmRNAs_results_clean$pvalue < 0.05)] = '#4477AA'
names(keyvals)[which(DEmRNAs_results_clean$log2FoldChange < -0.58 & DEmRNAs_results_clean$pvalue < 0.05)] = 'Downregulated'

#Determining top 5 upregulated mRNAs
a = rownames(DEmRNAs_results_clean)[which(DEmRNAs_results_clean$log2FoldChange > 0.58 & DEmRNAs_results_clean$pvalue < 0.05)]
b = order(DEmRNAs_results_clean[c(a),]$pvalue, decreasing = F)
c = a[b[1:5]]

#Determining top 5 downregulated mRNAs
x = rownames(DEmRNAs_results_clean)[which(DEmRNAs_results_clean$log2FoldChange < -0.58 & DEmRNAs_results_clean$pvalue < 0.05)]
y = order(DEmRNAs_results_clean[c(x),]$pvalue, decreasing = F)
z = x[y[1:5]]

#Make volcano plot
EnhancedVolcano(DEmRNAs_results_clean,
                lab = rownames(DEmRNAs_results_clean),
                axisLabSize = 20,
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(1,5),
                xlim = c(-15,15),
                shape = 19,
                ylab = expression("-Log"[10]*" (p-value)"),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                pointSize = 4,
                pCutoff = 5e-2,
                cutoffLineWidth = 0.75,
                FCcutoff = 0,
                colAlpha = 0.5,
                legendLabSize = 12,
                colCustom = keyvals,
                gridlines.minor = F,
                border = 'full',
                borderWidth = 1,
                selectLab = rownames(DEmRNAs_results_clean[c(c(z), c(c)),]),
                labSize = 3,
                boxedLabels = T,
                drawConnectors = T,
                widthConnectors = 0.5,
                lengthConnectors = unit(0.01, 'npc'),
                typeConnectors = 'open')

#Plotting heatmap for differentially expressed mRNAs between Mon-AM and Fri-PM
library(pheatmap)
library(RColorBrewer)
library(grid)
d = cbind(CountMatrix.Final[c(DEmRNAs_results_clean$gene_id),])
colnames(d) = c("S1_Pre_RESP", "S2_Pre_RESP", "S3_Pre_RESP", "S4_Pre_RESP","S5_Pre_RESP","S1_Post_RESP", "S2_Post_RESP", "S3_Post_RESP", "S4_Post_RESP", "S5_Post_RESP")
rownames(d) = DEmRNAs_results_clean$gene_id

annotation_col = data.frame(Subjects = factor(c(rep("Pre_RESP", 5), rep("Post_RESP", 5))))
rownames(annotation_col) = colnames(d)

annotation_colors = list(Subjects = c(Pre_RESP = "red3",
                           Post_RESP = "gold2"))

par(mar = c(0.5,0.5,0.5,0.5), bg = "#F8F8F8")
ph = pheatmap(d,
         scale = "row",
         cluster_rows = T,
         clustering_distance_rows = "manhattan",
         cluster_cols = T,
         clustering_distance_cols = "manhattan",
         clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(8,'RdYlBu')))(256),
         show_rownames = F,
         show_colnames = T,
         kmeans_k = NA,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cutree_rows = 2,
         cutree_cols = 2,
         silent = F,
         fontsize = 11,
         legend = T,
         treeheight_row = 0,
         annotation_legend = T)
