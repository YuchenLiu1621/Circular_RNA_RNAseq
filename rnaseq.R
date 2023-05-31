#-----------------------------------------------------------
# File: Circular_RNA_RNA-seq.R
# Author: https://github.com/YuchenLiu1621
# Created: 2023-04-14
# Last Modified: 2023-05-31
# Description: This script performs RNA-seq analysis for circular RNA data. It includes the generation of volcano plots and scatter plots.
#-----------------------------------------------------------

#-----------------------------------------------------------
# 01. Preparing

# clean
rm(list=ls()) ##empty environment variable
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE="en")

# library
library(tidyverse)
library(stringr)
library(dplyr) 
library(corrplot) 
library(pheatmap)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(limma)
library(dplyr)
library(biomaRt)

p1 <- '/Path/to/file/' 
file_count = 'raw.anno.csv'
file_fpkm = 'fpkm.anno.csv'
gen1 <- 'NF2'

setwd(p1); getwd()
count1 <- read.csv(file_count)
fpkm1 <- read.csv(file_fpkm)

#-----------------------------------------------------------
# 02. Different Expression

DESeq.res <- function(arg.col1,arg.col2,arg.condition){
  exprSet <- count1[,c(arg.col1,arg.col2)]
  exprSet <- exprSet[rowSums(exprSet)!=0,]
  
  countData <- exprSet
  condition <- factor(arg.condition)
  
  colData <- data.frame(row.names=colnames(countData), condition)
  dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition )
  head(dds) 
  dds <- DESeq(dds)
  dds2 <- dds
  resultsNames(dds2)
  res <- results(dds2)
  
  name1 <- resultsNames(dds2)[2]
  out <- list(res = res, dds2 = dds2)
  return(out)
  }

DESeq.analysis <- function(arg.out){
  res <- arg.out$res
  dds2 <- arg.out$dds2
  
  res <- res[order(res$padj),]
  diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  diff_gene_deseq2 <- row.names(diff_gene_deseq2)
  resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
  DEG <- resdata
  
  DEG$change.padj = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1,ifelse(DEG$log2FoldChange > 1 ,'upregulated','downregulated'),'not_change'))
  return(DEG)
  }

NTvsPRE <- DESeq.res(7:9,4:6,c(rep("PRE",3),rep('NT',3))) 
NTvsL <- DESeq.res(1:3,4:6,c(rep("L1",3),rep('NT',3))) 
PREvsL <- DESeq.res(1:3,7:9,c(rep("L1",3),rep('PRE',3))) 

condition_NT_vs_PRE <- DESeq.analysis(NTvsPRE)
condition_NT_vs_L <- DESeq.analysis(NTvsL)
condition_PRE_vs_L <- DESeq.analysis(PREvsL)

#-----------------------------------------------------------
# 03. Heatmap

heatmap.res <- function(arg.condition,arg.color1,arg.title,arg.color2){
  dat1 <- na.omit(arg.condition)
  dat1$label <- ifelse(dat1$Row.names == gen1,gen1, "")
  label_data <- filter(dat1, Row.names == gen1)
  
  p <- ggplot(dat1, aes(x=log2FoldChange, y=-log10(padj),color=change.padj)) +
    scale_color_manual(values=arg.color) +
    geom_point(alpha=1, size=0.1) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold", angle = 90)) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2 (Fold Change)",
         y="-log10 (q value)",
         fontface = "bold") +
    ggtitle(arg.title) +
    geom_point(size = 2, alpha=1, color = arg.color2,data = label_data)
  
  return(p)
}

heatmap_PRE_vs_L <- heatmap.res(PRE_vs_L,c('#0942FD', '#D2DAE2','#ff2D04'),'PRE_vs_L1','#0942FD')
heatmap_NT_vs_L <- heatmap.res(NT_vs_L,c('#0942FD', '#D2DAE2','#ff2D04'),'NT_vs_L','#0942FD')
heatmap_NT_vs_PRE <- heatmap.res(NT_vs_PRE,'#D2DAE2','NT_vs_PRE','#D2DAE2')

#-----------------------------------------------------------
# 04. Scatter plot

dat1 <- fpkm1
dat1$gene <- rownames(dat1)
dat1$gene <- factor(dat1$gene)
highlight_gene <- dat1[dat1$gene == gen1, ]

gene_labels_nl <- data.frame(
  gene = gen1,
  label = gen1,
  x = highlight_gene$NT + 3,
  y = highlight_gene$L1 - 4
)

gene_labels_np <- data.frame(
  gene = gen1,
  label = gen1,
  x = highlight_gene$NT + 3,
  y = highlight_gene$PRE - 4
)

gene_labels_pl <- data.frame(
  gene = gen1,
  label = gen1,
  x = highlight_gene$PRE + 3,
  y = highlight_gene$L1 - 4
)

model.lm.nl <- lm(NT ~ L1, data = fpkm1[is.finite(rowSums(fpkm1)),])
model.lm.np <- lm(NT ~ PRE, data = fpkm1[is.finite(rowSums(fpkm1)),])
model.lm.pl <- lm(PRE ~ L1, data = fpkm1[is.finite(rowSums(fpkm1)),])

r2.nl = format(summary(model.lm.nl)$r.squared, digits = 4)
r2.np = format(summary(model.lm.np)$r.squared, digits = 4)
r2.pl = format(summary(model.lm.pl)$r.squared, digits = 4)

scatterplot_NT_vs_L <- ggplot(dat1, aes(x = NT, y = L1)) +
  geom_point(size = 0.1) +
  geom_point(
    data = highlight_gene,
    aes(x = NT, y = L1),
    color = "red",
    size = 1
  ) +
  geom_text(
    data = gene_labels_nl,
    aes(x = x, y = y, label = label),
    color = "black",
    size = 3,
    fontface = "bold.italic"
  ) +
  geom_text(
    aes(x = 0, y = 10, label = paste0("R^2=",r2.nl)),
    color = "black",
    size = 2.5,
    fontface = "bold"
  ) +
  labs(x=paste0('Log2(FPKM),NT'),
       y=paste0('Log2(FPKM),L1'))

scatterplot_NT_vs_PRE <- ggplot(dat1, aes(x = NT, y = PRE)) +
  geom_point(size = 0.1) +
  geom_point(
    data = highlight_gene,
    aes(x = NT, y = PRE),
    color = "red",
    size = 1
  ) +
  geom_text(
    data = gene_labels_np,
    aes(x = x, y = y, label = label),
    color = "black",
    size = 3,
    fontface = "bold.italic"
  ) +
  geom_text(
    aes(x = 0, y = 10, label = paste0("R^2=",r2.np)),
    color = "black",
    size = 2.5,
    fontface = "bold"
  ) +
  labs(x=paste0('Log2(FPKM),NT'),
       y=paste0('Log2(FPKM),PRE'))

scatterplot_PRE_vs_L <- ggplot(dat1, aes(x = PRE, y = L1)) +
  geom_point(size = 0.1) +
  geom_point(
    data = highlight_gene,
    aes(x = PRE, y = L1),
    color = "red",
    size = 1
  ) +
  geom_text(
    data = gene_labels_pl,
    aes(x = x, y = y, label = label),
    color = "black",
    size = 3,
    fontface = "bold.italic"
  ) +
  geom_text(
    aes(x = 0, y = 10, label = paste0("R^2=",r2.pl)),
    color = "black",
    size = 2.5,
    fontface = "bold"
  ) +
  labs(x=paste0('Log2(FPKM),PRE'),
       y=paste0('Log2(FPKM),L1'))


