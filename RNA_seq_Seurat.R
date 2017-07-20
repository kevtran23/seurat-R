#R code for analysis of single-cell RNA seq data

#load libraries
library(Seurat)
library(dplyr)
library(Matrix)

#set working directory
setwd("~/Downloads/Seurat_raw_data")

#read in data from different samples
TP335 <- read.table("TP335.txt",head=T)
TP336 <- read.table("TP336.txt",head=T)
TP328 <- read.table("TP328.txt",head=T)

#merge dataframes
mets <- cbind(TP335,TP336)
mets <- cbind(mets,TP328)

#create seurat object to analyze the data
object <- new("seurat", raw.data = mets)
object <- Setup(object, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = "Brain-Met Analysis")

#open pdf file for plotting
pdf("RNA_Seq_Plots.pdf")

#MeanVarPlot function calculates highly variable genes and loads it into the variable "object@var.genes"
object <- MeanVarPlot(object ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F, do.spike=T)
length(object@var.genes)     

#Perform PCA using genes that were highly dispersed(object@var.genes)
object <- PCA(object,pc.genes = object@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
object <- ProjectPCA(object)

#plot out results of PCA for visualization
PCA(object, pcs.print = 1:5, genes.print = 5, use.full = TRUE)

#use the results of the important principle components to group cells into clusters(that may represent the same cell type)
object <-FindClusters(object, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T,k.param = 10)

#perform tsNE to visualize the cell clusters made
object <- RunTSNE(object,dims.use = 1:10,do.fast=T,perplexity = 25)
TSNEPlot(object)

#Create a heatmap to visualize the genes that are most enriched in each cluster
object.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object, genes.use = top10$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)

#close pdfFile
dev.off()



