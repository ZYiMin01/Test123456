library(Seurat)   
library(ggplot2)
library(Matrix)
library(reshape2)
library(dplyr)
library(plyr)
library(rlang)
library(pheatmap)
library(SingleR)
library(clusterProfiler)

source("./fun/runseuratsingler.R")
source("./fun/scStatistics.R")
source("./fun/seuratCombination.R")
source("./fun/celltypeanno.R")


args=commandArgs(T)
dataPath0 <- args[1]
savePath <- args[2]
##BD file name like  *_RSEC_MolsPerCell
samples <-  as.vector(strsplit(args[3],',')[[1]]) 
species <- args[4]    #### human  or mouse
combmethod <- args[5]   ####"harmony" (one sample use None)
pctmt <- as.numeric(args[6])
coef <- as.numeric(args[7])
annodata <- as.vector(strsplit(args[8],',')[[1]]) 
anno_filter <-  as.vector(strsplit(args[9],',')[[1]])

# Basic Quality Control (QC)
for (name in samples){
    dataPath <- paste(dataPath0,"/",name,"/outs/",sep='')    

    stat.results <- runScStatistics(
        datatype = datatype,
        dataPath = dataPath,
        savePath = savePath,
        sampleName = name,
        species = species,
        coef = coef,
        genReport = F
    )
    # Single cell RNA-seq analysis using Seurat
	Seurat.result <- runSeurat2(
	dataPath, savePath, savePath,
	datatype = datatype,
	sampleName = name,
	species = species,
	anno.filter = anno_filter,
	nCell.min = 3, 
	bgPercent.max = 1,
	hg.mm.mix = F,
	pct.mt= pctmt,
	vars.to.regress = NULL
)

}
# ScRNA-seq Dataset Integration
combin.data <- runScCombination(samples, savePath,
               comb.method = combmethod,
               vars.to.regress = NULL,
               pc.use = 30,
               features = 2000
               )

# Run diffgene
diff.expr.genes <- runDiffExpr(combin.data,savePath,"all",n.markers=3)


# Celltype anno
annocelltype(combin.data,annodata,species=species)