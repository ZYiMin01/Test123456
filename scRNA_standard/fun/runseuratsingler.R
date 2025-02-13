filterCell <- function(cell.manifest, filter.thres){
    rownames(filter.thres) <- filter.thres$Index
    ix.cr <- which(cell.manifest$droplet.type == "cell")
    ix.thres <- getCellix(
        cell.manifest, filter.thres,
        arg = c("nUMI", "nGene", "mito.percent", "ribo.percent")
    )
    return(cell.manifest$barcodes[intersect(ix.cr, ix.thres)])
}



filterGene <- function(gene.manifest,
                       anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                       nCell.min = 3,
                       bgPercent.max = 1){
    ix.type <- which(!(gene.manifest$Annotation %in% anno.filter))
    ix.nCell <- which(gene.manifest$nCell >= nCell.min)
    ix <- intersect(ix.type, ix.nCell)

    if("bg.percent" %in% colnames(gene.manifest)){
        ix.bg <- which(gene.manifest$bg.percent <= bgPercent.max)
        ix <- intersect(ix, ix.bg)
    }else{
        if(bgPercent.max < 1){
            cat("- Warning in 'filterGene': Can not filter gene by 'bg.percent' due to the lack of background distribution.\n")
        }
    }
    return(gene.manifest$Symbol[ix])
}


rmContamination <- function(expr.data, cell.manifest, contamination.fraction){
    bg.bars <- subset(cell.manifest, nUMI >= 1 & nUMI <= 10)$barcodes
    bg.sum <- Matrix::rowSums(expr.data[, bg.bars])
    bg.percent <- bg.sum / sum(bg.sum)
    soupProfile <- data.frame(row.names = rownames(expr.data),
                         est = bg.percent,
                         counts = bg.sum)

    cell.manifest <- subset(cell.manifest, droplet.type == "cell")
    sc <- list(toc = expr.data[, cell.manifest$barcodes],
               metaData = data.frame(row.names = cell.manifest$barcodes,
                                     nUMIs = cell.manifest$nUMI,
                                     rho = contamination.fraction),
               soupProfile = soupProfile)
    class(sc) = c("list", "SoupChannel")
    out = adjustCounts(sc, roundToInt = TRUE)

    return(out)
}


#' runSeurat2
#'
#' According to the QC results of scStatistics, filter cells and genes.
#' Prepare a Seurat object.
#'
#' @param dataPath A path containing the cell ranger processed data.
#' Under this path, folders 'filtered_feature_bc_matrix' and 'raw_feature_bc_matrix' exist generally.
#' @param statPath A path containing the results files of step 'runScStatistics'.
#' @param savePath A path to save the results files. If NULL, the 'statPath' will be used instead.
#' @param datatype A character string giving a type of data.
#' @param sampleName A character string giving a label for this sample.
#' @param bool.filter.cell A logical value indicating whether to filter the cells
#' according to the QC of 'scStatistics'.
#' @param bool.filter.gene A logical value indicating whether to filter the genes
#' according to the QC of 'scStatistics'.
#' @param anno.filter A vector indicating the types of genes to be filtered.
#' Must be some of c("mitochondrial", "ribosome", "dissociation")(default) or NULL.
#' @param nCell.min An integer number used to filter gene. The default is 3.
#' Genes with the number of expressed cells less than this threshold will be filtered.
#' @param bgPercent.max A float number used to filter gene. The default is 1 (no filtering).
#' Genes with the background percentage larger than this threshold will be filtered.
#' @param bool.rmContamination A logical value indicating whether to remove ambient RNA contamination based on 'SoupX'.
#' @param contamination.fraction A float number between 0 and 1 indicating the estimated contamination fraction.
#' The default is NULL and the result of scStatistics will be used.
#' @param vars.add.meta A vector indicating the variables to be added to Seurat object's meta.data.
#' The default is c("mito.percent", "ribo.percent", "diss.percent").
#' @param vars.to.regress A vector indicating the variables to regress out in R package Seurat.
#' The default is c("nUMI", "mito.percent", "ribo.percent").
#'
#' @export
#'
#' @import Matrix  ggplot2 Seurat
#'

runSeurat2 <- function(dataPath, statPath, savePath,
                          datatype = "10X",
                          sampleName = "sc",species = "human",
                          bool.filter.cell = T,
                          bool.filter.gene = T,
                          anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                          nCell.min = 3, bgPercent.max = 1,
                          hg.mm.mix = F,pct.mt=0.25,
                          bool.rmContamination = T,
                          contamination.fraction = NULL,
                          vars.add.meta = c("mito.percent", "ribo.percent", "diss.percent"),
                          vars.to.regress = c("nUMI", "mito.percent", "ribo.percent")){
    if(datatype == "10X"){
        raw.data = F
        data.path <- get10Xpath(dataPath, raw.data = raw.data)
        if(is.null(data.path)){
            raw.data = F
            data.path <- get10Xpath(dataPath, raw.data = raw.data)
            if(is.null(data.path)){
                stop("Cannot find the raw data or filtered data.\n")
            }else{
                bool.rmContamination <- F
                cat("- Warning in 'prepareSeurat': Cannot find the raw data, and use the filtered data instead.\n")
            }
        }
        message("[", Sys.time(), "] -----: data preparation")
        expr.data <- Read10Xdata(data.dir = data.path, only.expr = T)
    }else{ 
        raw.data = F
        bool.rmContamination <- F
        message("[", Sys.time(), "] -----: data preparation")
        expr.data <- t(read.table(file = paste(dataPath,"/",sampleName,"_RSEC_MolsPerCell.csv",sep=""),row.names=1, header =TRUE,sep = ","))
        colnames(expr.data) <- paste(colnames(expr.data),"_",sampleName,sep="")
    }
    
    rownames(expr.data) <- gsub("_", "-", rownames(expr.data))

    gene.manifest <- read.table(file.path(statPath, paste(sampleName,'_qc/geneManifest.txt',sep="")),
                                header = T, sep = "\t", stringsAsFactors = F)
    rownames(gene.manifest) <- gene.manifest$Symbol

    cell.manifest <- read.table(file.path(statPath, paste(sampleName,'_qc/cellManifest-all.txt',sep="")),
                                header = T, stringsAsFactors = F)
    rownames(cell.manifest) <- cell.manifest$barcodes

    filter.thres <- read.table(file.path(statPath, paste(sampleName,'_qc/cell.QC.thres.txt',sep="")),
                               header = T, stringsAsFactors = F)
    filter.thres[3,3] = as.numeric(gsub(Inf, 0.5, filter.thres[3,3]))

    if(hg.mm.mix){
        rownames(expr.data) <- substr(rownames(expr.data), 6, 60)
    }

    if(bool.rmContamination & is.null(contamination.fraction)){
        if(file.exists(file.path(statPath, 'ambientRNA-SoupX.txt'))){
            soupX.file <- readLines(file.path(statPath, 'ambientRNA-SoupX.txt'))
            contamination.fraction <- as.numeric(soupX.file[length(soupX.file)])
        }else{
            bool.rmContamination <- F
            cat("- Warning in removing ambient RNAs contamination: Cannot find the estimated contamination fraction and skip this step.\n")
        }
    }

    if(bool.rmContamination){
        message("[", Sys.time(), "] -----: contamination removing")
        expr.data <- rmContamination(expr.data, cell.manifest, contamination.fraction)
    }

    if(bool.filter.cell){
        cells.select <- filterCell(cell.manifest,
                                   filter.thres = filter.thres)
    }else{
        cells.select <- subset(cell.manifest, droplet.type == "cell")$barcodes
    }
    if(bool.filter.gene){
        genes.select <- filterGene(gene.manifest,
                                   anno.filter = anno.filter,
                                   nCell.min = nCell.min,
                                   bgPercent.max = bgPercent.max)
    }else{
        genes.select <- gene.manifest$Symbol
    }
    cells.select <- intersect(cells.select,cell.manifest[cell.manifest[,'HB.percent']<0.01,'barcodes'])
    cells.select <- intersect(cells.select,cell.manifest[cell.manifest[,'mito.percent']< pct.mt,'barcodes'])

    ##scrublet doublets 
    doublets.file <- read.table(file.path(statPath, paste("scrublet/",sampleName,".doublet_prediction.txt",sep="")),sep=",",
                               header = T, stringsAsFactors = F)
    cells.select <- intersect(cells.select,doublets.file$barcodes[doublets.file$prediction == "False"])
    expr.data <- expr.data[genes.select, cells.select]   

    message("[", Sys.time(), "] -----: Seurat object creation")
    expr = CreateSeuratObject(counts = expr.data,
                              min.cells = 0,
                              min.features = 0,
                              project = sampleName)

    for(mv in vars.add.meta){
        if(!(mv %in% colnames(cell.manifest))){
            cat("- Warning in 'prepareSeurat': the", mv, " column is not in cell.manifest.\n")
        }else if(mv %in% colnames(expr@meta.data)){
            cat("- Warning in 'prepareSeurat': the", mv, " column has been in the Seurat object meta.data.\n")
        }else{
            tmp <- cell.manifest[cells.select, mv]
            names(tmp) <- cells.select
            expr[[mv]] <- tmp
        }
    }
    expr$stim <- sampleName
    saveRDS(expr, file = file.path(savePath, paste(sampleName,".expr.RDS",sep="")))
    expr <- NormalizeData(object = expr,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000,
                          verbose = F)

    message("[", Sys.time(), "] -----: highly variable genes")
    expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000, verbose = F)

    message("[", Sys.time(), "] -----: data scaling")
    expr <- ScaleData(object = expr,
                      vars.to.regress = vars.to.regress,
                      verbose = F)

    gene.select.type <- rep("filter", dim(gene.manifest)[1])
    names(gene.select.type) <- gene.manifest$Symbol
    gene.select.type[genes.select] <- "keep"
    names(gene.select.type) <- gsub("_", "-", names(gene.select.type))
    gene.select.type[VariableFeatures(expr)] <- "hvg"
    names(gene.select.type) <- gene.manifest$Symbol
    gene.manifest$Select <- gene.select.type

    top10 <- head(VariableFeatures(expr), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(expr)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    pdf(file=file.path(savePath, paste(sampleName,"_VariableFeaturePlot.pdf",sep="")),width=8,height=5,onefile=F)
    print(plot2)
    dev.off()
	png(file=file.path(savePath, 
	paste(sampleName,"_VariableFeaturePlot.png",sep="")),width=800,height=500)
    print(plot2)
    dev.off()
    ############cell cycle analysis
    s.genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
    g2m.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2", "G2E3", "GAS2L3","CBX5", "CENPA")
    if(species == "human"){
        s.genes <- s.genes
        g2m.genes <- g2m.genes
    }else if(species == "mouse" | species == "rat"){
        s.genes <- getMouseGene(s.genes)
        g2m.genes <- getMouseGene(g2m.genes)
    }
    expr <- CellCycleScoring(expr, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    write.table(expr[[]],file.path(savePath, paste(sampleName,"cell_cycle_score.csv",sep="_")),sep=',',quote=F,col.names=NA)
    #CellCycle gene
    expr <- RunPCA(expr, features = c(s.genes, g2m.genes))
    pdf(file=file.path(savePath, paste(sampleName, "cell_cycle_assessment_CellCyclegene.pdf",sep="_")),width=5,height=4,onefile=F)
    print(DimPlot(expr))
    dev.off()
	png(file=file.path(savePath, paste(sampleName, "cell_cycle_assessment_CellCyclegene.png",sep="_")),width=500,height=400)
    print(DimPlot(expr))
    dev.off()
    #all gene
    expr <- RunPCA(expr)
    pdf(file=file.path(savePath, paste(sampleName, "cell_cycle_assessment_allgene.pdf",sep="_")),width=5,height=4,onefile=F)
    print(DimPlot(expr))
    dev.off()    
    png(file=file.path(savePath, paste(sampleName, "cell_cycle_assessment_allgene.png",sep="_")),width=500,height=400)
    print(DimPlot(expr))
    dev.off()
    write.table(gene.manifest, file = file.path(savePath, paste(sampleName,"geneManifest.txt",sep="")),
                quote = F, sep = "\t", row.names = F)

}