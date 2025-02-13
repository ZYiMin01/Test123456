preDEheatmap <- function(expr, cell.annotation, genes = NULL, cells = NULL,
                         sel.col.anno = c("Cluster", "Cell.Type"),
                         slot = "scale.data",
                         min.value = -2.5, max.value = 2.5){

    
    if(slot == "scale.data"){
        expr.data <- GetAssayData(object = expr, slot = slot)
    }else{
        expr.data0 <- expr[["RNA"]]@data
        ####z-score
        scaletpm = apply(expr.data0,1,scale)
        expr.data = t(scaletpm)
        colnames(expr.data) = colnames(expr.data0)
    }
    if(!is.null(genes)){
        genes <- intersect(genes, rownames(expr.data))
        expr.data <- expr.data[genes, ]
    }
    if(!is.null(cells)){
        cells <- intersect(cells, colnames(expr.data))
        expr.data <- expr.data[, cells]
    }

    rownames(cell.annotation) <- cell.annotation$barcodes
    cell.annotation <- cell.annotation[colnames(expr.data), ]
    cell.cluster <- cell.annotation[order(cell.annotation$Cluster), sel.col.anno, drop = FALSE]
    expr.data <- expr.data[, rownames(cell.cluster)]

    ## limitData
    expr.data <- limitData(expr.data, min = min.value, max = max.value)

    ## gaps_col
    num.cluster <- table(cell.cluster$Cluster)
    num.cluster <- num.cluster[as.character(0 : (length(num.cluster)-1))]
    gaps_col <- cumsum(num.cluster)

    return(list(expr.data = expr.data,
                cell.cluster = cell.cluster,
                gaps_col = gaps_col))
}

#' plot_markergene_heatmap
#'
#' Plot marker gene heatmap.
#'
#' @param expr A Seurat object.
#' @param savePaths A vecotr of paths containing the results files.
#'
#' @export
#'

plot_markergene_heatmap <- function(expr, genes,name,slot = "data"){
    cell.annotation <- data.frame(barcodes = colnames(x = expr), stringsAsFactors = F)
    cell.annotation <- cbind(cell.annotation, expr@reductions$tsne@cell.embeddings)
    cell.annotation$Cluster <- factor(Idents(object = expr))
    clusters <- unique(cell.annotation$Cluster)
    clusters <- sort(clusters)

    def.colors <- getDefaultColors(n = length(clusters))
    cluster.colors <- c()
    for(i in 1:length(clusters)){
        cluster.colors[as.character(clusters[i])] = def.colors[i]
    }  
   
    de.pre <- preDEheatmap(expr = expr,
                               cell.annotation = cell.annotation,
                               genes = genes,
                               sel.col.anno = c("Cluster"),
                               slot = slot,
                               min.value = -2.5, max.value = 2.5)
    ann_colors = list(Cluster = cluster.colors)
    p.de.heatmap <-
        pheatmap(de.pre$expr.data,
                 color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
                 annotation_col = de.pre$cell.cluster,
                 annotation_colors = ann_colors,
                 fontsize = 7,
                 gaps_col = de.pre$gaps_col,
                 cluster_rows = F, cluster_cols = F,
                 show_colnames = F,
                 # legend = F,
                 silent = T)
    DEplot.height <- 0.5 + 0.1 * length(genes)
    ggsave(filename = paste(name,"_marker_heatmap.pdf",sep=""),p.de.heatmap, width = 8, height = DEplot.height, dpi = 500)
	ggsave(filename = paste(name,"_marker_heatmap.png",sep=""),p.de.heatmap, width = 8, height = DEplot.height, dpi = 500)
}



#' plot_pca_cor
#'
#' Perform sample PCA and cluster cor analysis.
#'
#' @param expr A Seurat object.
#' @param savePaths A vecotr of paths containing the results files.
#'
#' @export
#'
plot_pca_cor <- function(expr,savePath){
    nsample <- length(unique(expr[["stim"]][,1]))
    data = t(as_matrix(expr[["RNA"]]@data))
    data = data.frame(data,cluster=paste("C",Idents(expr),sep=""),sample=expr[["stim"]])
    cluster_gene_mean <- aggregate(x = data[,1:(ncol(data)-2)],
                            by = list(data$cluster),
                            FUN = mean)
    rownames(cluster_gene_mean) <- cluster_gene_mean$Group.1
    cluster_gene_mean <- t(cluster_gene_mean[,-1])
    t_cor=cor(cluster_gene_mean[,colnames(cluster_gene_mean)],method="pearson")
    
    p.de.heatmap0 <- pheatmap(t_cor,kmeans_k=NA,color=colorRampPalette(c("blue", "white", "red"))(50), scale="none",
        breaks=NA, border_color=NA,
        legend=TRUE, drop_levels = FALSE,
        show_rownames=F,show_colnames=T,
        main=NA, fontsize=10,
        cluster_rows = T, cluster_cols = T
    )
    ggsave(filename = file.path(savePath, "cluster_cor_heatmap.pdf"),p.de.heatmap0, width = 8, height =8, dpi = 500)
	ggsave(filename = file.path(savePath, "cluster_cor_heatmap.png"),p.de.heatmap0, width = 8, height =8, dpi = 500)
    if(nsample>2 & nsample<30){
        sample_gene_mean <- aggregate(x = data[,1:(ncol(data)-2)],
                                by = list(data$stim),
                                FUN = mean)
        rownames(sample_gene_mean) <- sample_gene_mean$Group.1
        sample_gene_mean <- sample_gene_mean[,-1]
        bb<-rownames(sample_gene_mean) 
        data.pca <- prcomp(sample_gene_mean,scale=F)
        a <- summary(data.pca)
        tmp <- a$importance
        pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
        pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
        pca_reuslt<-as.data.frame(a$x)  
        pca_reuslt$sample <- bb  
        n<- int(length(bb))/26+2 
        vals <- rep(seq(0,25),n)
        p1<-ggplot(pca_reuslt)+geom_point(aes(x=pca_reuslt[,1],y=pca_reuslt[,2],color=sample,shape = sample),size=4)+scale_shape_manual(values = vals)
        p1<-p1+ theme_bw() +theme(legend.title =element_blank())+labs(x=paste("PC1(",pro1,"%)",sep=""),y=paste("PC2(",pro2,"%)",sep=""))+ theme(axis.title=element_text(size=18),axis.text.x=element_text(face="bold",size=15),axis.text.y=element_text(face="bold",size=15),legend.text=element_text(size=15))
        ggsave(file.path(savePath,paste("sample_pca.pdf",sep="")),  width=8, height=6)
		ggsave(file.path(savePath,paste("sample_pca.png",sep="")),  width=8, height=6)
    } 
}

plot_cluster_by_sample <- function(expr,savePath,sampleNames,name){
    nsample <- length(sampleNames)
    ncluster <- length(Idents(expr)[!duplicated(Idents(expr))])
    if(ncluster > 16){
        m=0.5
    }else{
        m=0
    }
    if(nsample == 2){
        w = 7+m;h=4;n=2
    }else if(nsample == 3){
        w = 10+m;h=4;n=3
    }else if(nsample == 4){
        w = 6+m;h=6;n=2
    }else if(nsample <= 6){
        w = 10+m;h=7;n=3
    }else if(nsample > 6 & nsample <= 9){
        w = 10+m;h=10;n=3
    }else if(nsample >= 10 & nsample <= 12){
        w = 10+m;h=7;n=4
    }else if(nsample >= 13 & nsample <= 16){
        w = 10+m;h=10;n=4
    }else if(nsample > 16 ){
        n = int(nsample-16)/4
        w = 10+m;h=10+n*3;n=4
    }
    pdf(file.path(savePath,paste(name,"_split_sample_by_tsne.pdf",sep="")),width=w,height=h,onefile=F)
    print(DimPlot(expr, reduction = "tsne", split.by = "stim",ncol=n))
    dev.off()
	png(file.path(savePath,paste(name,"_split_sample_by_tsne.png",sep="")),width=w*100,height=h*100)
    print(DimPlot(expr, reduction = "tsne", split.by = "stim",ncol=n))
    dev.off()
    pdf(file.path(savePath,paste(name,"_split_sample_by_umap.pdf",sep="")),width=w,height=h,onefile=F)
    print(DimPlot(expr, reduction = "umap", split.by = "stim",ncol=n))
    dev.off()
	
	png(file.path(savePath,paste(name,"_split_sample_by_umap.png",sep="")),width=w*100,height=h*100)
    print(DimPlot(expr, reduction = "umap", split.by = "stim",ncol=n))
    dev.off()
}

#' run_cluster
#'
#' Perform cell cluster.
#'
#' @param expr A Seurat object.
#' @param sampleNames sample name list.
#' @param savePaths A vecotr of paths containing the results files.
#' @param sampleNames A vector of labels for all samples.
#' @param resolution A float number used in function 'FindClusters' in Seurat. The default is 0.8.
#' @param clusterStashName A character string used as the name of cluster identies. The default is "default".
#' @inheritParams runScCombination
#'
#' @return A seurat objects .
#' @export
#'

run_cluster <- function(expr,sampleNames, savePath,
                             resolution = 0.8,
                             clusterStashName = "default"
                             ){
    expr <- FindClusters(expr, resolution = resolution, verbose = F) 
    ncluster <- length(unique(Idents(expr)))
    Idents(expr)=factor(Idents(expr),levels=0:(ncluster-1))
    expr[[clusterStashName]] <- Idents(object = expr)

    resolution <- resolution
    if(ncluster < 15 ){
        resolution <- resolution + 0.2
        expr <- FindClusters(expr, resolution = resolution, verbose = F)       
        ncluster <- length(Idents(expr)[!duplicated(Idents(expr))])
        Idents(expr)=factor(Idents(expr),levels=0:(ncluster-1))
        expr[[clusterStashName]] <- Idents(object = expr)
    }else if(ncluster > 30){
        resolution <- resolution - 0.2
        expr <- FindClusters(expr, resolution = resolution, verbose = F)
        ncluster <- length(Idents(expr)[!duplicated(Idents(expr))])
        Idents(expr)=factor(Idents(expr),levels=0:(ncluster-1))
        expr[[clusterStashName]] <- Idents(object = expr)
    }
    write.table(resolution, file = file.path(savePath, "resolution.txt"),
                quote = F, sep = "\t", row.names = F)
    pdf(file.path(savePath,"all_cluster_by_umap.pdf"),width=6.2,height=4.5,onefile=F)
    print(pumap2 <- DimPlot(expr, reduction = "umap", label = TRUE))
    dev.off()
    pdf(file.path(savePath,"all_cluster_by_tsne.pdf"),width=6.2,height=4.5,onefile=F)
    print(ptsne2 <- DimPlot(expr, reduction = "tsne", label = TRUE))
    dev.off()
    png(file.path(savePath,"all_cluster_by_umap.png"),width=620,height=450)
    print(pumap2 <- DimPlot(expr, reduction = "umap", label = TRUE))
    dev.off()
    png(file.path(savePath,"all_cluster_by_tsne.png"),width=620,height=450)
    print(ptsne2 <- DimPlot(expr, reduction = "tsne", label = TRUE))
    dev.off()   	
    meta.data = data.frame(expr[[c("nFeature_RNA","nCount_RNA","stim")]],Idents(expr))
    colnames(meta.data) <- c("nGene","nUMI","Samples","clusters")
    
    
    write.table(meta.data,file=file.path(savePath,"all_cells_features.csv"),sep=',',quote=F,col.names=NA)
    
    sample_count = ddply(meta.data,.(Samples),summarize,Number=length(Samples))
    median_gene = ddply(meta.data,.(Samples),summarize,Number=median(nGene))
    mean_UMI = ddply(meta.data,.(Samples),summarize,Number=mean(nUMI))
    meta.data2 = cbind(sample_count,median_gene[,'Number'],mean_UMI[,'Number'])
    colnames(meta.data2) <- c("Samples",'cell.count',"gene.median","UMI.mean")
    write.table(meta.data2,file=file.path(savePath,"sample_cell_gene_filter.csv"),sep=',',quote=F,row.names=F)

    data_tsne = cbind(barcodes=rownames(expr[['stim']]),expr@reductions$tsne@cell.embeddings)
    data_tsne[,'barcodes'] = gsub("_","-",data_tsne[,'barcodes'])
    write.table(data_tsne,"data_tsne.csv",sep=",",row.names=F)
    data_umap = cbind(expr@reductions$umap@cell.embeddings)
    rownames(data_umap) = gsub("_","-",rownames(data_umap))
    write.table(data_umap,"data_umap.csv",sep=",",col.names=NA)
    data_cluster = data.frame(cluster=Idents(expr))
    rownames(data_cluster) = gsub("_","-",rownames(data_cluster))
    write.table(data_cluster,"data_cluster.csv",sep=",")
    
    saveRDS(expr, file = file.path(savePath, "combin.data.RDS"))
    plot_pca_cor(expr,savePath)
    ###plot by samples split
    if(length(sampleNames)>1){
        plot_cluster_by_sample(expr,savePath,sampleNames,"all_cluster")        
    }
    return(expr)
}


#' runScCombination
#'
#' Perform multi-samples analyses.
#'
#' @param savePaths A vecotr of paths containing the results files.
#' @param sampleNames A vector of labels for all samples.
#' @param combName A label for the combined samples.
#' @param comb.method The method to combine samples. The default is "NormalMNN". "NormalMNN", "SeuratMNN", "Raw" and "Regression" are optional.
#'
#' @return A seurat objects .
#' @export
#'
#' @import plyr
#' @importFrom pheatmap pheatmap
#' @importFrom rlang int
#'
 
runScCombination <- function(sampleNames, savePath,
                             comb.method = "NormalMNN",
                             vars.to.regress = c("nUMI", "mito.percent", "ribo.percent"),
                             pc.use = 30,
                             features = 2000
                             ){

    message("[", Sys.time(), "] -----: sample data combination")
    expr.list <- list()
    sample.ident <- c()
    for(i in 1:length(sampleNames)){
        sampleName <- sampleNames[i]
        print(sampleName)
        expr.list[[sampleName]] <- readRDS(paste0(savePath,"/",sampleName,".expr.RDS"))
        sample.ident <- c(sample.ident, rep(sampleName, dim(expr.list[[sampleName]])[2]))
    }
    sample.ident <- as.factor(sample.ident)

    bool.plotHVG = T
    if(comb.method == "SeuratMNN"){
        message("[", Sys.time(), "] -----: combine data by Seurat MNN")
        suppressWarnings( expr.anchors <- FindIntegrationAnchors(object.list = expr.list,
                                                                 dims = 1:pc.use,anchor.features = features))
        expr <- IntegrateData(anchorset = expr.anchors,
                              dims = 1:pc.use, verbose = F)
        expr <- ScaleData(expr, verbose = FALSE)
        DefaultAssay(expr) <- "integrated"
        expr[["sample.ident"]] <- sample.ident
        bool.plotHVG = F

        saveRDS(expr.anchors@anchors, file = file.path(savePath, "anchors.RDS"))

    }else if(comb.method == "None"){
        message("[", Sys.time(), "] -----: get one sample source")
        suppressWarnings( expr <- expr.list[[1]] )
        expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = features, verbose = F)
        expr <- ScaleData(object = expr, vars.to.regress = vars.to.regress, verbose = F)
        expr[["sample.ident"]] <- sample.ident
    }else if(comb.method == "harmony"){
	    message("[", Sys.time(), "] -----: combine data by harmony")
		library(harmony)
		suppressWarnings( expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]) )
		expr <- NormalizeData(expr) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
        ##整合
        expr <- RunHarmony(expr, group.by.vars = "stim")

	}
	if(comb.method != "harmony"){
    message("[", Sys.time(), "] -----: PCA")
    expr <- RunPCA(expr, verbose = F)
	}
    ####PCA plot
    pdf(file.path(savePath,paste("PCA_DimHeatmap_1-6.pdf",sep="")),width=8,height=8,onefile=F)
    print(DimHeatmap(expr, dims = 1:6, cells = 500, balanced = TRUE))
    dev.off()
    pdf(file.path(savePath,paste("PCA_DimHeatmap_7-12.pdf",sep="")),width=8,height=8,onefile=F)
    print(DimHeatmap(expr, dims = 7:12, cells = 500, balanced = TRUE))
    dev.off()
    pdf(file.path(savePath,paste("PCA_Dim_distribution.pdf",sep="")),width=8,height=7,onefile=F)
    print(ElbowPlot(expr))
    dev.off()
	png(file.path(savePath,paste("PCA_DimHeatmap_1-6.png",sep="")),width=800,height=800)
    print(DimHeatmap(expr, dims = 1:6, cells = 500, balanced = TRUE))
    dev.off()
    png(file.path(savePath,paste("PCA_DimHeatmap_7-12.png",sep="")),width=800,height=800)
    print(DimHeatmap(expr, dims = 7:12, cells = 500, balanced = TRUE))
    dev.off()
    png(file.path(savePath,paste("PCA_Dim_distribution.png",sep="")),width=800,height=700)
    print(ElbowPlot(expr))
    dev.off()
    message("[", Sys.time(), "] -----: clustering")
    if (comb.method != "harmony") {
        expr <- FindNeighbors(expr, dims = 1:pc.use, verbose = F)
        message("[", Sys.time(), "] -----: tSNE")
        expr <- RunTSNE(object = expr, dims = 1:pc.use)
        expr <- RunUMAP(expr, dims = 1:pc.use, verbose = F)
        message("[", Sys.time(), "] -----: UMAP")
    }else{
    expr <- FindNeighbors(expr, dims = 1:pc.use, verbose = F,reduction = "harmony")
    message("[", Sys.time(), "] -----: tSNE")
    expr <- RunTSNE(object = expr, dims = 1:pc.use,reduction = "harmony")
    message("[", Sys.time(), "] -----: UMAP")
    expr <- RunUMAP(expr, dims = 1:pc.use, verbose = F,reduction = "harmony")
    }
    #saveRDS(expr, file = file.path(savePath, "combin.data.RDS"))
    expr <- run_cluster(expr,sampleNames, savePath,
                resolution = 0.5,
                clusterStashName = "default")
           
    return(expr)
}

#' runDiffExpr
#'
#' Perform cluster diff gene analyses.
#'
#' @param expr A Seurat object.
#' @param savePaths A vecotr of paths containing the results files.
#' @param name A label for diff analysis.
#' @param n.markers A number for top marker gene for singler cluater.
#' @param clusterStashName A character string used as the name of cluster identies. The default is "default".
#'
#' @return A seurat objects .
#' @export
#'
#' @importFrom dplyr "%>%" top_n group_by
#'

runDiffExpr <- function(expr,savePath,name,n.markers=3,clusterStashName = "default"){
    cell.annotation <- data.frame(barcodes = colnames(x = expr), stringsAsFactors = F)
    cell.annotation <- cbind(cell.annotation, expr@reductions$tsne@cell.embeddings)
    cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]])
    clusters <- unique(cell.annotation$Cluster)
    clusters <- sort(clusters)

    def.colors <- getDefaultColors(n = length(clusters))
    cluster.colors <- c()
    for(i in 1:length(clusters)){
        cluster.colors[as.character(clusters[i])] = def.colors[i]
    }
    diff.expr.genes <- FindAllMarkers(expr, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
    diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")]
    diff.expr.genes[,"cluster"] = factor(diff.expr.genes[,"cluster"],levels=0:(length(clusters)-1))
    diff.expr.genes <- diff.expr.genes[order(diff.expr.genes[,"cluster"]),]
    save(diff.expr.genes,file=file.path(savePath, paste(name,".markers.Rdata",sep="")))
    write.table(diff.expr.genes,paste(name,"_ConservedMarkers.csv",sep=""),sep=',',quote=F,row.names=F) 
    top.genes <- diff.expr.genes %>% group_by(cluster) %>% top_n(n = n.markers, wt = avg_logFC)
    #top.genes <- top.genes[order(top.genes$cluster, top.genes$avg_logFC, decreasing = c(F, T)), ]
    
    de.pre <- preDEheatmap(expr = expr,
                               cell.annotation = cell.annotation,
                               genes = top.genes$gene,
                               sel.col.anno = c("Cluster"),
                               slot = "scale.data",
                               min.value = -2.5, max.value = 2.5)
    ann_colors = list(Cluster = cluster.colors)
    p.de.heatmap <-
        pheatmap(de.pre$expr.data,
                 color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
                 annotation_col = de.pre$cell.cluster,
                 annotation_colors = ann_colors,
                 fontsize = 7,
                 gaps_col = de.pre$gaps_col,
                 cluster_rows = F, cluster_cols = F,
                 show_colnames = F,
                 # legend = F,
                 silent = T)
    DEplot.height <- 0.5 + 0.1 * n.markers * length(unique(cell.annotation$Cluster))
    ggsave(filename = file.path(savePath, paste(name,"_top_marker_heatmap.pdf",sep="")),p.de.heatmap, width = 8, height = DEplot.height, dpi = 500)
    ggsave(filename = file.path(savePath, paste(name,"_top_marker_heatmap.png",sep="")),p.de.heatmap, width = 8, height = DEplot.height, dpi = 500)
    ncluster  = length(clusters)    
    nn = dim(top.genes)[1]
    if(ncluster<10){
        nn2=8
    }else{
        nn2=12
    }
    if(ncluster>17){nn2=16}
    pdf(file.path(savePath, paste(name,"_top_marker_violin_nodot.pdf",sep="")),width=nn2+1,height=nn/1.7,onefile=F) 
    print(VlnPlot(expr,features = top.genes$gene, group.by = "seurat_clusters",pt.size = 0))
    dev.off()
    png(file.path(savePath, paste(name,"_top_marker_violin_nodot.png",sep="")),width=(nn2+1)*100,height=nn/1.7*100) 
    print(VlnPlot(expr,features = top.genes$gene, group.by = "seurat_clusters",pt.size = 0))
    dev.off()
    pdf(file.path(savePath, paste(name,"_top_marker_tsne.pdf",sep="")),width=12,height=nn/1.7,onefile=F)
    print(FeaturePlot(expr, features = top.genes$gene, cols = c("grey", "red"),reduction = "tsne"))
    dev.off()
    png(file.path(savePath, paste(name,"_top_marker_tsne.png",sep="")),width=12*100,height=nn/1.7*100)
    print(FeaturePlot(expr, features = top.genes$gene, cols = c("grey", "red"),reduction = "tsne"))
    dev.off()
    
    #pdf(file.path(savePath, paste(name,"_top_marker_heatmap1.pdf",sep="")),width=nn2,height=nn/4.5,onefile=F)
    #print(DoHeatmap(expr, features = top.genes$gene,size = 2) + NoLegend())
    #dev.off()
    return(diff.expr.genes)   
}

