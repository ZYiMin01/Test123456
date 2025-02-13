
get_singleR_result <- function(data,anno,name){
    data[["SingleR.labels"]] <- as.character(anno$labels)
    cluster_type <- tapply(Idents(data),data[["SingleR.labels"]],table)
    cluster_type <- as.data.frame(bind_cols(cluster_type))
    rownames(cluster_type)=names(cluster_type[[1]])
    write.table(cluster_type,paste("celltype_anno_cluster_",name,".xls",sep=""),sep="\t")
    cluster_anno <- apply(cluster_type, 1, function(t) colnames(cluster_type)[which.max(t)])
    cluster_anno <- as.data.frame(cluster_anno)
    cluster.labels_main <- cluster_anno$cluster_anno[match(Idents(data), rownames(cluster_anno))]
    
    cluster.labels <- as.factor(paste("C",Idents(data),":",cluster_anno$cluster_anno[match(Idents(data), rownames(cluster_anno))],sep=""))
    cluster.labels <- factor(cluster.labels,levels=paste("C",rownames(cluster_anno),":",cluster_anno$cluster_anno,sep=""))
    data[["clusters"]] = cluster.labels
    anno_result <- data.frame(data[["SingleR.labels"]],cluster.labels_main,cluster.labels)
    anno_result1 <- data.frame(rownames(data[['stim']]),Idents(data),data[['stim']],cluster.labels_main,cluster.labels)
    colnames(anno_result) <- c('SingleR.labels','cluster.combin.labels','cluster.labels')
    colnames(anno_result1) <- c('barcodes','cluster','sample','cluster.combin.labels','cluster.labels')
    write.table(anno_result,paste("SingleR_anno_",name,"_result.xls",sep=""),sep="\t")
    write.table(anno_result1,paste("cell_sample_cluster_annotation_",name,".xls",sep=""),sep="\t",row.names=F)
	anno_result2 <- anno_result1
	anno_result2[,'barcodes'] <- gsub('_','-',anno_result2[,'barcodes'])
	write.table(anno_result2,paste("loupe_annotation_",name,".csv",sep=""),sep=",",row.names=F)
    ####plot
    data2 = data
    Idents(data2) <- anno$labels
    w <- 5+get_lable_len(Idents(data2))
    pdf(paste("celltype_anno_cell_",name,"_by_tsne.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data2, reduction = "tsne" ))
    dev.off()
    pdf(paste("celltype_anno_cell_",name,"_by_umap.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data2, reduction = "umap" ))
    dev.off()
	png(paste("celltype_anno_cell_",name,"_by_tsne.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data2, reduction = "tsne" ))
    dev.off()
    png(paste("celltype_anno_cell_",name,"_by_umap.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data2, reduction = "umap" ))
    dev.off()
    Idents(data2) <- cluster.labels
    lables=levels(Idents(data2))
    w <- 5+get_lable_len(Idents(data2))
    pdf(paste("celltype_anno_cluster_",name,"_by_tsne.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data, reduction = "tsne", label =TRUE)+scale_colour_discrete(labels=c(lables)))
    dev.off()
    pdf(paste("celltype_anno_cluster_",name,"_by_umap.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data, reduction = "umap", label =TRUE)+scale_colour_discrete(labels=c(lables)))
    dev.off()
	png(paste("celltype_anno_cluster_",name,"_by_tsne.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data, reduction = "tsne", label =TRUE)+scale_colour_discrete(labels=c(lables)))
    dev.off()
    png(paste("celltype_anno_cluster_",name,"_by_umap.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data, reduction = "umap", label =TRUE)+scale_colour_discrete(labels=c(lables)))
    dev.off()
    Idents(data2) <- cluster.labels_main
    w <- 5+get_lable_len(Idents(data2))
    pdf(paste("celltype_anno_cluster_combin_",name,"_by_tsne.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data2, reduction = "tsne", label =TRUE))
    dev.off()
    pdf(paste("celltype_anno_cluster_combin_",name,"_by_umap.pdf",sep=""),width=w,height=5,onefile=F)
    print(ptsne <- DimPlot(data2, reduction = "umap", label =TRUE))
    dev.off()
	png(paste("celltype_anno_cluster_combin_",name,"_by_tsne.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data2, reduction = "tsne", label =TRUE))
    dev.off()
    png(paste("celltype_anno_cluster_combin_",name,"_by_umap.png",sep=""),width=w*100,height=500)
    print(ptsne <- DimPlot(data2, reduction = "umap", label =TRUE))
    dev.off()
    return(data)
}

#' get_lable_len
#'
#' get the length of lable.
#'
#' @export
#'
get_lable_len <- function(labels){
    len = length(unique(labels))
    charmax = max(nchar(as.character(unique(labels))))
    nn = floor((len-1)/20)
    nn2 = (nn+1)*charmax*0.085
    return (nn2)
}

anno_one_database <- function(combin.data,type,species="human"){
    data = combin.data[["RNA"]]@data
    if(type  %in% c('RNA_type','Immune','Blueprint','Noversht','Monaco')){
        if(species == "human"){
            ref.fine = readRDS(system.file("humandb",paste(type,".RDS",sep=""), package = "singlercell"))
        }else if(species == "mouse"){ 
            ref.fine = readRDS(system.file("mousedb",paste(type,".RDS",sep=""), package = "singlercell"))
        }
        anno.main=SingleR(test =data , ref = ref.fine, labels = ref.fine$label.main)
        anno.fine=SingleR(test =data , ref = ref.fine, labels = ref.fine$label.fine)
    }else{
        if(species == "human"){
            ref.fine = readRDS(system.file("humandb",paste(type,".fine.RDS",sep=""), package = "singlercell"))
            celltype = read.table(system.file("humandb",paste(type,"_celltype.txt",sep=""), package = "singlercell"),header=T,sep="\t")
        }else if(species == "mouse"){ 
            ref.fine = readRDS(system.file("mousedb",paste(type,".fine.RDS",sep=""), package = "singlercell"))
            celltype = read.table(system.file("mousedb",paste(type,"_celltype.txt",sep=""), package = "singlercell"),header=T,sep="\t")
        }
        anno.fine=SingleR(test =data , ref = ref.fine, labels = ref.fine$label)
        anno.main=data.frame(labels = celltype$Celltype2[match(anno.fine$labels, celltype$Celltype1)])
    }
    name1 = paste(type,".main",sep="")
    combin.data = get_singleR_result(combin.data,anno.main,name1)
    combin.data = get_singleR_result(combin.data,anno.fine,paste(type,".fine",sep=""))
    return(combin.data)
}


#' annocelltype
#'
#' Annotation of celltype with database .
#'
#' @param combin.data A seurat object.
#' @param type  A character string  or a vector giving the type of database.
#' @param species A character string  giving the species of database.
#'
#' @export
#'
#' @import SingleR  Seurat 
#' @importFrom dplyr bind_cols
#'

annocelltype <- function(combin.data,type,species="human"){
    if(length(type)==1){
        combin.data = anno_one_database(combin.data,type,species=species)   
    }
    return(combin.data)
}

#' marker_gene_anno
#'
#' Annotation of celltype with marker gene .
#'
#' @param expr A seurat object.
#' @param genefile  marker gene file.
#' @param tissue the tissue of sample.
#'
#' @export
#'
#' @import ggplot2  Seurat 
#'

