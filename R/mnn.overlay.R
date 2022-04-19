library(batchelor)
library(RANN)
library(SingleCellExperiment)
library(rlang)
library(Seurat)


#matrix1 is the new data, rows are cells and columns are markers
#matrix2 is the reference data, rows are cells and columns are markers
#k.param number of neighbors used for the overlay
#ref.position 2-column header of dim reduction, UMAP, TSNE
#ref.assignments cluster or annotation vector 
#markers slection of markers to integrate. If left NULL, the common markers between matrix1 and matrix2 will be used

#This will take awhile depending on the number of cells - testd on matrix1 = 4e6, matrix2 = 1e6

MNN.overlay <- function(matrix1, 
                        matrix2, 
                        k.param = 50, 
                        ref.positions = NULL,
                        ref.assignments = NULL,
                        markers = NULL) {
    if(!is.null(markers)) {
        markers.use <- markers
    } else {
        markers.use <- intersect(colnames(matrix1), colnames(matrix2))
    }
    #Removing FSC/SSC, if not removed, this ruins the batch correction
    markers.use <- markers.use[!grepl("FSC|SSC", markers.use)]
    matrix.list <- list(matrix1, matrix2)
    #mmatrix.list adds barcode and filtrs for specific markers in order
    for (i in seq_along(matrix.list)) {
        barcode <- seq_len(nrow(matrix.list[[i]]))
        len <- nchar(nrow(matrix.list[[i]]))
        barcode<-sprintf(paste0("%0", len, "d"), barcode)
        barcode <- paste0("matrix.", i, ":", barcode)
        rownames(matrix.list[[i]]) <- barcode
        matrix.list[[i]] <- matrix.list[[i]][,match(markers.use, colnames(matrix.list[[i]]))]
        matrix.list[[i]] <- t(matrix.list[[i]])
    }
    print("Correcting for Batch Effect...")
    sce <- suppressWarnings(fastMNN(matrix.list[[1]], matrix.list[[2]], k = k.param))
    #Pull coredcted dimensions
    corrected.matrix <- reducedDim(sce)
    corrected.matrix1 <- corrected.matrix[!grepl("matrix.2:", rownames(corrected.matrix)),]
    corrected.matrix2 <- corrected.matrix[grepl("matrix.2:", rownames(corrected.matrix)),]
    
    print("Finding cross-batch Neighbors...")
    #Nearest Neighbor calculation based on RANN package
    nn.output <- Seurat:::NNHelper(corrected.matrix2,corrected.matrix1, k = k.param, method = "rann")
    rownames(nn.output@nn.dist) <- rownames(corrected.matrix1)
    rownames(nn.output@nn.idx) <- rownames(corrected.matrix1)
    cellnames <- colnames(matrix.list[[1]])
    nproj <- matrix(data = NA, nrow = length(cellnames), ncol = 4,
                         dimnames = list(cellnames, c("dim_1","dim_2", "pred.cluster", "pred.conf")))
    print("Calculating cell assignments and weighted positions...")
    #Pulled from Projectil to implement UMAP overlay
    for (i in 1:length(cellnames)) {
        row <- exp(-nn.output@nn.dist[cellnames[i],])  #calculate exp(-dist) as weights for nearest neighbors
        weights = row/sum(row)
        top.k <- nn.output@nn.idx[cellnames[i],]
        scores <- sort(table(ref.assignments[top.k]), decreasing = T)/k.param
        pred.type <- names(scores)[1]
        pred.conf <- scores[1]
        nproj[i,] = c(c(weights %*% (ref.positions[nn.output@nn.idx[cellnames[i],],])), pred.type, pred.conf)  #assign UMAP coordinates of (weighted) neighbors
    }
    return(nproj)
}
    
        
    