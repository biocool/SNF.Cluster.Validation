# SNF_Louvain.R
#
# 1) Perform similarity fusion network analysis with 6 ROI variables (temporal ROI 
# activation from three language paradigms) and 14 clinical variables; 
# 2) Run Louvain clustering analysis with the similarity network (strongest 15% 
# connecting partners of each subject) 
# 
# INPUT
#

ROI_var <- 
  colnames(yaqoing.input.data)[c(5:8)]
clinic_var <- 
  colnames(yaqoing.input.data)[c(9:14)]

temp <- 
  SNF_Louvain(dat = yaqoing.input.data, 
              ROI_var = ROI_var, 
              clinic_var = clinic_var)

SNF_Louvain <- function(dat, ROI_var, clinic_var) {
	
	# load libraryies
	library(SNFtool)
	library(ggplot2)
	library(factoextra)
	library(effsize)
	library(data.table)
	library(igraph)
	
	# normalization variables for each data type
	ROI_dat <- standardNormalization(dat[,ROI_var])
	clinc_data <- standardNormalization(dat[,clinic_var])
	
	# pair-wise distance
	DistROIs = (dist2(as.matrix(ROI_dat),as.matrix(ROI_dat)))^(1/2)
	DistClinic = (dist2(as.matrix(clinc_data),as.matrix(clinc_data)))^(1/2)
	
	# add row and column names
	rownames(DistROIs) <- dat$subj
	rownames(DistClinic) <- dat$subj
	
	colnames(DistROIs) <- dat$subj
	colnames(DistClinic) <- dat$subj
	
	# similarity graphs
	WROIs <- affinityMatrix(DistROIs, K = 20, sigma = 0.5)
	WClinic <- affinityMatrix(DistClinic, K = 20, sigma = 0.5)
	
	# fuse/merge two matricies into one similarity network
	W = SNF(list(WROIs, WClinic), 20, 20)

	# visualization
	dat$index[dat$Dx == "ASD"] <- 1
	dat$index[dat$Dx == "ASDFeat"] <- 2
	dat$index[dat$Dx == "DD/Other"] <- 3
	dat$index[dat$Dx == "LD"] <- 4
	dat$index[dat$Dx == "TD"] <- 5
	
	fviz_cluster(object = list(data = W, cluster = dat$index), geom = "point",
			  ellipse = F, show.clust.cent = F, shape = 19,
			  pointsize = 3, ellipse.type = "norm", ggtheme = theme_void())
	
	# reorganize similarity network
	ind <- which(lower.tri(W,diag=F) , arr.ind = TRUE)
	
	snf_matrix <- as.data.frame(ind)
	for (x in 1:dim(ind)[1]) {
		snf_matrix[x, 3:5] <- cbind(colnames(W)[ind[x, "col"]],
					    colnames(W)[ind[x, "row"]],
					    round(W[colnames(W)[ind[x, "col"]],
					    	colnames(W)[ind[x, "row"]]],5))
	}
	colnames(snf_matrix)[3:5] <- c("Source", "Target", "Weight")
	snf_matrix$name <- paste(as.character(snf_matrix$Source), 
				 as.character(snf_matrix$Target), sep = ".")
	
	# threshold the similarity network
	snf_matrix_thr <- as.data.frame(matrix(0, 0, 6))
	for (i in 1:dim(W)[1]) {
		tmp_weight <- snf_matrix[snf_matrix$col == i, ]
		top <- tmp_weight[as.numeric(tmp_weight$Weight) > quantile(as.numeric(tmp_weight$Weight), prob = 1 - 15/100),]
		snf_matrix_thr <- rbind(top, snf_matrix_thr)
	}
	
	dim(snf_matrix_thr)
	
	# thresholded matrix

	matrix_thr <- as.data.frame(matrix(0, nrow(dat), nrow(dat)))
	rownames(matrix_thr) <- dat$subj
	colnames(matrix_thr) <- dat$subj
	
	for (i in seq_along(snf_matrix_thr$Source)) {
		matrix_thr[snf_matrix_thr$Source[i], snf_matrix_thr$Target[i]] <- 
			as.numeric(snf_matrix_thr$Weight[i])
	}
	
	
	# Louvain clustering analysis
	gg <- as.matrix(matrix_thr)
	G1 <- graph_from_adjacency_matrix(gg, mode = "undirected", weighted = T, 
					  diag = F)

	clusterlouvain <- igraph::cluster_louvain(G1)
	nclust <- max(clusterlouvain$membership)
	
	print(paste0("There are ",  nclust, " clusters")) 
	
	# clustering results
	W_Clustering_thr <- as.data.frame(matrix(0, dim(W)[1], 3))
	colnames(W_Clustering_thr) <- c("subj", "Clustering", "group")
	
	W_Clustering_thr$Clustering <- clusterlouvain$membership
	W_Clustering_thr$subj <- dat$subj
	W_Clustering_thr$group <- dat$group
	
	return(list(snf_matrix_thr, W_Clustering_thr))
}

