ARI_matrix <- function(Seurat_Object, resolutions){
  
  # Input : A Seurat object and a list of resolution parameter values
  # Output : The ARI matrix
  # Goal : Compute the ARI matrix for each resolution of the list of resolution
  
  library(Seurat)
  library(mclust)
  
  # Extract cluster identifiers for each resolution
  clusters <- lapply(resolutions, function(r) {
    FetchData(Seurat_Object, vars = paste0("Clusters_", r))[, 1]
  })
  
  # Initialise the ARI matrix
  ari_matrix <- matrix(0, nrow = length(clusters), ncol = length(clusters))
  diag(ari_matrix) <- 1
  
  # Calculate the adjusted Rand index for each pair of clusters
  for (i in 1:(length(clusters) - 1)) {
    for (j in (i + 1):length(clusters)) {
      ari_matrix[i, j] <- adjustedRandIndex(clusters[[i]], clusters[[j]])
      ari_matrix[j, i] <- ari_matrix[i, j]
    }
  }
  
  rm(clusters)
  return(ari_matrix)
}

heatmapari <- function(ari_matrix,resolutions){ 
  
  # Input : An ARI matrix and a list of resolution parameter values
  # Output : The plot of the ARI matrix
  # Goal : Plot the ARI matrix
  
  library(reshape2)
  library(ggplot2)
  
  co <- melt(ari_matrix)
  co$Var1 <- factor(co$Var1, levels = seq_along(resolutions))  # Set levels based on the index of resolutions
  co$Var2 <- factor(co$Var2, levels = seq_along(resolutions))  # Set levels based on the index of resolutions
  
  # Create the heatmap
  gari <- ggplot(co, aes(Var1, Var2)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 2)), size = 2) +
    scale_fill_gradient2(low = "white", 
                         mid = "#e6bffb", 
                         high = "#4ab2f9", 
                         midpoint = 0.5) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12, face = "bold"),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold")) + 
    theme(legend.title = element_text(face = "bold", size = 14)) + 
    scale_x_discrete(name = "", 
                     labels = resolutions) +  # Use resolution values as labels
    scale_y_discrete(name = "", 
                     labels = resolutions) +  # Use resolution values as labels
    labs(fill = "ARI value")
  
  return(gari)
}

NbCluster_gg<-function(Seurat_Object,resolutions){
  df <- data.frame(
    "Resolution" = resolutions,
    "NombreClasses" = sapply(resolutions, function(r) {
      length(unique(FetchData(Seurat_Object, vars = paste0("Clusters_", r))[, 1]))
    })
  )
  
  # Creation of a graph showing the number of clusters
  cluster_plot <- ggplot(df, aes(x = Resolution, y = NombreClasses)) +
    geom_line() +
    geom_point() +
    labs(x = "Resolution", y = "Number of clusters")

  return(cluster_plot)  
}


EffPlot<-function(Seurat_Object,clustname,resolution){
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  #library(tidyverse)
  
  clusters<-as.factor(Seurat_Object[[clustname]][, 1])
  n_clusters <- length(unique(clusters))
  colors <- brewer.pal(min(n_clusters, 12), "Set3")
  if (n_clusters > 12) {
    colors <- colorRampPalette(colors)(n_clusters)
  }
  # Calculating the workforce for each cluster
  Nbcluster <- table(clusters)
  df_Nbcluster<-as.data.frame(Nbcluster)
  names(df_Nbcluster)=c("Cluster","Effectif")
  # Calculating overall proportions
  proportions <- prop.table(Nbcluster)
  df_proportions <- as.data.frame(proportions)
  colnames(df_proportions) <- c("Cluster", "Proportion")
  df_proportions<-merge(df_proportions,df_Nbcluster,by="Cluster",sort=F)
  
  geff<-ggplot(df_proportions, aes(x=Cluster,y=Proportion,fill=Cluster))+
    geom_bar(stat="identity")+
    geom_text(aes(label=Effectif),vjust=1.6,position=position_dodge(0.9),size=3,color="black")+
    theme_minimal()+
    ggtitle(paste("Total proportion of\n individuals by cluster (r =", resolution, ")"))+
    scale_fill_manual(values = colors)
  
  return(geff)
  
}