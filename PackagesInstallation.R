if(!require("pak")){
  install.packages("pak")
  require("pak")
}
# System Librairies required : libhdf5-dev, libgsl27
# Maybe libglpk-dev

if(!require("BPCells")){pak::pkg_install("bnprks/BPCells/r",dependencies = T)}
if(!require("presto")){pak::pkg_install("immunogenomics/presto",dependencies = T)}
if(!require("Seurat")){pak::pkg_install("satijalab/seurat",dependencies = T)}
if(!require("SeuratData")){pak::pkg_install("satijalab/seurat-data")}
if(!require("harmony")){pak::pkg_install("immunogenomics/harmony",dependencies = T)}
if(!require("tidyverse")){pak::pkg_install("tidyverse",dependencies = T)}
if(!require("RColorBrewer")){pak::pkg_install("RColorBrewer",dependencies = T)}
if(!require("reshape2")){pak::pkg_install("reshape2",dependencies = T)}
if(!require("FactoMineR")){pak::pkg_install("FactoMineR",dependencies = T)}
if(!require("corrplot")){pak::pkg_install("corrplot",dependencies = T)}
if(!require("mclust")){pak::pkg_install("mclust",dependencies = T)}
if(!require("clustree")){pak::pkg_install("clustree",dependencies = T)}
if(!require("reshape2")){pak::pkg_install("reshape2",dependencies = T)}
if(!require("tidyverse")){pak::pkg_install("tidyverse",dependencies = T)}
if(!require("DESeq2")){pak::pkg_install("DESeq2",dependencies = T)}
if(!require("metap")){pak::pkg_install("metap",dependencies = T)}
if(!require("nebula")){pak::pkg_install("lhe17/nebula",dependencies = T)}
if(!require("randomForest")){pak::pkg_install("randomForest",dependencies = T)}
if(!require("randomForestVIP")){pak::pkg_install("randomForestVIP",dependencies = T)}
if(!require("ggvenn")){pak::pkg_install("ggvenn",dependencies = T)}
if(!require("pheatmap")){pak::pkg_install("pheatmap",dependencies = T)}
if(!require("clusterProfiler")){pak::pkg_install("clusterProfiler",dependencies = T)}
if(!require("enrichplot")){pak::pkg_install("enrichplot",dependencies = T)}
if(!require("org.Hs.eg.db")){pak::pkg_install("org.Hs.eg.db",dependencies = T)}
if(!require("ggupset")){pak::pkg_install("ggupset",dependencies = T)}
if(!require("BayesSpace")){pak::pkg_install("BayesSpace",dependencies = T)}
# asked for SpatialFeaturePlot(), skip if not spatial
if(!require("glmGamPoi")){pak::pkg_install("glmGamPoi",dependencies = T)} 
