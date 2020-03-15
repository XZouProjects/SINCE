# Required packages: knitr, Matrix, scater, SC3, SingleCellExperiment, OGFSC and SINCE.

#### install knitr
# install it from cran

#### install Matrix
# install it from cran

#### install scater
# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("scater", version = "3.8")

#### install SC3
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SC3", version = "3.8")

#### install SingleCellExperiment
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SingleCellExperiment", version = "3.9")

#### install OGFSC 
# The OGFSC package is freely downloadable from github: https://github.com/XZouProjects/OGFSC_R.git

#### install SINCE
#install.packages("SINCE_1.0.0.tar.gz", repos = NULL, type="source")

#####################################################################################

#### SC3 clustering ####

library(knitr)
library(Matrix)
library(scater)
library(SC3)
library(SingleCellExperiment)
library(OGFSC)
library(SINCE)

rm(list = ls())

#### Step 1: Gene filtering by OGFSC ####

data <- read.table("Pollen_tpm.txt", sep="\t", header=TRUE)
rownames(data)<- data[,1]
data <- data[,-1]
data_matrix <- as.matrix(data)

dataBuffer <- log2(data_matrix+1)
OGF <- OGFSC(dataBuffer, plot_option = 1, alpha = c(0.99, 0.999, 0.9999)) 
OGFSC_idx <- OGF$OGFSC_idx
data_matrix <- data_matrix[OGFSC_idx,]

#### Step 2: SC3 based cell clustering ####

N <- c(2:14)

sce<-SingleCellExperiment(assays=list(counts=data_matrix, logcounts=log2(as.matrix(data_matrix)+1)))
rowData(sce)$feature_symbol <- unlist(rownames(data_matrix))

sce<-sc3(sce,ks=N,biology=TRUE,n_cores=6, gene_filter = FALSE)

# save the clustering results and create a results list file
cluster_res_list=list()
i <- 1
for (j in N)
{
  cluster<-eval(parse(text = paste0("sce@colData@listData$", paste0(paste0("sc3_", j),"_clusters")))) 
  write.table(cluster, paste0("SC3_Pollen_", j,".xls") ,sep="\t", row.names = FALSE, col.names = FALSE)
  cluster_res_list[[i]] <- paste0("SC3_Pollen_", j,".xls")
  i <- i+1
}
write.table(cluster_res_list, "cluster_res_list.txt" ,sep="\n", row.names = FALSE, col.names = FALSE)


#### Step 3: SINCE to evaluate all clustering results ####

CERS <- SINCE(data = t(dataBuffer[OGFSC_idx,]), cluster_res_list = 'cluster_res_list.txt', paralSize = 6)
