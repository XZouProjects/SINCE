# SINCE_demo
SINCE source code


The SINCE algorithm is designed to estimate the number of cell clusters by incorporating existing clustering methods. This repository contains 5 files

SINCE_0.1.0.tar.gz - the R package of SINCE.
run_SC3_SINCE.R - the R script to run a demo analysis.
Pollen_tpm.txt - the data for demo analysis.
Clustering results - the folder containing clustering results with different number of clusters.
cluster_res_list.txt - the directory list of clustering results.
In this demo, SC3 is incorporated with SINCE, although the users are free to explore any other clustering method. 
To run the demo script, the users should install the following packages: knitr, Matrix, scater, SC3, SingleCellExperiment, doSNOW, ClusterR, OGFSC and SINCE. 
The OGFSC package is freely downloadable from github: https://github.com/XZouProjects/OGFSC-R.git
