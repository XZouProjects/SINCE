\name{SINCE}
\alias{SINCE}
\title{Perform SINCE and Plot Analysis Result}
\description{The main function performs SINCE to evaluate clustering results with different parameters.}

\usage{
 SINCE(data, cluster_res_list, CERS_plot_option = 1, r2Cutoff = c(0.5,
	0.7, 0.8), paralSize = 6, miniClusterSize = 5)
	}

\arguments{
\item{data}{The log transformed expression data, cells vs. genes.}
\item{cluster_res_list}{The file contains the list of clustering result directories.}
\item{CERS_plot_option}{The option to plot SINCE results, based on which the users can determine the optimal clustering result.}
\item{r2Cutoff}{The list of candidate parameters to identify variable genes. The values should be within [0, 1]}
\item{paralSize}{The number of repeats SINCE performed.}
\item{miniClusterSize}{The minimum number of cells contained in each cell cluster.}
}

\value{

 \item{CERS}{The calculated CERS (Clustering Error Risk Scores) for different clustering resutls.}
 }
 
 \examples{

# load demo data 

url = "https://raw.githubusercontent.com/XZouProjects/SINCE_test/master/SC3_Pollen_8.xls"
download.file(url, "SC3_Pollen_8.xls")

url = "https://raw.githubusercontent.com/XZouProjects/SINCE_test/master/SC3_Pollen_9.xls"
download.file(url, "SC3_Pollen_9.xls")

url = "https://raw.githubusercontent.com/XZouProjects/SINCE_test/master/SC3_Pollen_10.xls"
download.file(url, "SC3_Pollen_10.xls")

# The clustering result directories are listed in the file named 'cluster_res_list.txt'
url = "https://raw.githubusercontent.com/XZouProjects/SINCE_test/master/cluster_res_list.txt"
download.file(url, "cluster_res_list.txt")

url = "https://raw.githubusercontent.com/XZouProjects/SINCE_test/master/data.rda"
download.file(url, "data.rda")
load("data.rda")

# SINCE analysis

CERS <- SINCE(data = data, cluster_res_list = 'cluster_res_list.txt', paralSize = 4)

}

