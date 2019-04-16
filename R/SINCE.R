SINCE <- function(data, cluster_res_list,  CERS_plot_option = 1, r2Cutoff =c(0.5, 0.7, 0.8),
                  paralSize = 6, miniClusterSize = 5)
{
  # SINCE to evaluate clustering resutls
  # Inputs: data - a log2 transformed scRNA-seq data matrix with rows being genes and columns being cells.
  #         cluster_res_list - a .txt file containing the list of clustering results.
  #         varargin - adjustable parameters, including
  #                    OGFSC_option - the option to perform OGFSC based gene filtering. By default 1.
  #                    CERS_plot_option - the option to show CERS plot. By default 1.
  #                    r2Cutoff - the OPLSDA r2 cutoff values in variable gene selection. By default[0.5, 0.7, 0.8].
  #                    paralSize - number of repeats when calculating CERS. By default 6.
  #                    miniClusterSize - the minimum size of each cell cluster. By default 5.
  # Outputs: CERS - the results matrix, with rows being the CERS with different
  #                 r2Cutoff and columns being different clustering results listed in 'cluster_res_list'.

  ## calculation of CERS

  resultPath <- read.table(file = cluster_res_list, sep = "\n", colClasses = "character")    #read entire file into string

  NumRes <- dim(resultPath)[1]   #number of rows
  flag <- array(0, dim = c(length(r2Cutoff), NumRes, paralSize))
  CERS <- matrix(0,length(r2Cutoff), NumRes)
  flagbuffer <-matrix(0,length(r2Cutoff),paralSize)


  #data = dataBuffer
  #groupIdx = clusterIdx
  #r2Cutoff =r2Cutoff
  #miniClusterSize = miniClusterSize


  cl<-makeCluster(paralSize, type = "SOCK")
  registerDoSNOW(cl)

  for (i in seq_len(NumRes))
  {
    cat ("i = ")
    print(i[1])
    cat ("\n")

    clusterIdx <- read.table(file = resultPath[i,],sep="\t")   #read entire file into string
    nrows <- dim(clusterIdx)[1]   #number of rows


    flagbuffer <-matrix(0,length(r2Cutoff),paralSize)

    flagbuffer <- foreach (paral = seq_len(paralSize),  .combine = 'cbind', .packages =c("ClusterR","e1071","SINCE")) %dopar% {
      Binarization(data = data, groupIdx = clusterIdx, r2Cutoff =r2Cutoff,  miniClusterSize = miniClusterSize)
    }

    for (paral in seq_len(paralSize))
    {
      flag[, i, paral] <- flagbuffer[,paral]
    }


    for (rc in seq_len(length(r2Cutoff)))
    {
      CERS[rc, i] = 1-sum(drop(!flag[rc, i, ]))/paralSize
    }
  }

  stopCluster(cl)

  ## plot CERS
  if (CERS_plot_option == 1)
  {
    mycolors <- brewer.pal(8, "Pastel1")

    for (i in seq_len(dim(CERS)[1]))
    {
      if (i == 1)
      {
        plot(seq_len(dim(CERS)[2]), CERS[i,],type = "o", col =  mycolors[i], lwd =5, cex = 2, cex.lab=1.5,
             pch=16,cex.axis = 1.5, xlab="Clustering results index",ylab="Clusering Error Risk Score",    ylim = c(0,1))
      } else {
        lines(seq_len(dim(CERS)[2]), CERS[i,],type = "o", col =  mycolors[i], lwd =5, cex = 2,pch=16+i-1)
      }
    }
    legend("topleft", legend = paste('Risk Score', r2Cutoff,sep=" "),
           col = mycolors[seq_len(dim(CERS)[1])], pch = 16:(16+dim(CERS)[1]-1),lwd =2,cex = 1)
  }
  return(CERS = CERS)
}

Binarization <- function(data, groupIdx, r2Cutoff = c(0.5, 0.7, 0.8),miniClusterSize = 5)
{
    # Binarization: binarization of gene expression data
    #              'data'- data matrix with cells vs genes
    #              'groupIdx' - the cell categorization results
    #              'miniClusterSize' - the minimum size of each cell cluster.
    # Output: flag - flag indicating clustering error happens. 1 means over-partitioning, otherwise 0.
    
    oplsda_para <- list(NULL)
    oplsda_para$nc  <- 1
    oplsda_para$ncox  <- 1
    oplsda_para$ncoy  <- 0
    oplsda_para$pcutoff  <- 0.05
    oplsda_para$preprocessing  <-'uv'
    
    cutoff <- r2Cutoff
    flag <- rep(-1,length(cutoff))
    temp_ori <- NULL
    ## identify variable genes using different r2Cutoff values
    for (rc in seq_len(length(cutoff)))
    {
        oplsda_para$r2cutoff <- cutoff[rc]
        biomarkersIdx <- NULL
        N <- max(groupIdx)
        # pairwise comparison between clusters
        C <- combn(seq_len(N),2)
        C <- t(C)
        var <- seq_len(dim(data)[2])
        
        for (i in seq_len(N))
        {
            tempData1 <- data[groupIdx==i, ]
            if (dim(tempData1)[1]<miniClusterSize)
            {
                flag[rc]  <- 1
                return (flag = flag)
            }
            
            tempData2  <- data[groupIdx!=i,]
            model  <- OPLSDA(tempData1, tempData2, var, oplsda_para)
            if (is.null(model))
            {
                next
            }
            
            biomarkersIdx  <- c(biomarkersIdx, model$sig_idx)
        }
        
        for (i in seq_len(dim(C)[1]))
        {
            
            tempData1  <- data[groupIdx==C[i,1],]
            tempData2  <- data[groupIdx==C[i,2],]
            model  <- OPLSDA(tempData1, tempData2, var, oplsda_para)
            if (is.null(model))
            {
                next
            }
            biomarkersIdx  <- c(biomarkersIdx, model$sig_idx)
        }
        
        biomarkersIdx  <- sort(unique(biomarkersIdx))
        if (is.null(biomarkersIdx) || length(biomarkersIdx) < 10)
        {
            flag[rc]  <- 1
            next
        }
        
        ## binarization
        tt <-NULL
        mu <- NULL
        X_binary  <- matrix(0,dim(groupIdx)[1], length(biomarkersIdx))
        for (i in seq_len(length(biomarkersIdx)))
        {
            temp  <- data[,biomarkersIdx[i]]
            # # k-means binarization
            IDX  <- kmeans(temp, 2,iter.max = 100, nstart = 20)$cluster#, 'replicates', 20)
            mu[1]  <- mean(temp[IDX == 1])
            mu[2]  <- mean(temp[IDX == 2])
            IX <- which.max(mu)
            thresh  <- min(temp[IDX == IX])
            tt[i] = thresh
            X_binary[temp>=thresh,i]  <- 1
            X_binary[temp<thresh,i] <- 0
        }
        
        # cluster genes and associate each gene cluster to a cell cluster
        
        P  <- vector('list',length(floor(log2(max(groupIdx))):2^max(groupIdx)))
        for (i in floor(log2(max(groupIdx))):2^max(groupIdx)) # categorize markers into numbers of groups which are <= max(groupIdx)
        {
            temp  <- NULL
            
            if (i == 1)
            {
                counts  <- matrix(0, 2,max(groupIdx))
                for (j in seq_len(max(groupIdx)))# for each cell category
                {
                    counts[1,j] <- length(which(X_binary[groupIdx==j,]==1))
                    counts[2,j] <- length(which(X_binary[groupIdx==j,]==0))
                }
                K  <- colSums(counts)
                N  <- rowSums(counts)
                M  <- sum(K)
                P[[i]] <- matrix(0,dim(counts)[2],dim(counts)[1])
                for (ki in seq_len(dim(counts)[2]))
                {
                    P[[i]][ki,seq_len(dim(counts)[1])] <- phyper(counts[,ki], m = N[1], n =N[2], k = K[ki],lower.tail = FALSE ,log.p = FALSE)
                }
                
                temp  <- P[[i]][1,]
                temp[which(P[[i]][1,]<0.05)]  <- 1
                temp[which(P[[i]][1,]>=0.05)] <- 0
                
            }
            else {
                P[[i]]  <- vector('list',i)
                
                IDX <- Cluster_Medoids(t(X_binary), clusters = i, distance_metric = 'hamming', swap_phase = FALSE, seed = 2)$clusters
                #IDX <- kmeans(t(X_binary), ,iter.max = 200, nstart = 10)$cluster#, 'replicates', 20, 'Distance', 'hamming')
                #IDX <- readMat("C:/Users/Jie/Desktop/SC3-SINCE/IDX.mat")$IDX
                geneNgroup <- matrix(0,1,i)
                for (jj in seq_len(i)) # for each marker subgroup
                {
                    
                    counts <- matrix(0, 2,max(groupIdx))
                    for (j in seq_len(max(groupIdx))) # for each cell category
                    {
                        counts[1,j] <- length(which(X_binary[groupIdx==j,IDX==jj]==1))
                        counts[2,j]  <- length(which(X_binary[groupIdx==j,IDX==jj]==0))
                    }
                    K  <- colSums(counts)
                    N  <- rowSums(counts)
                    M  <- sum(K)
                    P[[i]][[jj]] <- matrix(0,dim(counts)[1],dim(counts)[2])
                    
                    for (ki in seq_len(dim(counts)[2]))
                    {
                        for (k2 in seq_len(dim(counts)[1]))
                        {
                            P[[i]][[jj]][k2,ki] <- phyper(counts[k2,ki], m = N[1], n = N[2], k = K[ki],lower.tail = FALSE,log.p = FALSE)
                        }
                        
                    }
                    
                    temp_P  <- P[[i]][[jj]][1,]
                    temp_P[which(P[[i]][[jj]][1,]<0.05)]  <- 1
                    temp_P[which(P[[i]][[jj]][1,]>=0.05)]  <- 0
                    temp  <- rbind(temp, temp_P)
                    geneNgroup[jj] <- length(which(IDX==jj))
                }
                D <- hamming.distance(temp)
                D <- D[upper.tri(D)]/dim(temp)[2]
                if (min(D) == 0 || min(geneNgroup)<5)
                {
                    if (exists("temp_ori") && !is.null(temp_ori))
                    {
                        express_pattern <- temp_ori  # rows are gene clusters and columns are cells clusters
                    } else {
                        flag[rc] <- 1
                        express_pattern <- NULL
                        
                    }
                    break
                }
            }
            temp_ori <- temp
        }
        ## check if duplicate clusters exist
        if (!is.null(express_pattern))
        {
            if (is.matrix(express_pattern))
            {
                ttmp <- hamming.distance(t(express_pattern))
                ttmp <- ttmp[upper.tri(ttmp)]/dim(express_pattern)[1]
            } else {
                ttmp <- hamming.distance(t(t(express_pattern)))
                ttmp <- ttmp[upper.tri(ttmp)]
            }
            D <- min(ttmp)
            if (D>0)
            {
                flag[rc] <- 0
            } else # over-partitioning
            {
                flag[rc] <- 1
            }
        }
        
        
    }
    
    return(flag = flag)
}


mjrO2pls <- function(X, Y, pax, oax, oay,splitSize=300)
{
    #
    ## written by Dr Jie Hao & Dr Xin Zou (2016), Shanghai Jiaotong University, China
    
    Wo<-NULL
    Pyo<-NULL
    To<-NULL
    Us <- NULL
    
    ssx<-sum(X^2)
    ssy<-sum(Y^2)
    
    Tmps <- svd(t(t(Y)%*%X))
    W <- Tmps$u
    S <- matrix(0,length(Tmps$d),length(Tmps$d))
    diag(S) <- Tmps$d
    C <- Tmps$v
    W<-W[,seq_len(pax)]
    C<-C[,seq_len(pax)]
    S<-S[seq_len(pax),seq_len(pax)]
    
    Ts <- NULL
    Ts[[1]]<-X%*%W
    TT <-Ts[[1]]
    
    Exy<-X-TT%*%t(W)
    Xr<-X
    
    R2Xo<-0
    R2Xcorr<-sum((Ts[[1]]%*%t(W))^2)/ssx
    R2X<-1-sum(Exy^2)/ssx
    
    if (oax>=1)
    {for (i in seq_len(oax))
        {
            
            Tmps <-svd(t(Exy)%*%Ts[[i]])
            wo <- Tmps$u
            syo <- matrix(0,length(Tmps$d),length(Tmps$d))
            diag(syo) <- Tmps$d
            wor <- Tmps$v
            if (length(Tmps$d)>1)
            wo<-wo[,1] #/sqrt(wo'*wo)
            
            to<-Xr%*%wo  #/(wo'*wo)
            pyo<-t(Xr)%*%to%*%solve(t(to)%*%to)
            Xr<-Xr-to%*%t(pyo)
            
            Wo<-cbind(Wo,wo)
            Pyo<-cbind(Pyo,pyo)
            To<-cbind(To,to)
            
            Ts[[i+1]]<-Xr%*%W
            Exy<-X-Ts[[i+1]]%*%t(W)-To%*%t(Pyo)
            
            R2Xo<-c(R2Xo,sum((To%*%t(Pyo))^2)/ssx)
            R2Xcorr<-c(R2Xcorr,sum((Ts[[i+1]]%*%t(W))^2)/ssx)
            R2X<-c(R2X,1-sum(Exy^2)/ssx)
            TT<-Ts[[i+1]]
        }
    }
    
    
    Yr<-Y
    Us[[1]]<-Yr%*%C
    U<-Us[[1]]
    Fxy<-Y-(Us[[1]]%*%t(C))
    
    R2Yo<-0
    R2Ycorr<-sum((Us[[1]]%*%t(C))^2)/ssy
    R2Y<-1-sum(Fxy^2)/ssy
    
    Uo<-NULL
    Pxo<-NULL
    Co<-NULL
    
    if (oay >=1)
    {
        for(i in seq_len(oay))
        {
            Tmps <- svd(t(Fxy)%*%Us[[i]])
            co <- Tmps$u
            sxo <- matrix(0,length(Tmps$d),length(Tmps$d))
            diag(sxo) <- Tmps$d
            cor <- Tmps$v
            if (length(Tmps$d)>1)
            co<-co[,1]  #/sqrt(co'*co)
            uo<-Yr%*%co  #/(co'*co)
            pxo<-t(Yr)%*%uo%*%solve(t(uo)%*%uo)
            Yr<-Yr-uo%*%t(pxo)
            Co<-cbind(Co,co)
            Pxo<-cbind(Pxo,pxo)
            Uo<-cbind(Uo,uo)
            
            Us[[i+1]]<-Yr%*%C
            Fxy<-Y-Us[[i+1]]%*%t(C)-Uo%*%t(Pxo)
            
            R2Yo<-c(R2Yo,sum((Uo%*%t(Pxo))^2)/ssy)
            R2Ycorr<-c(R2Ycorr,sum((Us[[i+1]]%*%t(C))^2)/ssy)
            R2Y<-c(R2Y,1-sum(Fxy^2)/ssy)
            U<-Us[[i+1]]
        }
    }
    
    Bus<-NULL
    Bts<-NULL
    for(i in seq_len((oax+1)))
    {
        Bustemp <- NULL
        Btstemp <- NULL
        for(j in seq_len((oay+1)))
        {
            Bustemp[[j]]<- solve(t(Us[[j]])%*%Us[[j]])%*%t(Us[[j]])%*%Ts[[i]]
            Btstemp[[j]]<- solve(t(Ts[[i]])%*%Ts[[i]])%*%t(Ts[[i]])%*%Us[[j]]
        }
        Bus[[i]] <- Bustemp
        Bts[[i]] <- Btstemp
    }
    
    Bu<- solve(t(U)%*%U)%*%t(U)%*%TT
    Bt<- solve(t(TT)%*%TT)%*%t(TT)%*%U
    
    
    R2Yhat<-NULL
    R2Xhat<-NULL
    
    
    #This is OC style - keeping for now... :
    for(i in seq_len((oax+1)))
    {
        BtTmp<-solve(t(Ts[[i]])%*%Ts[[i]])%*%t(Ts[[i]])%*%U
        YhatTmp<-Ts[[i]]%*%BtTmp %*%t(C)
        R2Yhat<-c(R2Yhat, 1-sum((YhatTmp-Y)^2)/ssy)
    }
    
    for(i in seq_len((oay+1)))
    {
        BuTmp<-solve(t(Us[[i]])%*%Us[[i]])%*%t(Us[[i]])%*%TT
        XhatTmp<-Us[[i]]%*% BuTmp %*% t(W)
        R2Xhat<-c(R2Xhat, 1-sum((XhatTmp-X)^2)/ssx)
    }
    
    model <-list(NULL)
    
    model$TT<-TT
    model$Ts<-Ts
    model$W<-W
    model$Wo<-Wo
    model$Pyo<-Pyo
    model$To<-To
    
    model$U<-U
    model$Us<-Us
    model$C<-C
    model$Co<-Co
    model$Pxo<-Pxo
    model$Uo<-Uo
    
    model$Bt<-Bt
    model$Bu<-Bu
    
    
    model$Bts<-Bts
    model$Bus<-Bus
    
    model$R2X<-R2X
    model$R2Xcorr<-R2Xcorr
    model$R2Xo<-R2Xo
    model$R2Xhat<-R2Xhat
    
    model$R2Y<-R2Y
    model$R2Ycorr<-R2Ycorr
    model$R2Yo<-R2Yo
    model$R2Yhat<-R2Yhat
    
    model$ssx<-ssx
    model$ssy<-ssy
    
    
    return(model)
}


o2pls_m<-function(X, Y, prep,nc,ncox,ncoy,nrcv = 0, r2cutoff = 0.3, pCutoff = 0.05)
{
    # O2-PLS function used to build the model
    # For further predictions use the o2pls_pred function
    #
    # Input:
    # X - X matrix input
    # Y - Y matrix input
    # prep - preprocessing available
    #           - 'no' : no preprocessing
    #           - 'mc' : meancentering
    #           - 'uv' : Univariance scaling
    #
    # nc - number of correlated components, it can be determinated by PCA of Y'X
    # ncox - number of orthogonal components in X
    # ncoy - number of orthogonal component in Y
    # nrcv - number of fold in the cross-validation (full cross validation)
    #
    ## written by Dr Jie Hao & Dr Xin Zou (2016), Shanghai Jiaotong University, China
    
    
    opt <- 1
    idx <- NULL
    
    muClass1 <- colMeans(X[which(Y[,1]==1),])
    muClass2 <- colMeans(X[which(Y[,2]==1),])
    Signs <- sign(muClass2-muClass1)
    
    
    tmp <-dim(X)
    nsx <- tmp[1]
    nvx <- tmp[2]
    tmp <- dim(Y)
    nsy <- tmp[1]
    nvy <- tmp[2]
    
    model <- list(NULL)
    model$ns <- nsx
    
    model$nc<-nc
    model$ncox<-ncox
    model$ncoy<-ncoy
    
    if (nsx!=nsy)
    {
        model<-NULL
        stop('Number of samples are different in X and Y\n')
    }
    
    model$MeanX<-colMeans(X)
    model$MeanY<-colMeans(Y)
    
    model$StdX<-apply(X,2,sd)
    model$StdY<-apply(Y,2,sd)
    
    model$SSX<-sum(X^2)
    model$SSY<-sum(Y^2)
    model$CSX<-colSums(X^2)
    model$CSY<-colSums(Y^2)
    # Preprocessing
    
    model$modelType <- 'pls'
    
    if (nrcv == 0) # no CV
    {
        if (prep == 'no')
        {
            model$preprocessing<-'no'
        } else if (prep == 'mc')
        {
            #disp('Meancentering')
            model$preprocessing<-'mc'
            X<-X-kronecker(matrix(1,nsx,1), t(model$MeanX))
            Y<-Y-kronecker(matrix(1,nsy,1), t(model$MeanY))
        } else if (prep == 'uv')
        {
            # disp('Univariance scaling')
            model$preprocessing<-'uv'
            X<-X-kronecker(matrix(1,nsx,1), t(model$MeanX))
            Y<-Y-kronecker(matrix(1,nsy,1), t(model$MeanY))
            X<-X/kronecker(matrix(1,nsx,1), t(model$StdX))
            Y<-Y/kronecker(matrix(1,nsy,1), t(model$StdY))
        } else if (prep == 'pa')
        {
            model$preprocessing<-'pa'
            X<-X-kronecker(matrix(1,nsx,1), t(model$MeanX))
            Y<-Y-kronecker(matrix(1,nsy,1), t(model$MeanY))
            X<-X/kronecker(matrix(1,nsx,1), t(sqrt(model$StdX)))
            Y<-Y/kronecker(matrix(1,nsy,1), t(sqrt(model$StdY)))
        } else
        {
            model<-NULL
            stop('Unknown Preprocessing\n')
        }
        
        X[is.nan(X)] <- 0
        pax =nc; oax = ncox; oay = ncoy;splitSize=300
        M<-mjrO2pls(X,Y,nc,ncox,ncoy)
        M$preprocessing <- model$preprocessing
        #    M$W <- M$W
        M$ns <- model$ns
        
        M$nc <- model$nc
        M$ncox <- model$ncox
        M$ncoy <- model$ncoy
        
        M$MeanX <- model$MeanX
        M$MeanY <- model$MeanY
        
        M$StdX <- model$StdX
        M$StdY <- model$StdY
        M$X_out <- X
        
        M$SSX <- model$SSX
        M$SSY <- model$SSY
        
        M$CSX <- model$CSX
        M$CSY <- model$CSY
        M$Q2Yhatcum <- NULL
        model <- M
        
        if (is.null(dim(model$W)[2]))
        {
            nv <- 1
        } else
        {
            nv <-dim(model$W)[2]
        }
        
        
        if (prep == 'uv')
        {
            WC <- t(model$W)
        } else if (prep == 'pa')
        {
            for (k in seq_len(nv))
            WC[k,] <- t(model$W[,k])/sqrt(model$StdX)
        } else if (prep == 'mc')
        {
            for (k in seq_len(nv))
            WC[k,] <- t(model$W[,k])/model$StdX
        }
        
        tmp <-dim(X)
        nsx <- tmp[1]
        nvx <- tmp[2]
        
        buffer <- kronecker(matrix(1,1,nvx),model$TT[,1])
        CC <- colMeans((buffer-kronecker(matrix(1,nsx,1),t(colMeans(buffer))))*(X-kronecker(matrix(1,nsx,1),t(colMeans(X)))))/(nsx-1)*nsx
        s1 <- sd(model$TT[,1])
        S1 <- matrix(s1,1,nvx)
        S2 <- apply(X,2,sd)
        median_R <- CC/(S1*S2)
        median_R[is.nan(median_R)]<-0
        median_P <- tTest(median_R,nsx)
        
        IX_P <- which(median_P<pCutoff)
        
        if (prep == 'uv')
        {
            W2<-WC*WC
            for (k in seq_len(nv))
            W2[k,] <- W2[k,]/norm(as.matrix(W2[k,]), 'i')
        } else if (prep == 'pa')
        {
            for (k in seq_len(nv))
            {
                W2[k,]<-WC[k,]*WC[k,]
                W2[k,] <- W2[k,]/norm(as.matrix(W2[k,]), 'i')
            }
        } else if (prep == 'mc')
        {
            for (k in seq_len(nv))
            {
                W2[k,]<-WC[k,]*WC[k,]
                W2[k,] <- W2[k,]/norm(as.matrix(W2[k,]), 'i')
            }
        }
        median_W2_bt <- W2[1,]
        
        
        #Bool = exist('opt', 'var')
        #if ~isempty(pCutoff) && ~isempty(r2cutoff)
        if (opt == 0)
        {
            idx <- which(median_W2_bt>r2cutoff)
        } else
        {
            II <- which(median_W2_bt>r2cutoff)
            idx <- intersect(IX_P, II)
        }
        
        ## save output
        model$sig_idx <- idx   # the significant ppm points
        model$R <- median_R
        model$P <- median_P
        #     model$WC_unscaled <- WC.*model$StdX
        model$W2 <- median_W2_bt
        model$WC <- WC
        model$signs <- Signs
    }
    # else # CV
    # {
    #   # nrcv-fold cross validation
    #   block_num <- floor(nsx/nrcv)
    #   #     Q2Yhatcum <- zeros(1,nrcv)
    #
    #   for (cv in 1:nrcv) # nrcv iterations of CV
    #   {
    #     idx_test <- nrcv*((1:block_num)-1)+cv
    #     idx_tr <- 1:nsx
    #     idx_tr[idx_test] <- NULL
    #     X_test <- X[idx_test,]
    #     Y_test <- Y[idx_test,]
    #     X_tr <- X[idx_tr,]
    #     Y_tr <- Y[idx_tr,]
    #
    #     nsx_tr<-dim(X_tr)[1]
    #     nvx_tr<-dim(X_tr)[2]
    #     nsy_tr<-dim(Y_tr)[1]
    #     nvy_tr<-dim(Y_tr)[2]
    #     nsx_test<-dim(X_test)[1]
    #     nvx_test<-dim(X_test)[2]
    #     nsy_test<-dim(Y_test)[1]
    #     nvy_test<-dim(Y_test)[2]
    #
    #     if (prep == 'no')
    #     {
    #       model$preprocessing <-'no'
    #     } else if (prep == 'mc')
    #     {
    #       #                 disp('Meancentering')
    #       model$preprocessing <-'mc'
    #       X_tr<-X_tr-kronecker(matrix(1,nsx_tr,1),t(mean(X_tr)))
    #       Y_tr<-Y_tr-kronecker(matrix(1,nsy_tr,1),t(mean(Y_tr)))
    #       X_test<-X_test-kronecker(matrix(1,nsx_test,1),t(mean(X_test)))
    #       Y_test<-Y_test-kronecker(matrix(1,nsy_test,1),t(mean(Y_test)))
    #     } else if (prep == 'uv')
    #     {
    #       # disp('Univariance scaling')
    #       model$preprocessing <-'uv'
    #       X_tr<-X_tr-kronecker(matrix(1,nsx_tr,1),t(mean(X_tr)))
    #       Y_tr<-Y_tr-kronecker(matrix(1,nsy_tr,1),t(mean(Y_tr)))
    #       X_tr<-X_tr/kronecker(matrix(1,nsx_tr,1),t(sd(X_tr)))
    #       Y_tr<-Y_tr/kronecker(matrix(1,nsy_tr,1),t(sd(Y_tr)))
    #
    #       X_test<-X_test-kronecker(matrix(1,nsx_test,1),t(mean(X_test)))
    #       Y_test<-Y_test-kronecker(matrix(1,nsy_test,1),t(mean(Y_test)))
    #       X_test<-X_test/kronecker(matrix(1,nsx_test,1),t(sd(X_test)) )
    #       Y_test<-Y_test/kronecker(matrix(1,nsy_test,1),t(sd(Y_test)))
    #     } else if (prep == 'pa')
    #     {
    #       model$preprocessing <-'pa'
    #       X_tr<-X_tr-kronecker(matrix(1,nsx_tr,1),t(mean(X_tr)))
    #       Y_tr<-Y_tr-kronecker(matrix(1,nsy_tr,1),t(mean(Y_tr)))
    #       X_tr<-X_tr/kronecker(matrix(1,nsx_tr,1),t(sqrt(sd(X_tr))))
    #       Y_tr<-Y_tr/kronecker(matrix(1,nsy_tr,1),t(sqrt(sd(Y_tr))) )
    #
    #       X_test<-X_test-kronecker(matrix(1,nsx_test,1),t(mean(X_test)))
    #       Y_test<-Y_test-kronecker(matrix(1,nsy_test,1),t(mean(Y_test)) )
    #       X_test<-X_test/kronecker(matrix(1,nsx_test,1),t(sqrt(sd(X_test))))
    #       Y_test<-Y_test/kronecker(matrix(1,nsy_test,1),t(sqrt(sd(Y_test))))
    #     } else
    #     {
    #       model<-NULL
    #       stop('Unknown Preprocessing\n')
    #     }
    #     # training
    #     model<-mjrO2pls(X_tr,Y_tr,nc,ncox,ncoy,'standard')
    #     #yhat
    #     modelPredy<-mjrO2plsPred(X_test, Y_test, model,ncox,ncoy,'x')
    #     #         #xhat
    #     #         modelPredx<-mjrO2plsPred(cvSet$xTest, cvSet$yTest, model,ioax-1,ioay-1,'y')
    #
    #     # # the overall Q2
    #     SSY<-sum(Y_test^2)
    #     Q2Yhatcum[cv] <- 1-sum((modelPredy$Yhat-Y_test)*(modelPredy$Yhat-Y_test))/SSY
    #   }
    #   model$Q2Yhatcum <- colMeans(Q2Yhatcum(!is.na(Q2Yhatcum)))
    # }
    return (model)
}

OPLSDA <- function (data1, data2, ppm, oplsda_para)
{#  model
    # run OPLS-DA and display results
    # Input:
    # "data1" and "data2" are groups of spectra from two biological classes
    # "ppm" is the ppm variables
    # "oplsda_para" - parameters
    
    # Output:
    # "model" contains all OPLS-DA results, where the discriminatory variables are in the model.sig_idx.
    # To extract the ppm of biomarker variables, type in
    # 'biomarker_ppm=ppm(model.sig_idx)' in the command window.
    #
    # The color loading plot shows the discriminatory variables. The red color in the
    # loading plot means the variables have more contribution to the
    # discrimination of classes.
    
    # The interpretation of the red and blue dots above the loading corresponds to:
    # blue dot - upward loading peak - up-regulate in data2 compared to data1.
    # red dot - downward loading peak - down-regulate in data2 compared to data1.
    
    # The parameters are set in the OPLSDA para.txt file
    # nrcv - number of folds cross-validation
    # nc - number of correlated variables in X
    # ncox - number of orthogonal variables in X
    # ncoy - number of orthogonal variables in Y
    # r2cutoff - r2 cutoff value to identify discriminatory variables
    # p-cutoff - p cutoff value to identify discriminatory variables
    # np - number of permutation tests
    # preprocessing - pre processing methods, 'uv', 'pa', ''mc'
    #
    ## written by Dr Jie Hao & Dr Xin Zou (2016), Shanghai Jiaotong University, China
    
    
    nc <- oplsda_para$nc
    ncox <- oplsda_para$ncox
    ncoy <- oplsda_para$ncoy
    r2cutoff <- oplsda_para$r2cutoff
    pCutoff <- oplsda_para$pcutoff
    prep <- oplsda_para$preprocessing
    
    
    pCutoff <- pCutoff/length(ppm)
    
    ## format data
    X <- rbind(data1,data2)
    Y <- matrix(0, dim(X)[1],2)
    Y[seq_len(dim(data1)[1]),1] <- 1
    Y[(dim(data1)[1]+1):dim(Y)[1],2] <- 1
    
    ##
    model <- o2pls_m(X,Y,prep,nc,ncox,ncoy,0, r2cutoff, pCutoff)
    return(model)
    ##
}


tTest <- function(r, n)
{
    #
    ## written by Dr Jie Hao & Dr Xin Zou (2016), Shanghai Jiaotong University, China
    p <- pbeta((1-r^2),(n/2-1),1/2)
    return (p)
}
