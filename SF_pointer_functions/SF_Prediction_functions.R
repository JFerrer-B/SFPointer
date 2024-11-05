#' Splicing Factor Prediction
#' 
#' Methodology to predict context-specific splicing factors
#' 
#' @param P_value_PSI A data.frame with the p.values of the experiment.
#' @param ExS The ExS matrix biuldt in CreateExSmatrix function.
#' @param nSel Top ranked events to be considered as spliced events.
#' @param significance Threshold of P.value to consider which events are 
#' deferentially spliced. A vector of length equal to the number of contrasts.
#' If null it will consider the nSel top ranked events.
#' @param method methodology to apply: "Fisher" for Fisher's exact test 
#' (default) or "PoiBin" for Poisson Binomial test.
#' 
#' @return The function returs a list. This list has for each contrast a data.frame containing
#' the results of the prediction.
#' 
#' 
#' 
#' @import glmnet
#' @import poibin
#' @import Matrix
#' @importFrom stats binomial phyper qhyper pnorm dnorm gaussian coef
#' @importFrom speedglm control is.sparse
#' @importFrom IRanges IRanges

SF_Prediction <- function(P_value_PSI,ExS,nSel=1000,significance=NULL,
                              method="Fisher", valueRanking = "pvalue",PSI_table = NULL){
  
  
  resPred <- vector(mode="list",length = ncol(P_value_PSI))
  # P_value_PSI <- P_value_PSI[rownames(P_value_PSI) %in% rownames(ExS),,drop=FALSE]
  # to_add <- which(!rownames(ExS) %in% rownames(P_value_PSI))
  # to_add_matrix <- matrix(1,nrow=length(to_add),ncol=ncol(P_value_PSI))
  # rownames(to_add_matrix) <- rownames(ExS)[to_add]
  # colnames(to_add_matrix) <- colnames(P_value_PSI)
  # P_value_PSI <- rbind(P_value_PSI,to_add_matrix)
  ExS <- ExS[rownames(P_value_PSI), ]
  N <- nrow(ExS) #the same for each Fisher's test
  
  switch(method,
         Fisher={
           resPred <- hyperGeometricApproach(ExS, nSel, P_value_PSI, significance, resPred,N)
         },
         PoiBin={
           resPred <- poissonBinomialApproach(ExS, nSel, P_value_PSI, significance, resPred,N)
         },
         Wilcoxon={
           # nmTopEv <- significanceFunction (P_value_PSI, nSel=NULL, significance)
           # ExS <- ExS[nmTopEv, ]
           # resPred <- Wilcoxon.z.matrix(ExprT = t(abs(P_value_PSI[,4])),GeneGO = ExS)
           
           if (valueRanking == "PSI") {
             resPred <- WilcoxonApproach(P_value_PSI, ExS, significance, resPred, PSI_table)
           }else{
             resPred <- WilcoxonApproach(P_value_PSI, ExS, significance, resPred)
           }
           
         },
         Gsea ={
           if (valueRanking == "PSI") {
             resPred <- GseaApproach(P_value_PSI,ExS, significance, resPred, PSI_table)
           }else{
             resPred <- GseaApproach(P_value_PSI,ExS, significance, resPred)
           }
         } )
  
  return(resPred)
}




hyperGeometricApproach <- function(ExS, nSel, P_value_PSI, significance, resPred,N){
  
  for(cSel in 1:ncol(P_value_PSI)){
    
    nmTopEv <- significanceFunction (P_value_PSI, cSel, nSel, significance)
    nSel <- length(nmTopEv)
    hyperM <- hyperMatrixRes(cSel, nSel, ExS, P_value_PSI , significance,N)
    
    mynselevents <- matrix(0,nrow=nrow(ExS),ncol=1)
    rownames(mynselevents) <- rownames(ExS)
    mynselevents[nmTopEv,1] <- 1
    
    mix <- as.numeric(t(ExS) %*% mynselevents)
    # identical(hyperM$RBP,rownames(mix))
    mid <- hyperM$nHits
    miN_D <- N-mid
    hyperM[, "Pvalue_hyp_PSI"] <- myphyper(mix, mid, miN_D, nSel, lower.tail = FALSE)
    hyperM$x <- mix
    hyperM$qhyp_0.5 <- qhyper(0.5, mid, miN_D, nSel, lower.tail = FALSE)
    hyperM$Fc <- mix/hyperM$qhyp_0.5
    
    hyperM <- hyperM[order(hyperM$Pvalue_hyp_PSI), ]
    resPred[[cSel]]<-hyperM
    
  }
  
  return(resPred)
}


poissonBinomialApproach <- function(ExS, nSel, P_value_PSI, significance, resPred,N){
  myP <- getpij(ExS)

  for(cSel in 1:ncol(P_value_PSI)){
    
    nmTopEv <- significanceFunction(P_value_PSI, cSel, nSel, significance)
    nSel <- length(nmTopEv)
    hyperM <- hyperMatrixRes(cSel, nSel, ExS, P_value_PSI , significance,N)
    
    mynselevents <- matrix(0,nrow=nrow(ExS),ncol=1)
    rownames(mynselevents) <- rownames(ExS)
    mynselevents[nmTopEv,1] <- 1
    
    mix <- as.numeric(t(ExS) %*% mynselevents)
    
    myP2 <- myP[which(rownames(myP)%in%nmTopEv),]
    muk <- colSums(myP2)
    QQ <- 1 - myP2
    sigmak <- sqrt(diag(t(myP2)%*%QQ))
    MM <- QQ*(1-2*myP2)
    gammak <- diag(t(myP2)%*%MM)
    ind = gammak/(6 * sigmak^3 )
    
    kk1 = (mix + 0.5 - muk)/sigmak
    
    vkk.r = pnorm(kk1) + ind * (1 - kk1^2) * dnorm(kk1)
    vkk.r[vkk.r < 0] = 0
    vkk.r[vkk.r > 1] = 1
    
    hyperM$Pvalue_hyp_PSI <- 1-vkk.r
    
    weit_correction <- colSums(ExS)/colSums(ExS)
    if(any(is.na(weit_correction))){
      hyperM$Pvalue_hyp_PSI[which(is.na(weit_correction))] <- 1
    }
    
    hyperM$qhyp_0.5 <- t(myP) %*% mynselevents
    
    hyperM$x <- mix
    # hyperM$qhyp_0.5 <- qhyper(0.5, mid, miN_D, nSel, lower.tail = FALSE)
    
    hyperM$Fc <- mix/hyperM$qhyp_0.5
    if(any(is.na(hyperM$Pvalue_hyp_PSI))){
      hyperM$Pvalue_hyp_PSI[is.na(hyperM$Pvalue_hyp_PSI)] <- 1
    }
    hyperM <- hyperM[order(hyperM$Pvalue_hyp_PSI), ]
    resPred[[cSel]]<-hyperM
    
  }
  return(resPred)
}


significanceFunction <- function(P_value_PSI, cSel, nSel, significance ){
  
  # cSel <- 1
  
  if(is.null(significance)){
    if (!is.null(nSel)) {
      nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
    }else{
      nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])
    }
    
    #name of the top nSel events differentially spliced in this contrast
  }else{
    if(is.numeric(significance)){
      
      if(is.na(significance[cSel])){
        warning(paste0("no threshold selected for contrast ",cSel,". Top nSel events taken"))
        if (!is.null(nSel)) {
          nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
        }else{
          nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])
        }
      }else{
        nmTopEv <- rownames(P_value_PSI)[which(P_value_PSI[,cSel] < significance[cSel])]
      }
      
    }else{
      stop("significance must be numeric")
    }
  }
  return(nmTopEv)
}


hyperMatrixRes <- function(cSel, nSel, ExS, P_value_PSI , significance,N){
  hyperM <- data.frame(RBP = colnames(ExS),
                       nHits = colSums(ExS),
                       Pvalue_hyp_PSI =NA,
                       N = N,
                       d = colSums(ExS),
                       n = nSel,
                       x = NA,
                       qhyp_0.5 = NA,
                       Fc = NA,
                       stringsAsFactors = FALSE)
  return(hyperM)
}

GseaApproach <- function(P_value_PSI,ExS, significance, resPred, PSI_table=NULL){
  # P_value_PSI <- ResultBootstrap_kallisto_2$Pvalues
  # ExS <- miExS
  if(is.null(significance)){
    significance <- rep(1.1,ncol(P_value_PSI))
    message("There is no value for the variable significance. By default, a value of 1 has been taken for each contrast. 
    For the wilcoxon and GSEA methods, the significance variable is used to filter which events should be taken into account.
    By default, these two methods do not consider a filter to determine which events to consider, so a value of 1 is defaulted to each contrast. ")
  }
  listRes <- vector(mode="list",length = ncol(ExS))
  names(listRes) <- colnames(ExS)
  for (RBP in colnames(ExS)) {
    events_associatedRBP <- which(ExS[,RBP] == 1)
    listRes[[RBP]] <- names(events_associatedRBP)
    # listRes <- append(listRes, list(RBP=))
  }
  
  for(cSel in 1:ncol(P_value_PSI)){
    nmTopEv <- significanceFunction(P_value_PSI, cSel, nSel=NULL, significance)
    if (!is.null(PSI_table)) {
      ranks <- abs(PSI_table[nmTopEv,cSel])
      fgseaRes <- fgsea(pathways = listRes,
                        stats    = ranks)
    }else{
      ranks <- 1-P_value_PSI[nmTopEv,cSel]
      # fgseaRes <- try(fgsea(pathways = listRes,stats    = ranks,minSize=0),silent = TRUE)
      fgseaRes <- fgsea(pathways = listRes,stats    = ranks,minSize=0)
    }
    fgseaRes <- as.data.frame(fgseaRes[,c(1:7)])
    colnames(fgseaRes)[1] <- "RBP"
    fgseaRes <- fgseaRes[order(fgseaRes$pval),]
    resPred[[cSel]]<-fgseaRes
  }
  
  return(resPred)
  
}

WilcoxonApproach <- function(P_value_PSI,ExS, significance, resPred, PSI_table=NULL){
  if(is.null(significance)){
    significance <- rep(1.1,ncol(P_value_PSI))
    message("There is no value for the variable significance. By default, a value of 1 has been taken for each contrast. 
    For the wilcoxon and GSEA methods, the significance variable is used to filter which events should be taken into account.
    By default, these two methods do not consider a filter to determine which events to consider, so a value of 1 is defaulted to each contrast. ")
  }
  for(cSel in 1:ncol(P_value_PSI)){
    nmTopEv <- significanceFunction(P_value_PSI, cSel, nSel=NULL, significance)
    setExS <- ExS[nmTopEv,]
    if (!is.null(PSI_table)) {
      tabla_final <- Wilcoxon.z.matrix(ExprT = t(abs(PSI_table[nmTopEv,cSel])),GeneGO = setExS)
      
    }else{
      tabla_final <- Wilcoxon.z.matrix(ExprT = t(abs(P_value_PSI[nmTopEv,cSel])),GeneGO = setExS)
    }
    tabla_final <- data.frame(RBP=rownames(tabla_final), z_score = tabla_final[,1])
    tabla_final$Pvalue_hyp_PSI <- pnorm(tabla_final$z_score,lower.tail=TRUE)
    tabla_final <- tabla_final[order(tabla_final$Pvalue_hyp_PSI),]
    resPred[[cSel]] <- tabla_final
  }
  return(resPred)
}

Wilcoxon.z.matrix <- function(ExprT, GeneGO, 
                              alternative = c("two.sided", "less", "greater"),
                              mu = 0, paired = FALSE, 
                              exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95) {
  
  require(matrixStats)
  # Transp_ENSTxGO <- Matrix::t(GeneGO)
  RExprT <- rowRanks(ExprT, preserveShape = T)
  
  Prod <- t(RExprT %*% GeneGO)
  
  # Calculate the number of elements "0"(ny) and "1"(nx)
  nx <-  as.vector(rep(1,nrow(GeneGO)) %*% GeneGO)
  ny <- nrow(GeneGO) -nx
  
  # Calculate the estimated variance
  Var <- (nx*ny*(nx+ny+1)/12)
  
  # Calculate the estimated mean
  media <- nx*(nx + ny + 1)/2
  
  # Calculate the standard desviation
  std <- sqrt(Var)
  
  z <- (Prod-media)/std
  
  return((z))
}

