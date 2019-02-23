##########################
##   
##   iPASTIC: an online toolkit to calculate plant abiotic stress indices
##   
##   Authors:
##   Alireza Pour-Aboughadareh (a.poraboghadareh@gmail.com)
##   Mohsen Yousefian (contact@mohsenyousefian.com)
##   
##   Source: https://github.com/pour-aboughadareh/ipastic
##   
##########################
##   
##   Usage:
##   
##   
##   1. Load table data to a dataframe variable named "df"
##   2. Results <- Calculate(df)
##   3. Print(results$indices)
##   4. Print(results$ranks)
##   5. Correlation matrixes: out$correlations$pearson and out$correlations$spearman
##   6. Principal component analysis results: out$pca (out$pca$correlation_based and out$pca$covariance_based)
##   
##   Note 1: You can visualize correlation matrixes using tools such as corrplot or ggcorrplot
##   Note 2: You can draw biplot for pca results using built-in function biplot or similar third-party libraries
##   
##########################

(function()
{
  RC <- function(Ys, Yp, YsBar, YpBar)
  {
    return(((Yp - Ys) / Yp) * 100)
  }
  TOL <- function(Ys, Yp, YsBar, YpBar)
  {
    return(Yp - Ys)
  }
  MP <- function(Ys, Yp, YsBar, YpBar)
  {
    return((Yp + Ys) / 2)
  }
  GMP <- function(Ys, Yp, YsBar, YpBar)
  {
    return(sqrt(Ys * Yp))
  }
  HM <- function(Ys, Yp, YsBar, YpBar)
  {
    return(2 * (Ys * Yp) / (Ys + Yp))
  }
  SSI <- function(Ys, Yp, YsBar, YpBar)
  {
    return((1 - Ys / Yp) / (1 - YsBar / YpBar))
  }
  STI <- function(Ys, Yp, YsBar, YpBar)
  {
    return((Ys * Yp) / YpBar ^ 2)
  }
  YI <- function(Ys, Yp, YsBar, YpBar)
  {
    return(Ys / YsBar)
  }
  YSI <- function(Ys, Yp, YsBar, YpBar)
  {
    return(Ys / Yp)
  }
  RSI <- function(Ys, Yp, YsBar, YpBar)
  {
    return((Ys / Yp) / (YsBar / YpBar))
  }
  getranks_df <- function(df_orig)
  {
    descendings <- c(2,3,6,7,8,10,11,12,13)
    for (col in descendings)
      df_orig[col] = df_orig[col]* - 1
    
    df <- t(df_orig[, c(2, 3, 5, 6, 7, 8,9, 10, 11, 12, 13)])
    ranks <- apply(df, 1, rank, ties.method = "min")
    
    SR <- data.frame(apply(ranks, 1, sum))
    colnames(SR) <- "SR"
    
    AR <- SR / length(ranks[1,])
    colnames(AR) <- "AR"
    
    STD <- data.frame(apply(ranks, 1, sd))
    colnames(STD) <- "Std."
    
    ranks <- cbind(df_orig[1], ranks, SR, AR, STD)
    
    return(ranks)
  }
  
  setNumberPrecisionToMatchExcel<-function(num)
  {
    return(signif(num,digits = 15))
  }
  getSignChangingVec<-function(Matr)
  {
    mat <- apply(as.matrix(Matr), 2, setNumberPrecisionToMatchExcel)
    signChangingFlags <- c()
    for (i in 1:ncol(mat))
    {
      x <- mat[, i]
      max <- NULL
      maxInd <- 0
      for (j in 1:length(x))
      {
        if (is.null(max) || abs(x[j]) - max > (5e-16) || (abs(abs(x[j]) - max) < (5e-16) & x[j] < 0))
        {
          max <- abs(x[j])
          maxInd <- j
        }
      }
      colsign <- sign(x[maxInd])
      signChangingFlags <- c(signChangingFlags, colsign)
    }
    return(signChangingFlags)
  }
  getpcLabels <- function(count)
  {
    labels<-c()
    for(i in 1:count)
    {
      labels<-c(labels,paste("PC",toString(i),sep=""))
    }
    return(labels)
  }
  getpca <- function(df,iscov)
  {
    genonames<-rownames(df)
    factornames<-colnames(df)
    
    eigenvals<-NULL
    eigenvecs<-NULL
    loadings<-NULL
    contributions<-NULL
    scores<-NULL
    
    scoresignchangevec<-NULL
    
    if(!iscov)
    {
      eig<-eigen(cor(df))
      pca <- prcomp(df, scale = TRUE)
      scoresignchangevec<-getSignChangingVec(as.matrix(pca$rotation))[col(pca$x)]
      scores <- pca$x*scoresignchangevec
      
    }
    else
    {
      eig<-eigen(cov(df))
      pca <- prcomp(cov(df), scale = TRUE)
      scores <- predict(prcomp(df, scale = FALSE),df)
      scoresignchangevec<-getSignChangingVec(as.matrix(prcomp(df, scale = FALSE)$rotation))[col(scores)]
      scores <- scores*scoresignchangevec
      
    }
    
    
    eigenvals<-eig$values
    eigenvecs<-eig$vectors
    
    loadingCalcEigenvals<-eigenvals
    loadingCalcEigenvals[loadingCalcEigenvals<0]<-NaN
    
    loadings<-t(t(eigenvecs)*sqrt(loadingCalcEigenvals))
    
    loadings<-loadings[,complete.cases(t(loadings))]
    
    
    eigvecs_p2 <- eigenvecs^2
    contributions <- (eigvecs_p2)/rowSums(eigvecs_p2)*100
    
    
    
    signChangingFlags <- getSignChangingVec(loadings)
    
    
    
    eigenvecs<-eigenvecs[,1:length(signChangingFlags)]
    scores <- scores[, 1:min(nrow(scores), length(signChangingFlags))]
    
    eigenvals<-t(as.data.frame(eigenvals))
    importancedf<-summary(pca)$importance[2:3,]
    eigenvals<-rbind(eigenvals[,1:min(ncol(importancedf),ncol(eigenvals))],importancedf[,1:min(ncol(importancedf),ncol(eigenvals))])
    
    colnames(eigenvals)<-getpcLabels(ncol(eigenvals))
    rownames(eigenvals)<-c("Eigenvalue","Variability (%)","Cumulative %")
    
    eigenvecs<-as.data.frame(eigenvecs)
    colnames(eigenvecs)<-getpcLabels(ncol(eigenvecs))
    rownames(eigenvecs)<-factornames[1:nrow(eigenvecs)]
    
    loadings<-as.data.frame(loadings)
    colnames(loadings)<-getpcLabels(ncol(loadings))
    rownames(loadings)<-factornames[1:nrow(loadings)]
    
    contributions<-as.data.frame(contributions)
    colnames(contributions)<-getpcLabels(ncol(contributions))
    rownames(contributions)<-factornames[1:nrow(contributions)]
    
    if(!iscov)
      pca <- prcomp(df, scale = TRUE)
    else
      pca <- prcomp(df, scale = FALSE)
    
    pca$rotation <- pca$rotation * getSignChangingVec(pca$rotation[ , ! apply( pca$rotation , 2 , function(x) all(is.na(x)) ) ])[col(pca$rotation)]
    pca$x<-pca$x*getSignChangingVec(pca$x[ , ! apply( pca$x , 2 , function(x) all(is.na(x)) ) ])[col(pca$x)]
    
    return(list(eigenvals = eigenvals, eigenvecs = eigenvecs * signChangingFlags[col(eigenvecs)], loadings = loadings * signChangingFlags[col(loadings)], contributions = contributions[, 1:length(signChangingFlags)], scores = scores, pca_obj=pca ))
  }
  enablePcaVarCountWarning <- FALSE
  Calculate <<- function(table_original)
  {
    table <- cbind(table_original)[, -1]
    a <- nrow(table)
    b <- ncol(table)
    
    
    data <- data.matrix(table)
    
    Yp <- table[1]
    Ys <- table[2]
    YpBar <- apply(Yp, 2, mean)
    YsBar <- apply(Ys, 2, mean)
    
    runFunc <- function(func)
    {
      return(func(Ys, Yp, YsBar, YpBar));
    }
    stats_df <- data.frame(table_original[, 1], table_original[, 2], table_original[, 3], runFunc(RC), runFunc(TOL), runFunc(MP), runFunc(GMP), runFunc(HM), runFunc(SSI), runFunc(STI), runFunc(YI), runFunc(YSI), runFunc(RSI))
    colnames(stats_df) <- c("Species", "Yp", "Ys", "RC", "TOL", "MP", "GMP", "HM", "SSI", "STI", "YI", "YSI", "RSI")
    
    ranks_df <- getranks_df(stats_df)
    
    pcadf<-stats_df[2:ncol(stats_df)][-3]
    
    rownames(pcadf)<-stats_df[,1]
    
      if (enablePcaVarCountWarning)
        if(ncol(pcadf)>nrow(pcadf))
        {
            #warning("There are more variables than observations in output for dataset provided for principal component analysis.")
        }
    
    correlation_based_pca<-getpca(pcadf,FALSE)
    covariance_based_pca<-getpca(pcadf,TRUE)
    
    rownames(correlation_based_pca$scores) <- stats_df[,1]
    rownames(covariance_based_pca$scores) <- stats_df[,1]
    
    rownames(stats_df) <- stats_df[,1]
    stats_df<-stats_df[,-1]
    
    rownames(ranks_df) <- ranks_df[,1]
    ranks_df<-ranks_df[,-1]
    
    output <- list(indices = stats_df, ranks = ranks_df,correlations=list(pearson=cor(data.matrix(stats_df[-3][, 1:length(stats_df[-4])])),spearman=cor(data.matrix(ranks_df[, 1:(length(ranks_df) - 3)]))),pca=list(correlation_based=correlation_based_pca, covariance_based=covariance_based_pca))
    return(output)
  }
})()
