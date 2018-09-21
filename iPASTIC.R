##########################
##
##   iPASTIC: online plant abiotic stress index calculator. 
##   
##   Authors:
##   Alireza Pour-Aboughadareh (a.poraboghadareh@gmail.com)
##   Mohsen Yousefian (contact@mohsenyousefian.com)
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
##   You can visualize correlation matrixes using tools such as corrplot or ggcorrplot
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
    K1STI <- function(Ys, Yp, YsBar, YpBar)
    {
        return((Yp / YpBar) * STI(Ys, Yp, YsBar, YpBar))
    }
    K2STI <- function(Ys, Yp, YsBar, YpBar)
    {
        return((Ys / YsBar) * STI(Ys, Yp, YsBar, YpBar))
    }
    YI <- function(Ys, Yp, YsBar, YpBar)
    {
        return(Ys / YsBar)
    }
    YSI <- function(Ys, Yp, YsBar, YpBar)
    {
        return(Ys / Yp)
    }
    RDI <- function(Ys, Yp, YsBar, YpBar)
    {
        return((Ys / Yp) / (YsBar / YpBar))
    }
    SSPI <- function(Ys, Yp, YsBar, YpBar)
    {
        return((Yp - Ys) / (2 * YpBar) * 100)
    }
    ATI <- function(Ys, Yp, YsBar, YpBar)
    {
        return((Yp - Ys) / (YpBar / YsBar) * (sqrt(Yp * Ys)))
    }
    SNPI <- function(Ys, Yp, YsBar, YpBar)
    {
        return(((Yp + Ys) / (Yp - Ys)) ^ (1 / 3) * (Yp * Ys * Ys) ^ (1 / 3))
    }
    getranks_df <- function(df_orig)
    {
        descendings <- c(2,3,6,7,8,10,11,12,13,14,15,18)
        for (col in descendings)
            df_orig[col] = df_orig[col]* - 1

        df <- t(df_orig[, c(2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)])
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
        stats_df <- data.frame(table_original[, 1], table_original[, 2], table_original[, 3], runFunc(RC), runFunc(TOL), runFunc(MP), runFunc(GMP), runFunc(HM), runFunc(SSI), runFunc(STI), runFunc(K1STI), runFunc(K2STI), runFunc(YI), runFunc(YSI), runFunc(RDI), runFunc(SSPI), runFunc(ATI), runFunc(SNPI))
        colnames(stats_df) <- c("Species", "Yp", "Ys", "RC", "TOL", "MP", "GMP", "HM", "SSI", "STI", "K1STI", "K2STI", "YI", "YSI", "RDI", "SSPI", "ATI", "SNPI")

        ranks_df <- getranks_df(stats_df)

        output <- list(indices = stats_df, ranks = ranks_df,correlations=list(pearson=cor(data.matrix(results$indices[, 2:length(results$indices)])),spearman=cor(data.matrix(results$ranks[, 2:(length(results$ranks) - 3)]))))
        return(output)
    }

})()