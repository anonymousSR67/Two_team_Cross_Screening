library(haven)
library(caret)
library(xtable)
library(lsr)
library(ggplot2)
library(sensitivitymv)


setwd("C:/Study_2023/Cross_Screen/Data_storage")
cath_match <- read.csv("Catholics_matched_pairs.csv")
cath_match <- cath_match[,-c(1)]

data <- read_dta(file = "WLS_master_data.dta")
data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: female

#====================PWB=====================
vv <- c("z_rn002red","z_rn004red","z_rn006red","z_rn008red","z_rn010red","z_rn012red","z_rn014red")

col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}


results <- numeric()
for(l in 1: length(col))
{
  SP <- numeric()
  for(i in 1: nrow(cath_match))
  {
    t1 <- data_f[,col[l]][data_f$idpub==cath_match[i,2]]
    t2 <- data_f[,col[l]][data_f$idpub==cath_match[i,3]]
    temp <- c(t1,as.numeric(cath_match[i,c(2,3)]),t2)
    SP <- rbind(SP,temp)
  }
  SP[,1][SP[,1]<0]<- NA
  SP[,4][SP[,4]<0]<- NA
  
  tt <- wilcox.test(SP[,1],SP[,4], paired = T, alternative="less")
  summary_PWB <- c(mean(SP[,1], na.rm=T), mean(SP[,4], na.rm = T),tt$p.value)
  results <- rbind(results,summary_PWB)
}
PWB_92 <- results

#PWB_92 contains the numbers in Table D.2.
