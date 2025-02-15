library(haven)
library(caret)
library(xtable)
library(lsr)
library(ggplot2)
library(sensitivitymv)
library(sensitivitymw)



setwd("C:/Study_2023/Cross_Screen/Data_storage")

cath_match <- read.csv("non_Catholics_matched_pairs.csv")
data <- read_dta(file = "WLS_master_data.dta")
data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: female

vv <- c("z_mu001rec","z_mu003rer","z_mu004rer","z_mu005rer","z_mu006rer","z_mu007rer",
        "z_mu008rer","z_mu009rer","z_mu010rer","z_mu011rer","z_mu012rer",
        "z_mu013rer","z_mu014rer","z_mu015rer","z_mu016rer","z_mu017rer",
        "z_mu018rer","z_mu019rer","z_mu020rer","z_mu021rer","z_mu022rer")
col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}

#======================= First part of the table===========

#============CESD===========
l <- 1
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
tt <- senmw(SP[,1]-SP[,4],gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(SP[,1]-SP[,4],gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(SP[,1]-SP[,4],gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(SP[,1]-SP[,4],gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(SP[,1]-SP[,4],gamma = 1.2, method = "l")
p5 <- tt$pval
cesd_p_val <- c(p1,p2,p3,p4,p5)



#=======================Sub scale scores========

col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}

dat <- list()
ind_exc <- list()
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
  SP <- data.frame(SP)
  i1 <- which(SP[,1] %in% c(-3,-2, NA))
  i2 <- which(SP[,4] %in% c(-3,-2, NA))
  ii <- union(i1,i2)
  dat[[l]] <- SP
  ind_exc[[l]] <- ii
}

#==Depressed effect=====
which <- c(3,13,8,10,12) - rep(1,5)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]] + dd[[4]] + dd[[5]]
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "l")
p5 <- tt$pval
depressed_p_val <- c(p1,p2,p3,p4,p5)


#==Low positive effect=====

which <- c(16,6,9) - rep(1,3)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]]
ddd$X1 <- 21-ddd$X1
ddd$X4 <- 21-ddd$X4
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "l")
p5 <- tt$pval
low_positive_p_val <- c(p1,p2,p3,p4,p5)




#==Somatic=====

which <- c(15,14,17,20,21,22) - rep(1,6)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]] + dd[[4]] + dd[[5]] + dd[[6]]
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "l")
p5 <- tt$pval
somatic_p_val <- c(p1,p2,p3,p4,p5)



#===Interpersonal====

which <- c(7,11) - rep(1,2)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] 
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(ddd[,1]-ddd[,4],gamma = 1.2, method = "l")
p5 <- tt$pval
Interpersonal_p_val <- c(p1,p2,p3,p4,p5)


first_part_p_values <- rbind(cesd_p_val,depressed_p_val,low_positive_p_val,somatic_p_val, Interpersonal_p_val)




#===================== Second part of the table=================

vv <- c("z_mn001rec","z_mn010rec","z_mn019rec","z_mn028rec", "z_mn037rec", "z_mn046rec")

temp_1 <- rep(0,383)
temp_2 <- rep(0,383)

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
  
  
  temp_1 <- temp_1 + as.numeric(SP[,1])
  temp_2 <- temp_2 + as.numeric(SP[,4])
  
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "w")
  p1 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "t")
  p2 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "l")
  p3 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1.2, method = "t")
  p4 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1.2, method = "l")
  p5 <- tt$pval
  
  summary_PWB <- c(p1,p2,p3,p4,p5)
  
  results <- rbind(results,summary_PWB)
}

tt <- senmw(temp_2 - temp_1  ,gamma = 1, method = "w")
p1 <- tt$pval
tt <- senmw(temp_2 - temp_1,gamma = 1, method = "t")
p2 <- tt$pval
tt <- senmw(temp_2 - temp_1,gamma = 1, method = "l")
p3 <- tt$pval
tt <- senmw(temp_2 - temp_1,gamma = 1.2, method = "t")
p4 <- tt$pval
tt <- senmw(temp_2 - temp_1,gamma = 1.2, method = "l")
p5 <- tt$pval

overall_PWB <- c(p1,p2,p3,p4,p5)

Second_part_p_values <- rbind(overall_PWB,results)




#=================Third (last) part of the table=================

vv <- c("z_ixsf1rec", "z_ixsf2rec")


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
  
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "w")
  p1 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "t")
  p2 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1, method = "l")
  p3 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1.2, method = "t")
  p4 <- tt$pval
  tt <- senmw(SP[,4]-SP[,1],gamma = 1.2, method = "l")
  p5 <- tt$pval
  
  summary_SF12 <- c(p1,p2,p3,p4,p5)
  
  results <- rbind(results,summary_SF12)
}
Third_part_p_values <- results


