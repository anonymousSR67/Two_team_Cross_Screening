library(sensitivitymw)
library(sensitivitymv)
library(DOS)
library(haven)
library(caret)
setwd("C:/Study_2023/Cross_Screen/Data_storage")


#==This code reproduces all the p-values reported in the left panel of Table D.5.
#== The analysis plan prepared by team A (after they explored the Catholics data) is now being 
#== implemented on the non-Catholics data. 





#=We read four datasets for this exercise:

#== (i) Master data file (from the WLS) 

#== (ii) Matched pairs index for the non-Catholics 
#(we got this file after running our Risk-set matching algorithm on the master file).

#== (iii) Pregnancy intention information for the treated non-Catholics - extracted from the master file 
#== (iv) Pregnancy intention information for the control non-Catholics - extracted from the master file
 





data <- read_dta(file = "WLS_master_data.dta")
cath_match <- read.csv("non_Catholics_matched_pairs.csv")
trt_child <- read.csv("Pregnancy_Intention_info_for_treated_non_Catholics.csv")
cont_child <- read.csv("Pregnancy_Intention_info_for_control_non_Catholics.csv")






#== Extracting information for only the female graduates (as our study focuses on them only) 

data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: female




#Now we follow the steps mentioned in the left panel of Figure 6.1 and obtain the p-values one by one.  

#============Step 1: test CESD===========

vv <- c("z_mu001rec")
col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}
SP <- numeric()
for(i in 1: nrow(cath_match))
{
  t1 <- data_f[,col][data_f$idpub==cath_match[i,2]]
  t2 <- data_f[,col][data_f$idpub==cath_match[i,3]]
  temp <- c(t1,as.numeric(cath_match[i,c(2,3)]),t2)
  SP <- rbind(SP,temp)
}
SP[,1][SP[,1]<0]<- NA
SP[,4][SP[,4]<0]<- NA


senmwCI(SP[,1]-SP[,4], gamma = 1, one.sided = F)




#==========Step 2: test Self-Acceptance=====================

vv <- c("z_mn046rec")

col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}


SP <- numeric()
for(i in 1: nrow(cath_match))
{
  t1 <- data_f[,col][data_f$idpub==cath_match[i,2]]
  t2 <- data_f[,col][data_f$idpub==cath_match[i,3]]
  temp <- c(t1,as.numeric(cath_match[i,c(2,3)]),t2)
  SP <- rbind(SP,temp)
}
SP[,1][SP[,1]<0]<- NA
SP[,4][SP[,4]<0]<- NA

senmwCI(SP[,1]-SP[,4], gamma = 1, one.sided = F)


#===============Step 3: Additional number of Children======

vv <- c("z_cmkdb1","z_cmkdb2",
        "z_cmkdb3","z_cmkdb4",
        "z_cmkdb5","z_cmkdb6",
        "z_cmkdb7","z_cmkdb8",
        "z_cmkdb9","z_cmkdb0")

col <- numeric()
for(kk in 1: length(vv))
{
  temp <- which(colnames(data_f)==vv[kk])
  col <- c(col,temp)
}

child_data <- numeric()
for(s in 1: length(cath_match[,1]))
{
  which_indx <- which(data_f$idpub==cath_match[s,2])
  trt_add_child <- as.numeric(data_f$z_kidsno[which_indx])- which(cath_match[s,1]==data_f[which_indx,col])[1]
  c1 <- c(which(cath_match[s,1]==data_f[which_indx,col])[1],trt_add_child)
  
  which_indx <- which(data_f$idpub==cath_match[s,3])
  cont_add_child <- data_f$z_kidsno[which_indx]- which(cath_match[s,4]==data_f[which_indx,col])[1]
  c2 <- c(cont_add_child,which(cath_match[s,4]==data_f[which_indx,col])[1])
  
  temp <- c(c1,c2)
  child_data <- rbind(child_data,temp)
  
}

cc <- numeric()
for(kk in 1:length(child_data[,1]))
{
  ind1 <- child_data[kk,1]
  c1 <- length(which(trt_child[kk,(2+ind1+1):ncol(trt_child)]==3))
  
  ind1 <- child_data[kk,4]
  c2 <- length(which(cont_child[kk,(2+ind1+1):ncol(cont_child)]==3))
  cc <- rbind(cc,c(c1,c2))
}
SP <- cc

senmwCI(SP[,1]-SP[,2], gamma = 1, one.sided = F)




#============= Step 4: Divorces=========================

aa <- which(colnames(data_f)=="z_cmslfm") 
bb <- which(colnames(data_f)=="z_cmslsm") 
cc <- which(colnames(data_f)=="z_cmsltm") 
marriage_data <- numeric()
marriage_data_2 <- numeric()
for(kk in 1: length(cath_match[,1]))
{
  #trtdata
  t1 <- data_f[data_f$idpub==cath_match[kk,2],c(aa,bb,cc)]
  t2 <- (cath_match[kk,1]-1900)*12+6
  t3 <- length(which(t1>t2))
  
  #Contdata
  t1 <- data_f[data_f$idpub==cath_match[kk,3],c(aa,bb,cc)]
  t2 <- (cath_match[kk,4]-1900)*12+6
  t4 <- length(which(t1>t2))
  
  marriage_data <- rbind(marriage_data,c(t3,t4))
}
SP <- marriage_data

senmwCI(SP[,1]-SP[,2], gamma = 1, one.sided = F)



#============== Step 5: No. of Job spells=========

emp_hist <- numeric()
for(i in 1: nrow(cath_match))
{
  w1 <- which(colnames(data_f)=="rf001j1c")
  w2 <- which(colnames(data_f)=="rf002j1c")
  t111 <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  t121 <- data_f[,w2][data_f$idpub==cath_match[i,2]]
  
  
  w1 <- which(colnames(data_f)=="rf001j2c")
  w2 <- which(colnames(data_f)=="rf002j2c")
  t112 <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  t122 <- data_f[,w2][data_f$idpub==cath_match[i,2]]
  
  w1 <- which(colnames(data_f)=="rf001j3c")
  w2 <- which(colnames(data_f)=="rf002j3c")
  t113 <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  t123 <- data_f[,w2][data_f$idpub==cath_match[i,2]]
  
  w1 <- which(colnames(data_f)=="rf001j4c")
  w2 <- which(colnames(data_f)=="rf002j4c")
  t114 <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  t124 <- data_f[,w2][data_f$idpub==cath_match[i,2]]
  
  w1 <- which(colnames(data_f)=="rf001jcd")
  w2 <- which(colnames(data_f)=="rf002jcd")
  t115 <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  t125 <- data_f[,w2][data_f$idpub==cath_match[i,2]]
  
  w1 <- which(colnames(data_f)=="rf006jsd")
  no_sp <- data_f[,w1][data_f$idpub==cath_match[i,2]]
  
  t1 <- c(t125,t115,t124,t114,t123,t113,t122,t112,t121,t111,(cath_match[i,1]-1900)*12+6, no_sp, cath_match[i,2])
  
  
  #================================================
  
  w1 <- which(colnames(data_f)=="rf001j1c")
  w2 <- which(colnames(data_f)=="rf002j1c")
  t111 <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  t121 <- data_f[,w2][data_f$idpub==cath_match[i,3]]
  
  
  w1 <- which(colnames(data_f)=="rf001j2c")
  w2 <- which(colnames(data_f)=="rf002j2c")
  t112 <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  t122 <- data_f[,w2][data_f$idpub==cath_match[i,3]]
  
  w1 <- which(colnames(data_f)=="rf001j3c")
  w2 <- which(colnames(data_f)=="rf002j3c")
  t113 <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  t123 <- data_f[,w2][data_f$idpub==cath_match[i,3]]
  
  w1 <- which(colnames(data_f)=="rf001j4c")
  w2 <- which(colnames(data_f)=="rf002j4c")
  t114 <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  t124 <- data_f[,w2][data_f$idpub==cath_match[i,3]]
  
  w1 <- which(colnames(data_f)=="rf001jcd")
  w2 <- which(colnames(data_f)=="rf002jcd")
  t115 <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  t125 <- data_f[,w2][data_f$idpub==cath_match[i,3]]
  
  w1 <- which(colnames(data_f)=="rf006jsd")
  no_sp <- data_f[,w1][data_f$idpub==cath_match[i,3]]
  
  t2 <- c(cath_match[i,3], no_sp,(cath_match[i,4]-1900)*12+6,t111,t121,t112,t122,t113,t123,t114,t124,t115,t125)
  emp_hist <- rbind(emp_hist,c(t1,t2))
}

aa <- which(is.na(emp_hist[,12])==T)
bb <- which(is.na(emp_hist[,15])==T)
ab <- union(aa,bb)
emp_hist <- emp_hist[-ab,]

aa <- which(emp_hist[,12]==-1)
bb <- which(emp_hist[,15]==-1)
ab <- union(aa,bb)
emp_hist <- emp_hist[-ab,]

aa <- which(emp_hist[,12]==0)
bb <- which(emp_hist[,15]==0)
ab <- union(aa,bb)
emp_hist <- emp_hist[-ab,]

SP <- emp_hist

senmwCI(SP[,12]-SP[,15], gamma = 1, one.sided = F)












