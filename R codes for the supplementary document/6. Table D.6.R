library(DOS)
library(haven)
library(caret)
library(xtable)
library(lsr)
library(ggplot2)
library(sensitivitymw)

setwd("C:/Study_2023/Cross_Screen/Data_storage")


cath_match <- read.csv("Catholics_matched_pairs.csv")
cath_match <- cath_match[,-c(1)]
trt_child <- read.csv("Pregnancy_Intention_info_for_treated_Catholics.csv")
cont_child <- read.csv("Pregnancy_Intention_info_for_control_Catholics.csv")
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
capture_m_vals <- numeric()
capture_method <- numeric()



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
box_dd <- data.frame(val=c(SP[,1], SP[,4]),class =c(rep("Trt",339),rep("Cont",339)), variable=rep("CESD",678))

#==== Wilcoxon===
tt <- wilcox.test(SP[,1],SP[,4], paired = T, alternative="greater")
wp <- tt$p.value

#===== U Statistics====
temp <- SP[,1]-SP[,4]
temp <- temp[-which(is.na(temp)==T)]
testing <- senU(temp,gamma=1,m=8,m1=5,m2=8)
pv1 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=6,m2=7)
pv2 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=7,m2=8)
pv3 <- testing$pval

cesd_res <- c(wp,pv1,pv2,pv3)



#======== Self Acceptance============

vv <- c("z_rn012red")

col <- numeric()
for(k in 1: length(vv))
{
  col <- c(col, which(colnames(data_f)==vv[k]))
}


SP <- numeric()
for(i in 1: nrow(cath_match))
{
  t1 <- data_f[,col[1]][data_f$idpub==cath_match[i,2]]
  t2 <- data_f[,col[1]][data_f$idpub==cath_match[i,3]]
  temp <- c(t1,as.numeric(cath_match[i,c(2,3)]),t2)
  SP <- rbind(SP,temp)
}
SP[,1][SP[,1]<0]<- NA
SP[,4][SP[,4]<0]<- NA

#==== Wilcoxon===
tt <- wilcox.test(SP[,1],SP[,4], paired = T, alternative="less")
wp <- tt$p.value

#===== U Statistics====
temp <- -SP[,1]+SP[,4]
temp <- temp[-which(is.na(temp)==T)]
testing <- senU(temp,gamma=1,m=8,m1=7,m2=8)
pv1 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=6,m2=7)
pv2 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=5,m2=8)
pv3 <- testing$pval

self_accep_res <- c(wp,pv1,pv2,pv3)


#=================== Additional number of Children======

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
#==== Wilcoxon======
tt <- wilcox.test(SP[,1],SP[,2], paired = T, alternative="greater")
wp <- tt$p.value

#===== U Statistics====
temp <- SP[,1]-SP[,2]
testing <- senU(temp,gamma=1,m=8,m1=5,m2=8)
pv1 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=6,m2=7)
pv2 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=7,m2=8)
pv3 <- testing$pval

add_child_res <- c(wp,pv1,pv2,pv3)

#=================== Additional number of Children======

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
  tt1 <- c(t2,t1,t3)
  #Contdata
  t1 <- data_f[data_f$idpub==cath_match[kk,3],c(aa,bb,cc)]
  t2 <- (cath_match[kk,4]-1900)*12+6
  t4 <- length(which(t1>t2))
  tt2 <- c(t2,t1,t4)
  marriage_data <- rbind(marriage_data,c(t3,t4))
  marriage_data_2 <- rbind(marriage_data_2,c(tt1,tt2))
}
SP <- marriage_data
#==== Wilcoxon======
tt <- wilcox.test(SP[,1],SP[,2], paired = T, alternative="greater")
wp <- tt$p.value

#===== U Statistics====
temp <- SP[,1]-SP[,2]
testing <- senU(temp,gamma=1,m=8,m1=5,m2=8)
pv1 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=6,m2=7)
pv2 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=7,m2=8)
pv3 <- testing$pval

marriage_res <- c(wp,pv1,pv2,pv3)



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
#==== Wilcoxon======
tt <- wilcox.test(SP[,12],SP[,15], paired = T, alternative="greater")
wp <- tt$p.value

#===== U Statistics====
temp <- SP[,12]-SP[,15]
testing <- senU(temp,gamma=1,m=8,m1=5,m2=8)
pv1 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=6,m2=7)
pv2 <- testing$pval
testing <- senU(temp,gamma=1,m=8,m1=7,m2=8)
pv3 <- testing$pval

job_res <- c(wp,pv1,pv2,pv3)

res <- rbind(cesd_res,self_accep_res,add_child_res,marriage_res, job_res)

print("Hence p-values in Table D.6 are as follows:")
print(res)





