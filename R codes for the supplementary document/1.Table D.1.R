library(haven)
library(caret)
library(xtable)
library(lsr)
library(ggplot2)
library(sensitivitymv)
library(here)

cath_match <- read.csv(here("Data", "Catholics_matched_pairs.csv"))
cath_match <- cath_match[,-c(1)]
data <- read_dta(here("Data", "WLS_master_data.dta"))
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


tt <- wilcox.test(SP[,1],SP[,4], paired = T, alternative="greater")
cesd__val <- c(mean(SP[,1], na.rm = T), mean(SP[,4], na.rm=T),tt$p.value)
  
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

tt <- wilcox.test(ddd[,1],ddd[,4], paired = T, alternative="greater")

depressed__val <- c(mean(ddd[,1],na.rm = T),mean(ddd[,4],na.rm = T), tt$p.value)


#==Low positive effect=====

which <- c(16,6,9) - rep(1,3)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]]
ddd$X1 <- 21-ddd$X1
ddd$X4 <- 21-ddd$X4
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- wilcox.test(ddd[,1],ddd[,4], paired = T, alternative="greater")

lowpositive__val <- c(mean(ddd[,1],na.rm = T),mean(ddd[,4],na.rm = T), tt$p.value)

#==Somatic=====

which <- c(15,14,17,20,21,22) - rep(1,6)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]] + dd[[4]] + dd[[5]] + dd[[6]]
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- wilcox.test(ddd[,1],ddd[,4], paired = T, alternative="greater")

somatic__val <- c(mean(ddd[,1],na.rm = T),mean(ddd[,4],na.rm = T), tt$p.value)

#===Interpersonal====

which <- c(7,11) - rep(1,2)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] 
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA

tt <- wilcox.test(ddd[,1],ddd[,4], paired = T, alternative="greater")
interpersonal__val <- c(mean(ddd[,1],na.rm = T),mean(ddd[,4],na.rm = T), tt$p.value)


print("Hence the values are")
print(rbind(cesd__val,depressed__val,lowpositive__val, somatic__val, interpersonal__val))


