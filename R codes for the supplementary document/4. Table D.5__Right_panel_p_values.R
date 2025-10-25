library(sensitivitymw)
library(sensitivitymv)
library(DOS)
library(haven)
library(caret)
library(here)


#==This code reproduces all the p-values reported in the right panel of Table D.5.
#== The analysis plan prepared by team B (after they explored the non-Catholics data) is now being 
#== implemented on the Catholics data. 





#=We read four datasets for this exercise:

#== (i) Master data file (from the WLS) 

#== (ii) Matched pairs index for the Catholics 
#(we got this file after running our Risk-set matching algorithm on the master file).

#== (iii) Pregnancy intention information for the treated Catholics - extracted from the master file 
#== (iv) Pregnancy intention information for the control Catholics - extracted from the master file
 


data <- read_dta(here("Data", "WLS_master_data.dta"))
cath_match <- read.csv(here("Data", "Catholics_matched_pairs.csv"))
cath_match <- cath_match[,-c(1)]
trt_child <- read.csv(here("Data", "Pregnancy_Intention_info_for_treated_Catholics.csv"))
cont_child <- read.csv(here("Data", "Pregnancy_Intention_info_for_control_Catholics.csv"))


#== Extracting information for only the female graduates (as our study focuses on them only) 

data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: female





#============Step 1: test low-positive score===========
vv <- c("z_mu001rec","z_mu003rer","z_mu004rer","z_mu005rer","z_mu006rer","z_mu007rer",
        "z_mu008rer","z_mu009rer","z_mu010rer","z_mu011rer","z_mu012rer",
        "z_mu013rer","z_mu014rer","z_mu015rer","z_mu016rer","z_mu017rer",
        "z_mu018rer","z_mu019rer","z_mu020rer","z_mu021rer","z_mu022rer")
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

which <- c(16,6,9) - rep(1,3)
dd <- dat[which]
ddd <- dd[[1]]+dd[[2]] +dd[[3]]
ddd$X1 <- 21-ddd$X1
ddd$X4 <- 21-ddd$X4
ddd[,1][ddd[,1]<0]<- NA
ddd[,4][ddd[,4]<0]<- NA
ddd[,2]<- dd[[1]][,2] 
ddd[,3]<- dd[[1]][,3] 


senmwCI(ddd[,1]-ddd[,4], gamma = 1, one.sided = F)




#============Step 2: test low-positive score at gamma =1.2===========

senmwCI(ddd[,1]-ddd[,4], gamma = 1.2, one.sided = F)


#============ Step 3: Effect modification by age====

age <- numeric()
for(k in 1: nrow(ddd))
{
  a1 <- data_f$z_age75[data_f$idpub==ddd[k,2]]-(1975- cath_match[k,1])
  a2 <- data_f$z_age75[data_f$idpub==ddd[k,3]]-(1975- cath_match[k,4])
  age <- rbind(age,c(a1,a2))
}

X=(age[,1]+age[,2])/2
Xp3 = quantile(X,0.3)
Xp4  = quantile(X,0.4)
Xp6 = quantile(X,0.6)
temp1 <- (ddd[,1]-ddd[,4])[X<=Xp3]
temp2 <- (ddd[,1]-ddd[,4])[X>=Xp4 & X<=Xp6]
temp11 <- temp1[1:length(temp2)]

senmwCI(temp11-temp2, gamma = 1, one.sided = F)



#============ Step 4:CES-D depression score====

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




#============ Step 5:CES-D depression score for Gamma =1.2====

senmwCI(SP[,1]-SP[,4], gamma = 1.2, one.sided = F)


#============ Step 6: effect modification by age for Gamma =1.2====
senmwCI(temp11-temp2, gamma = 1.2, one.sided = F)
















