library(DOS)
library(haven)
library(caret)
library(xtable)
library(lsr)
library(ggplot2)
library(sensitivitymw)
library(here)
#=================Data preparation============

data <- read_dta(here("Data", "WLS_master_data.dta"))
cath_match <- read.csv(here("Data", "Catholics_matched_pairs.csv"))
cath_match <- cath_match[,-c(1)]
data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: female



income_dat <- numeric()
for(k in 1: nrow(cath_match))
{
  hn <- data_f$z_rp001re[data_f$idpub==cath_match[k,2]]
  gn <- data_f$rp015ree[data_f$idpub==cath_match[k,2]]
  first_p <- c(hn,gn,cath_match[k,2])
  
  
  hn <- data_f$z_rp001re[data_f$idpub==cath_match[k,3]]
  gn <- data_f$rp015ree[data_f$idpub==cath_match[k,3]]
  second_p <- c(cath_match[k,2],gn,hn)
  
  income_dat <- rbind(income_dat,c(first_p,second_p))
}

#======Total Income boxplot======  
boxplot(income_dat[,2], income_dat[,5], main="total income")


#======Wages and salariesboxplot======  
boxplot(income_dat[,1], income_dat[,6], main="Wages and salaries")

