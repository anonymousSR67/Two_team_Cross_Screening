#==Library required======= 
library(haven)
library(caret)
library(optmatch)
library(ggplot2)
library(dplyr)
library(reshape2)
library(designmatch)
library(survival)
library(devtools)
library(here)



data <- read_dta(here("Data", "WLS_master_data.dta"))
data_f <- data.frame(data)
data_f <- data_f[which(data_f$rtype=="g"),] #filter: graduate
data_f <- data_f[which(data_f$z_sexrsp==2),] #filter: females


idp <- which(colnames(data_f)=="idpub")

#===========Pregnancy==========
a <- which(colnames(data_f)=="wantb1") #first birth
b <- which(colnames(data_f)=="wantb2") #second birth
c <- which(colnames(data_f)=="wantb3") #third birth
last <- which(colnames(data_f)=="wantbl") #last birth

#==============First depression age=====
nnn <- which(colnames(data_f)=="z_ru020re") #age at the time of first ever depression

#=========Childhood Measure====
d <- which(colnames(data_f)=="hsrankq") #School Rank
f <- which(colnames(data_f)=="gwiiq_bm") #IQ
dd <- which(colnames(data_f)=="ses57") #parental status

#===== Adulthood Measure=======
g <- which(colnames(data_f)=="z_edeqyr") #yrs of education
i <- which(colnames(data_f)=="pop57") # Town size


#=========== Personality Measures==================
o <- which(colnames(data_f)=="z_mh001rec") #extro
q <- which(colnames(data_f)=="z_mh009rec") #Agree
s <- which(colnames(data_f)=="z_mh017rec") #Cons
u <- which(colnames(data_f)=="z_mh025rec") #neuro
w <- which(colnames(data_f)=="z_mh032rec") #open


#================Religion===================
rlg <- which(colnames(data_f)=="relfml")



birth_year <- 1975 - data_f$z_age75
kids_no <- data_f$z_kidsno
first_birth_year <- data_f$z_cmkdb1
second_birth_year <- data_f$z_cmkdb2
third_birth_year <- data_f$z_cmkdb3
fourth_birth_year <- data_f$z_cmkdb4
fifth_birth_year <- data_f$z_cmkdb5
sixth_birth_year <- data_f$z_cmkdb6
seventh_birth_year <- data_f$z_cmkdb7
eighth_birth_year <- data_f$z_cmkdb8
ninth_birth_year <- data_f$z_cmkdb9
tenth_birth_year <- data_f$z_cmkdb0




data_try <- data_f[-which(data_f$relfml %in% c(1,-3)),]

#############====Missing Replacement===##################

data_try$z_ru020remis=is.na(data_try$z_ru020re) #Note: -2 is not a missing value
ind2 <- which(data_try$z_ru020re>0)
data_try$z_ru020re[data_try$z_ru020remis==1]=mean(data_try$z_ru020re[ind2])

data_try$hsrankq[data_try$hsrankq==-3]=NA
data_try$hsrankqmis=is.na(data_try$hsrankq)
data_try$hsrankq[data_try$hsrankqmis==1]=mean(data_try$hsrankq,na.rm=TRUE)


data_try$z_mh001rec[data_try$z_mh001rec==-3]=NA
data_try$z_mh001recmis=is.na(data_try$z_mh001rec)
data_try$z_mh001rec[data_try$z_mh001recmis==1]=mean(data_try$z_mh001rec,na.rm=TRUE)

data_try$z_mh009rec[data_try$z_mh009rec==-3]=NA
data_try$z_mh009recmis=is.na(data_try$z_mh009rec)
data_try$z_mh009rec[data_try$z_mh009recmis==1]=mean(data_try$z_mh009rec,na.rm=TRUE)

data_try$z_mh017rec[data_try$z_mh017rec==-3]=NA
data_try$z_mh017recmis=is.na(data_try$z_mh017rec)
data_try$z_mh017rec[data_try$z_mh017recmis==1]=mean(data_try$z_mh017rec,na.rm=TRUE)


data_try$z_mh025rec[data_try$z_mh025rec==-3]=NA
data_try$z_mh025recmis=is.na(data_try$z_mh025rec)
data_try$z_mh025rec[data_try$z_mh025recmis==1]=mean(data_try$z_mh025rec,na.rm=TRUE)

data_try$z_mh032rec[data_try$z_mh032rec==-3]=NA
data_try$z_mh032recmis=is.na(data_try$z_mh032rec)
data_try$z_mh032rec[data_try$z_mh032recmis==1]=mean(data_try$z_mh032rec,na.rm=TRUE)


#==Functions for rank based Mahalanobis Distance and adding caliper

smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    un <- numeric()
    for(j in 1:k)
    {
      un[j] <- length(unique(X[,j]))
    }
    if(length(which(un==1))>0)
    {
      X <- X[,-which(un==1)]
    }
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    icov<-ginv(cv)
    if(nrow(out)==1)
    {
      for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt,icov,inverted=T)
    }else{
      for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)  
    }
    out
  }

addcaliper=function(dmat,z,logitp,calipersd=.5,penalty=1000){
  # Pooled within group standard devation
  if(length(which(z==1))>1)
  {
    sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
  }else{
    sd.logitp <- sqrt((0+sd(logitp[z==0])^2)/2)
  }
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}



#===This section prepares the data for time-dependent Cox PH model===


data_check <- data.frame(id = data_f$idpub, age_at_75 = data_f$z_age75, birth_year, kids_no, first_birth_year, second_birth_year, third_birth_year,
                         fourth_birth_year, fifth_birth_year, sixth_birth_year, seventh_birth_year, eighth_birth_year, ninth_birth_year,
                         tenth_birth_year)


data_check <- data_check[-which(data_f$relfml %in% c(1,-3)),]


for(t in 1: nrow(data_check))
{
  if(data_check$kids_no[t]<10 & data_check$kids_no[t]>0)
  {
    data_check[t,(4+data_check$kids_no[t]+1):(ncol(data_check))] <- NA
  }
}

data_time <- data_check[,-c(2:4)]
melt_data <- melt(data_time, id.vars ="id")
melt_data_sorted <- melt_data[order(melt_data$id),]

preg_int <- data_try[,c(a,b,c)]
preg_int <- cbind(preg_int, matrix(0, nrow = nrow(data_check), ncol = 7))

for(p in 1: nrow(preg_int))
{
  
  if(data_try$z_kidsno[p] > 3 & data_try$z_kidsno[p] < 11)
  {
    preg_int[p,data_try$z_kidsno[p]] <- data_try$wantbl[p]
  }
}

preg_int <- as.matrix(preg_int)
preg_int_vec <- matrixcalc:: vec(t(preg_int))

melt_data_sorted$preg_int <- preg_int_vec



id_with_child_info <- data_check$id[-which(data_check$kids_no <0)]
melt_data_sorted <- melt_data_sorted[melt_data_sorted$id %in% id_with_child_info,]


melt_data_store <- melt_data_sorted[-which(is.na(melt_data_sorted$value)==T),]
melt_data_store$preg_int[which(melt_data_store$preg_int==0)] <- 1


zero_birth_data <- melt_data_store[melt_data_store$value<0,]
aa <- unique(zero_birth_data$id)
bb <- data_check$id[which(data_check$kids_no==0)]
zero_birth_data <- zero_birth_data[which(zero_birth_data$id %in% aa[which(aa %in% bb)]),]


melt_data_store <- melt_data_store[-which(melt_data_store$id %in% unique(zero_birth_data$id)),]


melt_data_store <- melt_data_store[which(melt_data_store$preg_int >0),]
melt_data_store <- melt_data_store[which(melt_data_store$value>0),]


no_birth_yet <- data.frame(id=unique(zero_birth_data$id), variable="zero_birth", value=0, preg_int=0)
melt_data_store <- rbind(melt_data_store,no_birth_yet)


child_no <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  child_no <- c(child_no,c(0:(pp-1)))
}



age_preg <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  age_75 <- data_f$z_age75[data_f$idpub==unique(melt_data_store$id)[j]]
  age_temp <- age_75 - (1975 - melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],3])
  age_preg <- c(age_preg, age_temp)
}


eq_ed_year <- numeric()
for(j in 1: length(unique(melt_data_store$id)))
{
  eq_ed_temp <- data_f$z_edeqyr[data_f$idpub==unique(melt_data_store$id)[j]]
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  eq_ed_year <- c(eq_ed_year,rep(eq_ed_temp,pp))
}


temp <- data.frame(t1=age_preg-6, t2= eq_ed_year)
years_of_education <- apply(temp,1,min)



first_dep <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  fi_dep_temp <- data_try$z_ru020re [data_try$idpub==unique(melt_data_store$id)[j]]
  first_dep <- c(first_dep,rep(fi_dep_temp,pp))
}


first_dep_mis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  fi_dep_mis_temp <- data_try$z_ru020remis [data_try$idpub==unique(melt_data_store$id)[j]]
  first_dep_mis <- c(first_dep_mis,rep(fi_dep_mis_temp,pp))
}






melt_data_store$child_no <- child_no
melt_data_store$age_preg <- age_preg
melt_data_store$years_of_edu <- years_of_education
melt_data_store$first_dep_age <- first_dep
melt_data_store$first_dep_age_mis <- first_dep_mis

first_dep_age_t_2 <- melt_data_store$first_dep_age < (melt_data_store$age_preg-2)
first_dep_age_t_2 <- as.numeric(first_dep_age_t_2)
first_dep_age_t_2[melt_data_store$first_dep_age==-2] <- 0
first_dep_age_t_2[melt_data_store$first_dep_age==-1] <- 0
melt_data_store$first_dep_age_t_2 <- first_dep_age_t_2


mar_status <- numeric()
mar_data1 <- which(colnames(data_f)=="z_cmfmbg")
mar_data2 <- which(colnames(data_f)=="z_cmenfm")
mar_data3 <- which(colnames(data_f)=="z_cmsmbg")
mar_data4 <- which(colnames(data_f)=="z_cmensm")
mar_data5 <- which(colnames(data_f)=="z_cmtmbg")
mar_data6 <- which(colnames(data_f)=="z_cmentm")

for(j in 1:length(unique(melt_data_store$id)))
{
  which_id <- unique(melt_data_store$id)[j]
  
  mar_hist <- data_f[data_f$idpub==which_id,c(mar_data1,mar_data2, mar_data3, mar_data4, mar_data5, mar_data6)]
  mar_hist[mar_hist==-2] <- 10000000
  id_match_year <- melt_data_store$value[melt_data_store$id==unique(melt_data_store$id)[j]]
  id_match_year <- data.frame(yr=(id_match_year-1900)*12+1)
  check_func <- function(cc,vv=mar_hist)
  {
    a1 <- cc>=vv[1] & cc<=vv[2]
    a2 <- cc>=vv[3] & cc<=vv[4]
    a3 <- cc>=vv[5] & cc<=vv[6]
    fin <- max(a1,a2,a3)
    return(fin)
  }
  
  mar_status <- c(mar_status,apply(id_match_year,1,check_func))
}

melt_data_store$marital_sta <- mar_status


rank <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  rank_temp <- data_try$hsrankq [data_try$idpub==unique(melt_data_store$id)[j]]
  rank <- c(rank,rep(rank_temp,pp))
}


rank_mis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  rank_mis_temp <- data_try$hsrankqmis [data_try$idpub==unique(melt_data_store$id)[j]]
  rank_mis <- c(rank_mis,rep(rank_mis_temp,pp))
}

melt_data_store$rank <- rank
melt_data_store$rankmis <- rank_mis

#==== Create "IQ"===
IQ <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  IQ_temp <- data_try$gwiiq_bm [data_try$idpub==unique(melt_data_store$id)[j]]
  IQ <- c(IQ,rep(IQ_temp,pp))
}
melt_data_store$IQ <- IQ

#==== Create "ses57"===
ses57 <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  ses57_temp <- data_try$ses57 [data_try$idpub==unique(melt_data_store$id)[j]]
  ses57 <- c(ses57,rep(ses57_temp,pp))
}
melt_data_store$ses57 <- ses57

#==== Create "pop57"===
pop57 <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  pop57_temp <- data_try$pop57 [data_try$idpub==unique(melt_data_store$id)[j]]
  pop57 <- c(pop57,rep(pop57_temp,pp))
}
melt_data_store$pop57 <- pop57


#==== Create "z_mh001rec"===
z_mh001rec <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh001rec_temp <- data_try$z_mh001rec [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh001rec <- c(z_mh001rec,rep(z_mh001rec_temp,pp))
}
melt_data_store$z_mh001rec <- z_mh001rec

z_mh001recmis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh001recmis_temp <- data_try$z_mh001recmis [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh001recmis <- c(z_mh001recmis,rep(z_mh001recmis_temp,pp))
}
melt_data_store$z_mh001recmis <- z_mh001recmis



#==== Create "z_mh009rec"===
z_mh009rec <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh009rec_temp <- data_try$z_mh009rec [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh009rec <- c(z_mh009rec,rep(z_mh009rec_temp,pp))
}
melt_data_store$z_mh009rec <- z_mh009rec

z_mh009recmis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh009recmis_temp <- data_try$z_mh009recmis [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh009recmis <- c(z_mh009recmis,rep(z_mh009recmis_temp,pp))
}
melt_data_store$z_mh009recmis <- z_mh009recmis



#==== Create "z_mh017rec"===
z_mh017rec <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh017rec_temp <- data_try$z_mh017rec [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh017rec <- c(z_mh017rec,rep(z_mh017rec_temp,pp))
}
melt_data_store$z_mh017rec <- z_mh017rec

z_mh017recmis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh017recmis_temp <- data_try$z_mh017recmis [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh017recmis <- c(z_mh017recmis,rep(z_mh017recmis_temp,pp))
}
melt_data_store$z_mh017recmis <- z_mh017recmis


#==== Create "z_mh025rec"===
z_mh025rec <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh025rec_temp <- data_try$z_mh025rec [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh025rec <- c(z_mh025rec,rep(z_mh025rec_temp,pp))
}
melt_data_store$z_mh025rec <- z_mh025rec

z_mh025recmis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh025recmis_temp <- data_try$z_mh025recmis [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh025recmis <- c(z_mh025recmis,rep(z_mh025recmis_temp,pp))
}
melt_data_store$z_mh025recmis <- z_mh025recmis


#==== Create "z_mh032rec"===
z_mh032rec <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh032rec_temp <- data_try$z_mh032rec [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh032rec <- c(z_mh032rec,rep(z_mh032rec_temp,pp))
}
melt_data_store$z_mh032rec <- z_mh032rec

z_mh032recmis <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  z_mh032recmis_temp <- data_try$z_mh032recmis [data_try$idpub==unique(melt_data_store$id)[j]]
  z_mh032recmis <- c(z_mh032recmis,rep(z_mh032recmis_temp,pp))
}
melt_data_store$z_mh032recmis <- z_mh032recmis


#======== Calculate tstart and tstop=======
tstart <- numeric()
for(j in 1:length(unique(melt_data_store$id)))
{
  pp <- nrow(melt_data_store[melt_data_store$id==unique(melt_data_store$id)[j],])
  temp <- melt_data_store$value[melt_data_store$id==unique(melt_data_store$id)[j]]
  tempp <- c(1957,temp[-c(pp)])
  tstart <- c(tstart,tempp)
}
melt_data_store$tstart <- tstart
melt_data_store$tstop <- melt_data_store$value

#====== treated and control======
melt_data_store$treated <- 0
melt_data_store$treated[melt_data_store$preg_int==3] <- 1 

#=======tstop for "no_birth"=======
melt_data_store$tstop[melt_data_store$preg_int==0] <- 1975

#===== years of education for "no_birth"========

id_collect <- melt_data_store$id[melt_data_store$preg_int==0]
for(j in 1: length(id_collect))
{
  melt_data_store$years_of_edu[melt_data_store$id==id_collect[j]] <- data_try$z_edeqyr[data_try$idpub==id_collect[j]]
}

mm <- melt_data_store[which(melt_data_store$tstart<melt_data_store$tstop),]


#===== Hazard Coefficients=======

propscore.model <- coxph( Surv(tstart, tstop, treated) ~ child_no+years_of_edu+first_dep_age_t_2
                          +first_dep_age_mis+rank+rankmis+IQ+ses57+pop57
                          +z_mh001rec+z_mh001recmis
                          +z_mh009rec+z_mh009recmis
                          +z_mh017rec+z_mh025rec+z_mh025recmis
                          +z_mh032rec+z_mh032recmis+age_preg+marital_sta, mm)

hazard_coef <- propscore.model$coefficients
hazard_coef[1] <- hazard_coef[1]*10
hazard_coef[19] <- hazard_coef[19]*10


matched_pair <- numeric()
trt_var <- numeric()
cont_var <- numeric()

hold_data <- numeric()
loveplot <- list()
counter <- 1

repeat {
  unintent_data <- melt_data_store[melt_data_store$preg_int==3,]
  temp_year <- sort(unique(unintent_data$value))[1]
  unintent_index <- unintent_data[unintent_data$value==temp_year,]
  cum_data <- melt_data_store[melt_data_store$value<(temp_year+1),]
  max_data <- cum_data %>% group_by(id) %>% 
    summarise_at(vars(preg_int), max)
  int_index <- max_data[max_data$preg_int<3,]
  sub_data_trt <- melt_data_store[which(melt_data_store$id %in% unintent_index$id),]
  sub_data_trt <- sub_data_trt[sub_data_trt$value==temp_year,]
  sub_data_trt$treated <- 1
  sub_data_cont <- melt_data_store[which(melt_data_store$id %in% int_index$id),]
  sub_data_cont <- sub_data_cont[sub_data_cont$value < (temp_year+1),]
  sub_data_cont$treated <- 0
  
  
  temp2 <- numeric()
  for(j in 1:length(unique(sub_data_cont$id)))
  {
    temp <- length(which(sub_data_cont$id==unique(sub_data_cont$id)[j]))
    temp2 <- c(temp2,temp) 
  }
  temp2 <- cumsum(temp2)
  sub_data_cont <- sub_data_cont[temp2,]
  
  temp1 <- numeric()
  for(j in 1:length(unique(sub_data_trt$id)))
  {
    tt <- sub_data_trt[sub_data_trt$id==unique(sub_data_trt$id)[j],]
    temp1 <- rbind(temp1, tt[1,])
  }
  sub_data_trt <- temp1
  
  sub_data <- rbind(sub_data_trt, sub_data_cont)
  
  
  
  hold_data <- rbind(hold_data,sub_data)
  
  #==============Capture the covariates ======
  
  temp_model=glm(treated ~child_no+years_of_edu+first_dep_age_t_2
                 +first_dep_age_mis+rank+rankmis+IQ+ses57+pop57
                 +z_mh001rec+z_mh001recmis
                 +z_mh009rec+z_mh009recmis
                 +z_mh017rec + z_mh017recmis
                 +z_mh025rec+z_mh025recmis
                 +z_mh032rec+z_mh032recmis+age_preg+marital_sta
                 ,family=binomial,x=TRUE,y=TRUE,data=sub_data)
  
  
  dmy=dummyVars(temp_model$formula,data=sub_data)
  Xmat=data.frame(predict(dmy,newdata=sub_data))  
  Xmatmahal=Xmat
  
  treated=sub_data$treated
  hazard <- numeric()
  for(j in 1: nrow(Xmatmahal))
  {
    temp <- sum(hazard_coef * Xmatmahal[j,-c(15)])
    hazard <- c(hazard,temp)
  }
  
  
  #sub_data$logit.ps=predict(temp_model)
  sub_data$logit.ps=hazard
  sub_data$treatment=sub_data$treated
  
  
  
  # Make the rownames in datatemp be 1:number of rows
  rownames(sub_data)=seq(1,nrow(sub_data),1) 
  
  
  
  
  if(nrow(unintent_index)==1)
  {
    distmat=smahal(sub_data$treated,Xmatmahal)
  }else{
    distmat = distmat(sub_data$treated,Xmatmahal)
  }
  # Add caliper
  distmat=addcaliper(distmat,sub_data$treated,sub_data$logit.ps,calipersd=.5)
  
  # Label the rows and columns of the distance matrix by the rownames in sub_data
  rownames(distmat)=rownames(sub_data)[sub_data$treated==1]
  colnames(distmat)=rownames(sub_data)[sub_data$treated==0]
  
  # Optimal matching to multiple controls
  
  # Number of controls to match to each treated unit
  nocontrols.per.match=1
  matchvec=pairmatch(distmat,controls=nocontrols.per.match,data=sub_data)
  sub_data$matchvec=matchvec
  
  ## Create a matrix saying which control units each treated unit is matched to
  ## Create vectors of the subject indices of the treatment units ordered by
  ## their matched set and corresponding control unit
  treated.subject.index=rep(0,sum(sub_data$treated==1))
  matched.control.subject.index.mat=matrix(rep(0,nocontrols.per.match*length(treated.subject.index)),ncol=nocontrols.per.match)
  matchedset.index=substr(matchvec,start=3,stop=10)
  matchedset.index.numeric=as.numeric(matchedset.index)
  for(i in 1:length(treated.subject.index)){
    matched.set.temp=which(matchedset.index.numeric==i)
    treated.temp.index=which(sub_data$treated[matched.set.temp]==1)
    treated.subject.index[i]=matched.set.temp[treated.temp.index]
    matched.control.subject.index.mat[i,]=matched.set.temp[-treated.temp.index]
  }
  matched.control.subject.index=matched.control.subject.index.mat
  
  id_trt <- sub_data$id[treated.subject.index]
  id_cont <- cbind(sub_data$id[matched.control.subject.index[,1]])
  
  temp_int_col <- numeric()
  for(pp in 1: length(id_cont[,1]))
  {
    temp_int_col <- c(temp_int_col,sub_data_cont$preg_int[sub_data_cont$id==id_cont[,1][pp]])
  }
  
  y1 <- numeric()
  y2 <- numeric()
  for(kk in 1: length(id_trt))
  {
    y1 <- c(y1, sub_data$value[sub_data$id==id_trt[kk]])
    y2 <- c(y2, sub_data$value[sub_data$id==id_cont[,1][kk]])
  }
  res_temp <- rbind(cbind(y1,id_trt, id_cont[,1],y2))
  
  matched_pair <- rbind(matched_pair,res_temp)
  t1 <- numeric()
  for(tt in 1: length(id_trt))
  {
    t1 <- c(t1,which(sub_data$id==id_trt[tt]))
  }
  
  t2 <- numeric()
  for(tt in 1: length(id_cont))
  {
    t2 <- c(t2,which(sub_data$id==id_cont[tt]))
  }
  t1_temp <- Xmatmahal[t1,]
  trt_var <- rbind (trt_var,t1_temp)
  t2_temp <- Xmatmahal[t2,]
  cont_var <- rbind(cont_var, t2_temp)
  
  exclude_ind <- as.numeric(matrixcalc:: vec(res_temp))
  #exclude_ind <- res_temp
  melt_data_store <- melt_data_store[-which(melt_data_store$id %in% exclude_ind),]
  
  if (sum(melt_data_store$preg_int==3)==0){
    break
  }
}
trt_var1 <- trt_var[-which(abs(trt_var$age_preg- cont_var$age_preg)>4),]
cont_var1 <- cont_var[-which(abs(trt_var$age_preg- cont_var$age_preg)>4),]
trt_mean <- apply(trt_var1, 2, mean)
cont_mean <- apply(cont_var1, 2, mean)
trt_variance <- apply(trt_var1,2,var)
cont_variance <- apply(cont_var1,2,var)
stand_diff_after <- (trt_mean - cont_mean)/sqrt((trt_variance+cont_variance)/2)
matched_pair1 <- matched_pair[-which(abs(trt_var$age_preg- cont_var$age_preg)>4),]
trt_var2 <- trt_var1[-which(abs(trt_var1$child_no- cont_var1$child_no)>1),]
cont_var2 <- cont_var1[-which(abs(trt_var1$child_no- cont_var1$child_no)>1),]
trt_mean <- apply(trt_var2, 2, mean)
cont_mean <- apply(cont_var2, 2, mean)
trt_variance <- apply(trt_var2,2,var)
cont_variance <- apply(cont_var2,2,var)
stand_diff_after <- (trt_mean - cont_mean)/sqrt((trt_variance+cont_variance)/2)
matched_pair2 <- matched_pair1[-which(abs(trt_var1$child_no- cont_var1$child_no)>1),]
inde <- which(abs(trt_var2$first_dep_age_t_2- cont_var2$first_dep_age_t_2)>0)
indee <- inde[-c(1,2)]
trt_var3 <- trt_var2[-indee,]
cont_var3 <- cont_var2[-indee,]
trt_mean <- apply(trt_var3, 2, mean)
cont_mean <- apply(cont_var3, 2, mean)
trt_variance <- apply(trt_var3,2,var)
cont_variance <- apply(cont_var3,2,var)
stand_diff_after <- (trt_mean - cont_mean)/sqrt((trt_variance+cont_variance)/2)
matched_pair3 <- matched_pair2[-indee,]


before_model=glm(treated ~child_no+years_of_edu+first_dep_age_t_2
                 +first_dep_age_mis+rank+rankmis+IQ+ses57+pop57
                 +z_mh001rec+z_mh001recmis
                 +z_mh009rec+z_mh009recmis
                 +z_mh017rec + z_mh017recmis
                 +z_mh025rec+z_mh025recmis
                 +z_mh032rec+z_mh032recmis+age_preg+marital_sta
                 ,family=binomial,x=TRUE,y=TRUE,data=hold_data)


dmy=dummyVars(before_model$formula,data=hold_data)
Xmat_before=data.frame(predict(dmy,newdata=hold_data))  
trt_var_before <- Xmat_before[which(hold_data$treated==1),]
cont_var_before <- Xmat_before[which(hold_data$treated==0),]
trt_mean_before <- apply(trt_var_before, 2, mean)
cont_mean_before <- apply(cont_var_before, 2, mean)
trt_variance_before <- apply(trt_var_before,2,var)
cont_variance_before <- apply(cont_var_before,2,var)
stand_diff_before <- (trt_mean_before - cont_mean_before)/sqrt((trt_variance_before+cont_variance_before)/2)



abs.stand.diff.before=abs(stand_diff_before)
abs.stand.diff.after=abs(stand_diff_after)
covariates=names(stand_diff_before)
covariates <- c("no_children","education","depressed_before","depressed_before_mis",
                "high_school_rank","high_school_rank_mis","IQ", "parents_socioec","town_pop_child",
                "extraversion","extraversion_mis","agreeableness","agreeableness_mis",
                "conscientious", "conscientious_mis","neuroticism","neuroticism_mis",
                "openness","openness_mis","age","married")
plot.dataframe=data.frame(abs.stand.diff=c(abs.stand.diff.before,abs.stand.diff.after),covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates))))
ggplot(plot.dataframe,aes(x=abs.stand.diff,y=covariates))+geom_point(size=5,aes(shape=factor(type)))+scale_shape_manual(values=c(4,1))+geom_vline(xintercept=c(.1,.2),lty=2)



