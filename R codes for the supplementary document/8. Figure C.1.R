
library(ggplot2)
library("gridExtra") 
library(dplyr)
library(latex2exp)

##################Creating the data frame###########
columns= c("Mu","Numrep","Power", "FWER", "Method", "Metric") 
muvec=seq(1, 5, by=0.1)
lenmuvec=length(muvec)
numrepvec=c(1,3,6,10,13,16)
lennumrepvec=length(c(1,3,6,10,13,16))
nrowtotal=lenmuvec*lennumrepvec*7
MyData=data.frame(matrix(nrow=nrowtotal, ncol=length(columns)))
colnames(MyData)=columns
lenmethod=lenmuvec*lennumrepvec
MyData$Method=c(rep("Fisher", lenmethod), rep("CombinedHolmGlobal", lenmethod), rep("Cross_glob", lenmethod), rep("HolmAllData", lenmethod), rep("CombinedHolmRep", lenmethod), rep("Holmmax", lenmethod), rep("Cross_rep", lenmethod))
MyData$Mu=rep(muvec, lennumrepvec*7)
MyData$Numrep=rep(c(rep(1, lenmuvec), rep(3, lenmuvec), rep(6, lenmuvec), rep(10, lenmuvec), rep(13, lenmuvec), rep(16, lenmuvec)),7) 
MyData$Metric= c(rep("Globalnull",  lenmethod*4), rep("Replicability", lenmethod*3))

#----- Here is the main function--------------------------------
Globalfunc<-function(N, muvec, Numrepvec, MyData, alpha){
  mulen=length(muvec)
  I1=200
  I2=200
  m=16
    for(r in 1:length(numrepvec)){
      Numrepval=numrepvec[r]
      print(Numrepval)
    #  Intializing all the estimator matrices
      fwerrepmax=matrix(NA, nrow=mulen, ncol=N)
      powerrepmax=matrix(NA, nrow=mulen, ncol=N)
      fwerfisher=matrix(NA, nrow=mulen, ncol=N)
      powerfisher=matrix(NA, nrow=mulen, ncol=N)
      fwerholm32global=matrix(NA, nrow=mulen, ncol=N)
      powerholm32global=matrix(NA, nrow=mulen, ncol=N)
      fwerholm32rep=matrix(NA, nrow=mulen, ncol=N)
      powerholm32rep=matrix(NA, nrow=mulen, ncol=N)
      fwerholmalldata=matrix(NA, nrow=mulen, ncol=N)
      powerholmalldata=matrix(NA, nrow=mulen, ncol=N)
      fwerglobalcrossscreen=matrix(NA, nrow=mulen, ncol=N)
      powerglobalcrossscreen=matrix(NA, nrow=mulen, ncol=N)
      fwerrepcrossscreen=matrix(NA, nrow=mulen, ncol=N)
      powerrepcrossscreen=matrix(NA, nrow=mulen, ncol=N)
      ###------------------------------------------
      for (i in 1:N){   
        #We run N iterations 
    datamat1=matrix(NA, nrow=I1, ncol=m)#datamat 1: Column j will contain the I1 differences from study 1 for outcome j
    datamat2=matrix(NA, nrow=I2, ncol=m)#datamat 2: Column j will contain the I2 differences from study 2 for outcome j
    y1=matrix(NA, nrow=I1, ncol=m)#These will be the base null differences in study 1 for iteration i
    y2=matrix(NA, nrow=I1, ncol=m)#These will be the base null differences in study 2 for iteration i
    
    for (j in 1:m){#We go over all the m=16 outcomes
      y1[ ,j]=rnorm(I1, 0, 1)#These are the base differences for study 1 under H_0 for outcome j, in iteration i
      y2[ ,j]=rnorm(I2, 0, 1)#These are the base differences for study 2 under H_0 for outcome j, in iteration i
    }
    #Now we will go over all the values of mu and for each value, will obtain the observation ("difference") for the false null hypothesis by 
    #adding this value of mu to the corresponding base observation
    for (u in 1:mulen){
      muval=muvec[u]
      resdata=Datagenfunc(Numrepval, muval)#####Generates the sets of indices of I_11, I_00, I_01, I_10, and the vector of mu values for the 16 hypotheses for the given value of Numrep (number of replicated signals)
      I11ind=resdata$I11ind
      I00ind=resdata$I00ind
      I10ind=resdata$I10ind
      I01ind=resdata$I01ind
      muvec1=resdata$muvec1
      muvec2=resdata$muvec2
      #Now we get the differences for study 1,2 for the given value of mu for each outcome j
      for (j in 1:m){
      datamat1[,j]=y1[ ,j]+muvec1[j]/sqrt(I1)
      datamat2[,j]=y2[ ,j]+muvec2[j]/sqrt(I2)
      }
    pvalres1=Computepval(datamat1, datamat2, m)
    pvalvec1left=pvalres1$pvalvec1left #left-sided p-values for study 1
    pvalvec2left=pvalres1$pvalvec2left #left-sided p-values for study 2
    pvalvec1right=pvalres1$pvalvec1right #right-sided p-values for study 1
    pvalvec2right=pvalres1$pvalvec2right #right-sided p-values for study 2
    
    ####Pvalues for Combined Holm (on 32 hypotheses): Two-sided p-values for all the 32 hypotheses
    pvalvec1ts=pvalres1$pvalvec1twosided
    pvalvec2ts=pvalres1$pvalvec2twosided
    
    ####Pvalues for HolmAllData (P-values based on data from both studies for each of the 16 hypotheses)
    pvalvecalldata=pvalres1$pvalalldata
    
    ####Computing the p-values for Holmmax (Holm on max pvalues)
    pvalmaxleft=pmax(pvalvec1left, pvalvec2left)
    pvalmaxright=pmax(pvalvec1right, pvalvec2right)
    pvalmax=2*pmin(pvalmaxleft, pvalmaxright)
  
    #---Holm on Max: Competitor for replicability------
    
    indrejmax=which(p.adjust(pvalmax, method="holm")<alpha)
    fwerrepmax[u,i]=1*(length(intersect(c(I00ind,I10ind, I01ind), indrejmax))>0)
    powerrepmax[u,i]=length(intersect(I11ind, indrejmax))/length(I11ind)
    ######################################################################
    #--- Holm on Fisher's combined p-values: Competitor for Global null----
    #######################################################################
    
    pvalfisherleft=1-pchisq(-2*((log(pvalvec1left)+log(pvalvec2left))), 4)
    pvalfisherright=1-pchisq(-2*((log(pvalvec1right)+log(pvalvec2right))), 4)
    pvalfisher=2*pmin(pvalfisherleft, pvalfisherright)#Fisher's p-values
    indrejfisher=which(p.adjust(pvalfisher, method="holm")<alpha)#Applying Holm at level alpha on Fisher's p-values
    fwerfisher[u,i]=1*(length(intersect(I00ind, indrejfisher))>0)#Global null FWER
    powerfisher[u,i]=length(intersect(c(I01ind, I10ind, I11ind), indrejfisher))/length(c(I01ind, I10ind, I11ind))#Global null power
    #####--Holm on 32 hypotheses for global nulls, replicability---
    
    indrejholm32=which(p.adjust(c(pvalvec1ts,pvalvec2ts), method="holm")<alpha)#Applying Holm at level alpha on 32 two-sided p-values (16 from study 1 and 16 from study 2)
    #Computing global null and replicability discoveries of this procedure
    indvec=rep(0,2*m)
    indvec[indrejholm32]=1
    indrejholmglobalind=rep(0, m)
    for (k in 1:m){
      indrejholmglobalind[k]=(indvec[k]>0)+(indvec[k+m]>0)
      }
    indrejholmglobal=which(indrejholmglobalind>0)#Global null discovery is made if the corresponding hypothesis is rejected in at least one of the studies
    indrejholmrep=which(indrejholmglobalind==2)#Replicability discovery is made if the corresponding hypothesis is rejected in both studies
    fwerholm32global[u,i]=1*(length(intersect(I00ind, indrejholmglobal))>0)
    powerholm32global[u,i]=length(intersect(c(I01ind, I10ind, I11ind), indrejholmglobal))/length(c(I01ind, I10ind, I11ind))
    fwerholm32rep[u,i]=1*(length(intersect(c(I00ind,I10ind, I01ind), indrejholmrep))>0)
    powerholm32rep[u,i]=length(intersect(I11ind, indrejholmrep))/length(I11ind)
    
    ####--Holm on 16 hypotheses with combined data from the two studies-------
    
    indrejholmalldata=which(p.adjust(pvalvecalldata, method="holm")<alpha)
    fwerholmalldata[u,i]=1*(length(intersect(I00ind, indrejholmalldata))>0)
    powerholmalldata[u,i]=length(intersect(c(I01ind, I10ind, I11ind), indrejholmalldata))/length(c(I01ind, I10ind, I11ind))
    
    ####--Cross-screening--------------------------------
    
    
    pvalvec1cross=rep(0, m)
    pvalvec2cross=rep(0, m)
    for (k in 1:m){
      if (pvalvec1right[k]<pvalvec1left[k]){pvalvec2cross[k]=pvalvec2right[k]}
      if (pvalvec1right[k]>=pvalvec1left[k]){pvalvec2cross[k]=pvalvec2left[k]}
      if (pvalvec2right[k]<pvalvec2left[k]){pvalvec1cross[k]=pvalvec1right[k]}
      if (pvalvec2right[k]>=pvalvec2left[k]){pvalvec1cross[k]=pvalvec1left[k]}
    }
    indselstudy1=which(pmin(pvalvec2right, pvalvec2left)<alpha/2)#Indices that are selected for study 1
    pvalstudy1sel=pvalvec1cross[indselstudy1]#P-values that are selected in study 1
    indrejcrosstemp=which(p.adjust(pvalstudy1sel, method="holm")<alpha/2)#Apply Holm's procedure at level alpha/2 on the selected study 1 p-values - these are indices among the selected p-values
    indrejcrossstudy1=indselstudy1[indrejcrosstemp]#Returning to the original indices
    
    indselstudy2=which(pmin(pvalvec1right, pvalvec1left)<alpha/2)#Indices selected for study 2
    pvalstudy2sel=pvalvec2cross[indselstudy2]#These are the selected p-values in study 2
    indrejcrosstemp=which(p.adjust(pvalstudy2sel, method="holm")<alpha/2)#Apply Holm's procedure at level alpha/2 on the selected study 2 p-values - these are indices among the selected p-values
    indrejcrossstudy2=indselstudy2[indrejcrosstemp]#Returning to the original indices
    
    indrep=intersect(indrejcrossstudy1, indrejcrossstudy2)#Replicability discoveries according to cross-screening
    indglobal=c(indrejcrossstudy1, indrejcrossstudy2)#Global null discoveries according to cross-screening
    #Calculating the realized FWER and power for global nulls and replicability according to cross screening
    fwerglobalcrossscreen[u,i]=1*(length(intersect(I00ind, indglobal))>0)
    powerglobalcrossscreen[u,i]=length(intersect(c(I01ind, I10ind, I11ind), indglobal))/length(c(I01ind, I10ind, I11ind))
    fwerrepcrossscreen[u,i]=1*(length(intersect(c(I00ind,I10ind, I01ind), indrep))>0)
    powerrepcrossscreen[u,i]=length(intersect(I11ind, indrep))/length(I11ind)
  }#End of iterations running on the values of mu
}#End of all the iterations for (i in 1:N)
    
   
    # Putting the results into the dataframe 
    ###Holm on Fisher - Competitor for global null
    MyData$FWER[MyData$Method=="Fisher"&MyData$Numrep==Numrepval]=rowMeans(fwerfisher)
    MyData$Power[MyData$Method=="Fisher"&MyData$Numrep==Numrepval]=rowMeans(powerfisher)
    
    ######Holm on max - Competitor for replicability
    MyData$FWER[MyData$Method=="Holmmax"&MyData$Numrep==Numrepval]=rowMeans(fwerrepmax)
    MyData$Power[MyData$Method=="Holmmax"&MyData$Numrep==Numrepval]=rowMeans(powerrepmax)
    
    ######Holm on 32 hypotheses for global nulls
    MyData$FWER[MyData$Method=="CombinedHolmGlobal"&MyData$Numrep==Numrepval]=rowMeans(fwerholm32global)
    MyData$Power[MyData$Method=="CombinedHolmGlobal"&MyData$Numrep==Numrepval]=rowMeans(powerholm32global)
    
    ######Holm on 32 hypotheses - replicability
    MyData$FWER[MyData$Method=="CombinedHolmRep"&MyData$Numrep==Numrepval]=rowMeans(fwerholm32rep)
    MyData$Power[MyData$Method=="CombinedHolmRep"&MyData$Numrep==Numrepval]=rowMeans(powerholm32rep)
    
    ######Holm on 16 hypotheses on the combined data from the two studies
    MyData$FWER[MyData$Method=="HolmAllData"&MyData$Numrep==Numrepval]=rowMeans(fwerholmalldata)
    MyData$Power[MyData$Method=="HolmAllData"&MyData$Numrep==Numrepval]=rowMeans(powerholmalldata)
    
    
    ####Cross screening for global null and replicability
    MyData$FWER[MyData$Method=="Cross_glob"&MyData$Numrep==Numrepval]=rowMeans(fwerglobalcrossscreen)
    MyData$Power[MyData$Method=="Cross_glob"&MyData$Numrep==Numrepval]=rowMeans(powerglobalcrossscreen)
    MyData$FWER[MyData$Method=="Cross_rep"&MyData$Numrep==Numrepval]=rowMeans(fwerrepcrossscreen)
    MyData$Power[MyData$Method=="Cross_rep"&MyData$Numrep==Numrepval]=rowMeans(powerrepcrossscreen)
    
  }#End of iterations for Numrep
  
  return(list(MyData=MyData, fwerrepcrossscreen=fwerrepcrossscreen, powerrepcrossscreen=powerrepcrossscreen))
}

##Datagenfunc returns the values of I_11, I_10, I_01, I_00, and the mu vectors for the hypotheses in study 1 and in study 2
Datagenfunc<-function(Numrep, mu){
  
  if (Numrep==1){
    CESD4=rbind(c(1,0,0,0), c(1,0,0,0))
    PWB6=rbind(rep(0,6), rep(0,6))
    EWB1=c(0,0)
    Alcohol3=rbind(rep(0,3), rep(0,3))
    Smoking1=c(0,0)
    Health1=c(0,0)
  }
  if (Numrep==3){
    CESD4=rbind(c(1,1,1,0), c(1,1,1,0))
    PWB6=rbind(rep(0,6), rep(0,6))
    EWB1=c(0,0)
    Alcohol3=rbind(rep(0,3), rep(0,3))
    Smoking1=c(0,0)
    Health1=c(0,0)
     }
  if (Numrep==6){
    CESD4=rbind(c(1,1,1,1), c(1,1,1,1))
    PWB6=rbind(rep(0,6), rep(0,6))
    EWB1=c(0,0)
    Alcohol3=rbind(rep(0,3), rep(0,3))
    Smoking1=c(1,1)
    Health1=c(1,1)
     }
  if (Numrep==10){
    CESD4=rbind(c(1,1,1,1), c(1,1,1,1))
    PWB6=rbind(rep(0,6), rep(0,6))
    EWB1=c(1,1)
    Alcohol3=rbind(rep(1,3), rep(1,3))
    Smoking1=c(1,1)
    Health1=c(1,1)
      }
  if (Numrep==13){
    CESD4=rbind(c(1,1,1,1), c(1,1,1,1))
    PWB6=rbind(rep(1,6), rep(1,6))
    EWB1=c(0,0)
    Alcohol3=rbind(rep(1,3), rep(1,3))
    Smoking1=c(0,0)
    Health1=c(0,0)
     }
  if (Numrep==16){
    CESD4=rbind(c(1,1,1,1), c(1,1,1,1))
    PWB6=rbind(rep(1,6), rep(1,6))
    EWB1=c(1,1)
    Alcohol3=rbind(rep(1,3), rep(1,3))
    Smoking1=c(1,1)
    Health1=c(1,1)
    }
 
  muvec1=mu*c(CESD4[1,], PWB6[1,], EWB1[1], Alcohol3[1,], Smoking1[1], Health1[1])
  muvec2=mu*c(CESD4[2,], PWB6[2,], EWB1[2], Alcohol3[2,], Smoking1[2], Health1[2])
  
  refvec1=c(CESD4[1, ], PWB6[1,], EWB1[1], Alcohol3[1,], Smoking1[1], Health1[1])
  refvec2=c(CESD4[2, ], PWB6[2,], EWB1[2], Alcohol3[2,], Smoking1[2], Health1[2])
  refvecsum=refvec1+refvec2
  refvecminus=refvec1-refvec2
  I11ind=which(refvecsum==2)#indices of replicated signals
  I00ind=which(refvecsum==0)#no signals in both studies
  I10ind=which(refvecminus==1)
  I01ind=which(refvecminus==-1)
  return(list(I11ind=I11ind, I00ind=I00ind, I10ind=I10ind, I01ind=I01ind, muvec1=muvec1, muvec2=muvec2))
}

#Computepval computes the p-values for all the methods
Computepval<-function(datamat1, datamat2, m){
  pvalvec1right=rep(NA, m)
  pvalvec1left=rep(NA, m)
  pvalvec2right=rep(NA, m)
  pvalvec2left=rep(NA, m)
  pvalalldata=rep(NA, m)
  for(j in 1:m){
    pvalvec1right[j]=wilcox.test(datamat1[,j], alternative="greater")$p.value
    pvalvec1left[j]=wilcox.test(datamat1[,j], alternative="less")$p.value
    pvalvec2right[j]=wilcox.test(datamat2[,j], alternative="greater")$p.value
    pvalvec2left[j]=wilcox.test(datamat2[,j], alternative="less")$p.value
  pvalalldata[j]=wilcox.test(c(datamat1[,j],datamat2[,j]), alternative="two.sided")$p.value
    }
  pvalvec1twosided=2*pmin(pvalvec1right, pvalvec1left)
  pvalvec2twosided=2*pmin(pvalvec2right, pvalvec2left)
  return(list(pvalvec1right=pvalvec1right, pvalvec1left=pvalvec1left, pvalvec2right=pvalvec2right, pvalvec2left=pvalvec2left, pvalvec1twosided=pvalvec1twosided, pvalvec2twosided=pvalvec2twosided, pvalalldata=pvalalldata))
}


resdata=Globalfunc(1000, muvec, numrepvec, MyData, 0.05)
MyData=resdata$MyData
MyData$Numrep <- factor(MyData$Numrep,labels=c('1'=parse(text=TeX('$K_{11}=1$')),
                                              '3'=parse(text=TeX('$K_{11}=3$')),
                                               '6'=parse(text=TeX('$K_{11}=6$')),
                                              '10'=parse(text=TeX('$K_{11}=10')),
                                               '13'=parse(text=TeX('$K_{11}=13')),
                                               '16'=parse(text=TeX('$K_{11}=16'))))


pglob2=ggplot(data = MyData%>%filter(Metric=="Globalnull"),
         mapping = aes(x = Mu, y = Power, color = Method)) +
  geom_line() +
  facet_wrap(facets =  vars(Numrep), labeller=label_parsed)+
  theme(legend.title=element_blank(), legend.position="bottom")+
  labs(x=parse(text=TeX('$mu$')))+
  scale_color_manual(breaks=c(names=c("Fisher", "CombinedHolmGlobal", "Cross_glob", "HolmAllData")),
                     labels=c("Holm on Fisher's p-values", "Combined Holm: global nulls", "Cross-screening: global nulls", "Holm - full data"),
                     values=c("purple", "blue", "red", "green"))
plot(pglob2)


prep2=ggplot(data = MyData%>%filter(Metric=="Replicability"),
              mapping = aes(x = Mu, y = Power, color = Method)) +
  geom_line() +
  facet_wrap(facets =  vars(Numrep), labeller=label_parsed)+
  theme(legend.title=element_blank(), legend.position="bottom")+
  labs(x=parse(text=TeX('$mu$')), y="Power")+
  scale_color_manual(breaks=c(names=c("CombinedHolmRep", "Holmmax", "Cross_rep")),
                     labels=c("Combined Holm: replicability", "Holm on max p-values", "Cross-screening: replicability"),
                     values=c("blue", "purple", "red"))
  plot(prep2)
 
  ####Graphs for FWER###############
  pfwer1=ggplot(data = MyData%>%filter(Metric=="Globalnull"),
                mapping = aes(x = Mu, y = FWER, color = Method)) +
    geom_line() +
    facet_wrap(facets =  vars(Numrep), labeller=label_parsed)+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(x=parse(text=TeX('$mu$')))+
  scale_color_manual(breaks=c(names=c("Fisher", "CombinedHolmGlobal", "Cross_glob", "HolmAllData")),
                     labels=c("Holm on Fisher's p-values", "Combined Holm: global nulls", "Cross-screening: global nulls", "Holm - full data"),
                     values=c("purple", "blue", "red", "green"))
  plot(pfwer1)
  
  pfwer2=ggplot(data = MyData%>%filter(Metric=="Replicability"),
                mapping = aes(x = Mu, y = FWER, color = Method)) +
    geom_line() +
    facet_wrap(facets =  vars(Numrep), labeller=label_parsed)+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(x=parse(text=TeX('$mu$')))+
    scale_color_manual(breaks=c(names=c("CombinedHolmRep", "Holmmax", "Cross_rep")),
                       labels=c("Combined Holm: replicability", "Holm on max p-values", "Cross-screening: replicability"),
                       values=c("blue", "purple", "red"))
  plot(pfwer2)
  
