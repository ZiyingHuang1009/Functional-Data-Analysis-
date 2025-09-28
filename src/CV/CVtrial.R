#y matrix
rm(list=ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(fda)
# Load packages
library(caret)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/X.csv")
Y <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/NFDB_point_stats/NFDB_point_txt/fireNOON.csv")
#date <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/NFDB_point_stats/NFDB_point_txt/cldweek.csv")

empty <- matrix(data = 0, ncol = 52)
yset1 = 1514:1908
yset2 = 1909:2303
yset3 = 2304:2698
yset4 = 2699:3093
yset5 = 3094:3488
yset = 1909:3488
groups <- c(quo(yset1), quo(yset2),quo(yset3),quo(yset4),quo(yset5))  # Create a quoture
#yset <- rbind(yset, yset1)
#yset <- rbind(yset, yset5) # 2-5
# 1975 = 5*395

for (i in seq_along(groups)) {
    print(i)
}
for(i in yset){
  for (i in i){
    new <-rbind(empty,df$precipitation[(i-1):(i-52)]) # previous 52 week 
    #print(new)
    empty <- new
  }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

Xmatrixtrain <- new[-1,] # matrix value - X train
firetrain <- Y[yset,2] # NOofFire - Y
ytrain = log(firetrain)
ytrain[ytrain == "-Inf"]<- 0
new <- NULL
empty <- NULL
argvals = seq(0,52, len = 52)
nbasis = 13
basisobj = create.bspline.basis(c(0,52),nbasis)
#plot(basisobj)
xstrain = smooth.basis(argvals=argvals, y=t(Xmatrixtrain), fdParobj=basisobj) 
xtrain<-xstrain$fd # fda function

xlisttrain = vector("list",2)
xlisttrain[[1]] = rep(1,length(yset))
xlisttrain[[2]] = xtrain

# test

# 1975 = 5*395
for(i in yset1){
  for (i in i){
    new <-rbind(empty,df$precipitation[(i-1):(i-52)]) # previous 52 week 
    #print(new)
    empty <- new
  }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

Xmatrixtest <- new # matrix value - X
firetest <- Y[yset1,2] # NOofFire - Y
ytest = log(firetest)
ytest[ytest == "-Inf"]<- 0
xstest = smooth.basis(argvals=argvals, y=t(Xmatrixtest), fdParobj=basisobj) 
xtest<-xstest$fd # fda function

xlisttest = vector("list",2)
xlisttest[[1]] = rep(1,length(yset1))
xlisttest[[2]] = xtest

# create a constant basis for the intercept
conbasis   = create.constant.basis(c(0,52))
#quartz()
#plot(conbasis)

#--------------- Select the optimal lambda ---------------# 

lambda = seq(3,13,length.out =10)
cv.lam = rep(NA, length(lambda)) 
for(iterlambda in 1:length(lambda)) 
{ 
  # Define the small number of basis functions for beta(t)
  # We do not use roughness penalty here
  nbasis = iterlambda # seq(3,100,by = 5) 
  betabasis5 = create.fourier.basis(c(0,52),nbasis)
  betalist1  = vector("list",2)
  betalist1[[1]] = conbasis
  betalist1[[2]] = betabasis5
  
  fRegressList1 =  fRegress(ytrain,xlisttrain,betalist1)
  names(fRegressList1)
  betaestlist1  = fRegressList1$betaestlist
  length(betaestlist1)
  # betaestlist1 has two elements. The first element is the intercept
  # The second element is the slope beta(t)
  
  # obtain beta(t)
  tempbetafd1   = betaestlist1[[2]]$fd
  
  #quartz()
  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(tempbetafd1, xlab="Precipitation", ylab="Beta for Fire")
  
  #predict 
  ytestfit = predict(fRegressList1, newdata=xlisttest$xtest, se.fit = FALSE)
  
  cv.lam[iterlambda] = sum((ytest-ytestfit)^2) # Square Prediction Error  
  print(iterlambda) 
} 
plot(log(lambda), cv.lam, type="l", col=2,lty=1, lwd=2,ylab="CV") 
cv.lam.opt = lambda[which.min(cv.lam)] 
cv.lam.opt 

#L = Xmatrixtrain%*%ginv(t(Xmatrixtrain)%*%Xmatrixtrain+cv.lam.opt)%*%t(Xmatrixtrain) 
#yhatcv = L%*%yset

#plot(X,yset,main="CV",xlab="log(Precipitation)",ylab="log(Fire)") 
#lines(x,yhatcv,col="green", lwd=2) 
