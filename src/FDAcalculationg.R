#y matrix
rm(list=ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(fda)
library(glmnet)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/X.csv")
Y <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/NFDB_point_stats/NFDB_point_txt/fireNOON.csv")
date <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/NFDB_point_stats/NFDB_point_txt/cldweek.csv")

empty <- matrix(data = 0, ncol = 52)
#value 52 min 
yset <- 1514:3488# 1975 = 25*79
#yset<- sample(10,yset)
for(i in yset){
  for (i in i){
    new <-rbind(empty,df$precipitation[(i-1):(i-52)]) # previous 52 week 
    # print(new)
    empty <- new
  }
}

Xmatrix <- new[-1,] # matrix value - X
fire <- Y[yset,2] # NOofFire - Y
y = log(fire)
y[y == "-Inf"]<- 0
mdf <- as.data.frame(Xmatrix)
x <- data.matrix(mdf[, c('V1','V2','V3','V4','V5',
                         'V6','V7','V8','V9','V10',
                         'V11','V12','V13','V14','V15',
                         'V16','V17','V18','V19','V20',
                         'V21','V22','V23','V24','V25',
                         'V26','V27','V28','V29','V30',
                         'V31','V32','V33','V34','V35',
                         'V36','V37','V38','V39','V40',
                         'V41','V42','V43','V44','V45',
                         'V46','V47','V48','V49','V50',
                         'V51','V52')])
cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)


#use fitted best model to make predictions
y_predicted <- predict(best_model, s = best_lambda, newx = x)

#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq


#y=Xmatrix[1,]# first week  week 300
argvals = seq(0,52, len = 52)
nbasis = 13
basisobj = create.bspline.basis(c(0,52),nbasis)
xs = smooth.basis(argvals=argvals, y=t(Xmatrix), fdParobj=basisobj) 
plot(xs)
x<-xs$fd # fda function

xlist = vector("list",2)
xlist[[1]] = rep(1,length(yset))
xlist[[2]] = x

# create a constant basis for the intercept
conbasis   = create.constant.basis(c(0,52))

#quartz()
#plot(conbasis)
# Define the small number of basis functions for beta(t)
# We do not use roughness penalty here
nbasis = 3
betabasis5 = create.fourier.basis(c(0,52),nbasis)
betalist1  = vector("list",2)
betalist1[[1]] = conbasis
betalist1[[2]] = betabasis5

fRegressList1 =  fRegress(y,xlist,betalist1)
names(fRegressList1)
betaestlist1  = fRegressList1$betaestlist
length(betaestlist1)
# betaestlist1 has two elements. The first element is the intercept
# The second element is the slope beta(t)

# obtain beta(t)
tempbetafd1   = betaestlist1[[2]]$fd

#quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempbetafd1, xlab="Week", ylab="β(t)")
x <- 1:52
y <- c(-0.0737061921,0,-0.0001000109,-0.0047236996,-0.0200412643,
       -0.0316449170,-0.0536387989,-0.0751471924,-0.1052093614,-0.0979901483,
       -0.1031142220,-0.1006618386,-0.1246870705,-0.1239538544,-0.1006120441,
       -0.1075900184,-0.0982814571,-0.1133192441,-0.1005294632,-0.1241472455,
       -0.1174767127,-0.1339048136,-0.1400092099,-0.1294226135,-0.1151585613,
       -0.0915167699,-0.0922951701,-0.0438828063,-0.0422427163,0,
       0.0108194370,0.0163352661,0.0205477270,0.0192191218,0.0175385084,
       0.0164200079,0.0258800617,0.0299715058,0.0361404753,0.0500885010,
       0.0507429211,0.0431320356,0.0930492908,0.1069075987,0.1084447673,
       0.1180122611,0.1390864073,0.1275272393,0.1084996437,0.1109395500,
       0.0738003552,0.0779780644)

plot(x, y)
par(new=TRUE)
lines(tempbetafd1,col='red', xlab="Week", ylab="β(t)")
abline(h = 0)
ny = length(y)
# obtain the intercept alpha
coef(betaestlist1[[1]])

# alpha is not equal to the mean of annualprec
mean(y)

# fitted value yhat
annualprechat1 = fRegressList1$yhatfdobj
# fitted residual
annualprecres1 = y - annualprechat1
# sum squared residuals
SSE1  = sum(annualprecres1^2)

# sum squared residuals for the null model y = alpha + \epsilon
SSE0    = sum((y - mean(y))^2)

# F test for the overall effect of x_i(t)
# H0: y = alpha + \epsilon
# H1: y = alpha + \int [beta(t)x(t)]dt + epsilon
(Fratio = ((SSE0-SSE1)/(nbasis-1)/(SSE1/(ny-nbasis))))

# 95% quantile of F(11,23)
qf(0.95,nbasis,ny-nbasis-1)

# Fratio >>qf(0.95,nbasis-1,ny-nbasis)
# indicating that x_i(t) has a significant effect on y_i

# calculate the p-value
1-pf(Fratio,nbasis-1,ny-nbasis)
# p-value indicating that x_i(t) has a significant effect on y_i



# Penalized Estimation

#Using the harmonic acceleration differential operator 
# to define roughness penalty on beta(t)
harmaccelLfd = int2Lfd(m=2)

# We use 35 Fourier basis functions to represent beta(t)
betabasis35 = create.fourier.basis(c(0, 52), 3)

# Choosing Smoothing Parameters using cross-validation

loglam = seq(3,13,1)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  print(paste("log lambda =", loglam[ilam]))
  lambda     = exp(loglam[ilam])
  # Define the small number of basis functions for beta(t)
  # We do not use roughness penalty here
  nbasis =ilam # seq(3,100,by = 5) 
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
  plot(tempbetafd1, xlab="Week",ylab="β(t)")
  x <- 1:52
  y <- c(-0.0737061921,0,-0.0001000109,-0.0047236996,-0.0200412643,
         -0.0316449170,-0.0536387989,-0.0751471924,-0.1052093614,-0.0979901483,
         -0.1031142220,-0.1006618386,-0.1246870705,-0.1239538544,-0.1006120441,
         -0.1075900184,-0.0982814571,-0.1133192441,-0.1005294632,-0.1241472455,
         -0.1174767127,-0.1339048136,-0.1400092099,-0.1294226135,-0.1151585613,
         -0.0915167699,-0.0922951701,-0.0438828063,-0.0422427163,0,
         0.0108194370,0.0163352661,0.0205477270,0.0192191218,0.0175385084,
         0.0164200079,0.0258800617,0.0299715058,0.0361404753,0.0500885010,
         0.0507429211,0.0431320356,0.0930492908,0.1069075987,0.1084447673,
         0.1180122611,0.1390864073,0.1275272393,0.1084996437,0.1109395500,
         0.0738003552,0.0779780644)
  
  plot(x, y)
  
  #predict 
  ytestfit = predict(fRegressList1, newdata=xlisttest$xtest, se.fit = FALSE)
  
  SSE.CV[ilam] = sum((ytest-ytestfit)^2) # Square Prediction Error  
  print(ilam) 
} 
plot(loglam,  SSE.CV, type="b", col=1,lty=1, lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex.lab=1,cex.axis=1)
#plot(log(lambda), cv.lam, type="l", col=2,lty=1, lwd=2,ylab="CV") 
SSE.CV.opt = lambda[which.min(SSE.CV )] 
SSE.CV.opt 

# Choose lambda which minimize SSE.CV
lambda      = 3
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.
# the first element of betalist2 is still the constant for the intercept

# do functional linear model
annPrecTemp    = fRegress(y, xlist, betalist2)

# get beta(t)
betaestlist2   = annPrecTemp$betaestlist

# get fitted value yhat
annualprechat2 = annPrecTemp$yhatfdobj

# get the effective degrees of freedom
print(annPrecTemp$df)

# do the F test 
# test the overall effect of x_i(t) on y_i
(SSE2 = sum((y-annualprechat2)^2))
(Fratio2 = ((SSE0-SSE2)/(annPrecTemp$df-1))/(SSE2/(ny-annPrecTemp$df)))
c(Fratio,Fratio2)
# Fratio2 < Fratio

# 95% quantile 
qf(0.95,annPrecTemp$df-1,ny-annPrecTemp$df)

# p-value
1-pf(Fratio2,annPrecTemp$df-1,ny-annPrecTemp$df)

# plot beta(t)
#quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(betaestlist2[[2]]$fd, xlab="Week",ylab="beta(t)",lwd=2,cex.lab=2,cex.axis=2)

# plot yhat vs. y
#quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(annualprechat2, y, lwd=2,cex.lab=2,cex.axis=2)
abline(lm(y~annualprechat2), lty='dashed', lwd=2)
abline(0,1,lty=1, lwd=2,col="red")





# Confidence Intervals

# fitted residuals
resid   = y - annualprechat2

# estimate sigma^2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))

# for smoothing temperature chat = a matrix * y
y2cMap  = xs$y2cMap

# obtain point-wise standard error for beta(t)
stderrList = fRegress.stderr(y, y2cMap, SigmaE)

betafdPar      = betaestlist2residbetafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

#quartz()
plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-0.2,0.1), lwd=2,cex.lab=2,cex.axis=2)
lines(betafd+2*betastderrfd, lty=2, lwd=2)
lines(betafd-2*betastderrfd, lty=2, lwd=2)
# The temperature has a significant positive effect on the annual precipitations from September to December. 
# There is no significant effect of the temperature on the annual precipitations in other time.






