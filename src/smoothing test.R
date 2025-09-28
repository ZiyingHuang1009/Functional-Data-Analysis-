#splines smoothing
rm(list=ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(fda)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/X.csv")
i<-c(1514:3488) # 1979/12/30——2017/11/4
j=sample(i, 10, replace=TRUE)
for (j in j) {
  
form = sprintf(paste0('Smoothing Week No.',j))
x <- seq(1,52)
y <- df$precipitation[j:(j+51)]

d <- data.frame(x = x, y = y)}


argvals = seq(0,52, len = 52)
nbasis = 13
basisobj = create.bspline.basis(c(0,52),nbasis)
plot(basisobj)
ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj)

# y matrix



# fit using ss
mod.ss <- ss(x, y, nknots = 10)
mod.ss

# fit using smooth.spline
mod.smsp <- smooth.spline(x, y, nknots = 50)
mod.smsp


# GCV selection
plot(x,y,xlab = form, ylab = "Precipitation") # I just do not want to show the plot
mod.ss <- ss(x, y, all.knots = TRUE,lambda = 1e-6)
line <- lines(mod.ss, ylim = c(0,10), lwd = 2,col="blue")# add estimated function with lambda chosen by GCV 

