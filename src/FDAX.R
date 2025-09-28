#lm smoothing
rm(list=ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(npreg)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/X.csv")
view(df)

i = 1514


  
form = sprintf(paste0('Smoothing Week No.',i))
x <- seq(1,52)
y <- df$precipitation[i:(i+51)]
d <- data.frame(x = x, y = y)


  Smoothing <- (ggplot(d, aes(x, y)) 
                + stat_smooth(method='lm', formula = y~poly(x,4), color="green",se = FALSE)
      )
print(Smoothing)
