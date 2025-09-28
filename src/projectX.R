rm(list=ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tibble)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/mean.csv")
dim(df)

new = df[3:33]
as.matrix(new)
final = t(new)
final<- final[!is.nan(final)]

prec = matrix(unlist(final), ncol = 7, byrow=T)
#write.csv(prec, "C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/weekly.csv", row.names=FALSE)

df <- read.csv("C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/weekly.csv")
precipitation = rowMeans(df) 
precipitation[3549] = 0.36 # only 1 value
precipitation = format(round(precipitation, 2), nsmall = 2)

weekX <- tibble(
  weeknumber = c(1:3549)
)
weekX <-
   add_column(weekX,precipitation, .after = 2)
write.csv(weekX, "C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/X.csv", row.names=FALSE)