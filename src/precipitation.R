rm(list=ls())
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)

#open all dataset files
temp = list.files(pattern="\\.txt$")
for (i in 1:length(temp)) assign(temp[i], read.delim(temp[i]))

for (i in 1:length(temp)) 
u[i] <- temp[i]
widths <- c(4,4,rep(9, 31)) 
data0[i] = read.fwf(u[i], skip = 1, widths) 

# data from 1950.1 - 2017.12
df[i]<-tail(data0[i],n=816)


df[i][df[i] == "-9999.99M"] <- "NA"

#remove T
for(i in 3:33){
  df[,i]<-gsub("T","",as.character(df[,i]))
}




#calculate each average
