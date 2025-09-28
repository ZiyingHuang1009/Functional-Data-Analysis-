rm(list=ls())
library(dplyr)
library(tidyverse)
library(stringr)

temp = list.files(pattern="\\.txt$")
for (j in 1:length(temp)) assign(temp[j], read.delim(temp[j]))

for (j in 1:57){
u <- temp[j]
widths <- c(4,4,rep(9, 31)) 
data0 = read.fwf(u, skip = 1, widths) 

# data from 1950.1 - 2017.12
df<-tail(data0,n=816)

#remove T\X\Y\A\C\L\Z\E\F
for(i in 3:33){
  df[,i]<-gsub("T","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("X","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("Y","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("A","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("C","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("L","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("Z","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("E","",as.character(df[,i]))
}
for(i in 3:33){
  df[,i]<-gsub("F","",as.character(df[,i]))
}


df[df == "-9999.99M"] <- "NA"
df[df == "-9999.99"] <- "NA"
form = sprintf('new//%s.csv', j)
write.csv(df, file = form)
}