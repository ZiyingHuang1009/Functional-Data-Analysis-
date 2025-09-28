rm(list=ls())
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)
library (plyr)

files <- list.files(pattern="\\.csv$")  #load the file names into the workspace
yourSum = matrix(0, nrow = 816, ncol = 34)  # create an empty means matrix
yourMeans = matrix(0, nrow = 816, ncol = 34)  # create an empty means matrix
df = matrix(0, nrow = 816, ncol = 33)  # create an empty means matrix
date = read.csv("1.csv")
# 1 - 3 

count = array(0, c(816, 34))
yourSum[,1] = NA
yourSum[,2] = date$V1
yourSum[,3] = date$V2
for(i in 1:length(files)){
  yourData <- read.csv(files[i])  # select the dataframe that needs to calculate the means
  for(iRow in 1:816)
  {
    for(iCol in 4:34)
    {
      
      #count[iRow, iCol] = 0
      if(is.na(yourData[iRow, iCol])) #df - yourData
        {
        count[iRow, iCol] = count[iRow, iCol] + 0
      yourSum[iRow, iCol] = yourSum[iRow, iCol] +  0 } else {count[iRow, iCol] = count[iRow, iCol] + 1
      yourSum[iRow, iCol] = yourSum[iRow, iCol] +  yourData[iRow, iCol]}
    }
    
    
  }
  print(i)

}
yourMeans[,4:34] = yourSum[,4:34]/count[,4:34]
yourMeans[,1] = date$V1
yourMeans[,2] = date$V2
df[,3:33] = format(round(yourMeans[,4:34], 2), nsmall = 2)
df[,1] = date$V1
df[,2] = date$V2
view(df) # final clean data
write.csv(df, "C:/Users/Jenny/Desktop/4F90/FDA/dataset/ON_DAILY/mean.csv", row.names=FALSE)

#———————— make a list for avg. of week number data ——————#


