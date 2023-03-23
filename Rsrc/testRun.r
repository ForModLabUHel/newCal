

####Run model in parallel
library(ggplot2)
# library(Metrics)
library(data.table)
library(parallel)
library(BayesianTools)
library(ggpubr)

####set project library if running on CSC
vLocal <- FALSE #### flag for runs on CSC(FALSE) or on laptop(TRUE)
if(!vLocal){
  setwd("/scratch/project_2000994/calibrations/srinet/newCal")
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

source('Rsrc/functions.r') ###run using url
source('Rsrc/settings.r') ###run using url

# load("inputs/initPrebas_old growth.rdata")
load("inputs/init_set1.rdata")
load("inputs/init_set2.rdata")
load("inputs/init_set3.rdata")
load("inputs/init_set4.rdata")
load("inputs/init_set5Flux.RData")

vapu_S<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_spruce.csv'))
nData_S <- length(vapu_S$plotNo)
vapu_P<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_pine.csv'))
nData_P <- length(vapu_P$plot)
#load("outCal/pMAP.rdata")
par<-read.csv('inputs/par_prebas_newCal.csv')
parmod<-par$parmod

# modOut <- list()
# 
###Run model 
startX <- Sys.time()
likX <- likelihood(parmod)
endX <- Sys.time()
timeX = endX- startX
print("parallel runs")
print(likX)
print(timeX)

startX <- Sys.time()
ll1 <- likelihood1(parmod)
ll2 <- likelihood2(parmod)
ll3 <- likelihood3(parmod)
ll4 <- likelihood4(parmod)
ll5 <- likelihood5Flux(parmod)
llS <- ll1 + ll2 +ll3 + ll4 + ll5
endX <- Sys.time()
timeX = endX- startX
print("sequential runs")
print(llS)
print(timeX)
