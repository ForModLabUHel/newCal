#1# specify packages folder if you are on puhti
vLocal <- FALSE
if(!vLocal){
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}

#1#
library(ggplot2)
library(Metrics)
library(data.table)
library(BayesianTools)
library(ggpubr)
library(dplyr)
library(ggpmisc)
library(broom)

#2# install package if needed
library(devtools)
#remove.packages("Rprebasso")
vPREB <- "newVersion"
devtools::install_github("ForModLabUHel/Rprebasso", ref=vPREB)
#devtools::install_github("ForModLabUHel/Rprebasso")
library(Rprebasso)

#3# setting working directory
if(vLocal){
  #you need to change this if you want to work on your local machine
  setwd("Z:/PREBAS_calibration/newCal/")
}else{
  setwd("/scratch/project_2000994/calibrations/srinet/newCal/")
}

###loading functions and settings
source('Rsrc/functions.r')
source('Rsrc/settings.r')
source("https://raw.github.com/ForModLabUHel/utilStuff/master/ErrorDecomposition/ErrorDecomposition.R")
###load data for initialization
load('inputs/init_set1.rdata')
load('inputs/init_set2.rdata')
load('inputs/init_set3.rdata')
load('inputs/init_set4.rdata')
load('inputs/init_set5Flux.RData')

vapu_S<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_spruce.csv'))
nData_S <- length(vapu_S$plotNo)
vapu_P<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_pine.csv'))
nData_P <- length(vapu_P$plot)

#load("outCal/pMAP.rdata")
par<-read.csv('inputs/par_prebas_newCal.csv')
parmod<-par$parmod

startX <- Sys.time()
modOut1 <- likelihood1(parmod,cal=F)
modOut2 <- likelihood2(parmod,cal=F)
modOut3 <- likelihood3(parmod,cal=F)
modOut4 <- likelihood4(parmod,cal=F)
modOut5 <- likelihood5Flux(parmod,cal=F)
endX <- Sys.time()
timeX = endX- startX
print(timeX)

speciesNam <- c("pine","spruce","birch")
fsiteNam<-c("Värriö","Lettosuo")

allData <- data.table(obs = modOut1$obsV,sim=modOut1$simV,cal_set='set 1',var="V", speciesID=speciesNam[extractSpecies(Vdata_s1,modOut1)], siteID = extractSite(Vdata_s1,modOut1))
allData <- rbind(allData, data.table(obs = modOut2$obsV,sim=modOut2$simV,cal_set='set 2',var="V",speciesNam[extractSpecies(Vdata_s2,modOut2)], extractSite(Vdata_s2,modOut2)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut3$obsV,sim=modOut3$simV,cal_set='set 3',var="V",speciesNam[extractSpecies(Vdata_s3,modOut3)], extractSite(Vdata_s3,modOut3)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$obsV,sim=modOut4$simV,cal_set='set 4',var="V",speciesNam[extractSpecies(Vdata_s4,modOut4)], extractSite(Vdata_s4,modOut4)), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut1$obsB,sim=modOut1$simB,cal_set='set 1',var="B",speciesNam[extractSpecies(Bdata_s1,modOut1)], extractSite(Bdata_s1,modOut1)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut2$obsB,sim=modOut2$simB,cal_set='set 2',var="B",speciesNam[extractSpecies(Bdata_s2,modOut2)], extractSite(Bdata_s2,modOut2)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut3$obsB,sim=modOut3$simB,cal_set='set 3',var="B",speciesNam[extractSpecies(Bdata_s3,modOut3)], extractSite(Bdata_s3,modOut3)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$obsB,sim=modOut4$simB,cal_set='set 4',var="B",speciesNam[extractSpecies(Bdata_s4,modOut4)], extractSite(Bdata_s4,modOut4)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsB,sim=modOut5$simB,cal_set='set 5',var="B",speciesNam[extractSpecies(Bdata_s5,modOut5)], extractSite(Bdata_s5,modOut5)), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut1$obsD,sim=modOut1$simD,cal_set='set 1',var="D",speciesNam[extractSpecies(Ddata_s1,modOut1)], extractSite(Ddata_s1,modOut1)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut2$obsD,sim=modOut2$simD,cal_set='set 2',var="D",speciesNam[extractSpecies(Ddata_s2,modOut2)], extractSite(Ddata_s2,modOut2)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut3$obsD,sim=modOut3$simD,cal_set='set 3',var="D",speciesNam[extractSpecies(Ddata_s3,modOut3)], extractSite(Ddata_s3,modOut3)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$obsD,sim=modOut4$simD,cal_set='set 4',var="D",speciesNam[extractSpecies(Ddata_s4,modOut4)], extractSite(Ddata_s4,modOut4)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsD,sim=modOut5$simD,cal_set='set 5',var="D",speciesNam[extractSpecies(Ddata_s5,modOut5)], extractSite(Ddata_s5,modOut5)), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut1$obsH,sim=modOut1$simH,cal_set='set 1',var="H",speciesNam[extractSpecies(Hdata_s1,modOut1)], extractSite(Hdata_s1,modOut1)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut2$obsH,sim=modOut2$simH,cal_set='set 2',var="H",speciesNam[extractSpecies(Hdata_s2,modOut2)], extractSite(Hdata_s2,modOut2)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut3$obsH,sim=modOut3$simH,cal_set='set 3',var="H",speciesNam[extractSpecies(Hdata_s3,modOut3)], extractSite(Hdata_s3,modOut3)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$obsH,sim=modOut4$simH,cal_set='set 4',var="H",speciesNam[extractSpecies(Hdata_s4,modOut4)], extractSite(Hdata_s4,modOut4)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsH,sim=modOut5$simH,cal_set='set 5',var="H",speciesNam[extractSpecies(Hdata_s5,modOut5)], extractSite(Hdata_s5,modOut5)), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut3$obsHc,sim=modOut3$simHc,cal_set='set 3',var="Hc",speciesNam[extractSpecies(Hcdata_s3,modOut3)], extractSite(Hcdata_s3,modOut3)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$obsHc,sim=modOut4$simHc,cal_set='set 4',var="Hc",speciesNam[extractSpecies(Hcdata_s4,modOut4)], extractSite(Hcdata_s4,modOut4)), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsHc,sim=modOut5$simHc,cal_set='set 5',var="Hc",speciesNam[extractSpecies(Hcdata_s5,modOut5)], extractSite(Hcdata_s5,modOut5)), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut5$obsGPP,sim=modOut5$simGPP,cal_set='set 5',var="GPP","all",fsiteNam[GPPdata_s5$outData[,1]]), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsET,sim=modOut5$simET,cal_set='set 5',var="ET","all",fsiteNam[ETdata_s5$outData[,1]]), use.names=FALSE)

allData<-na.omit(allData)

r2_rmse <- function(df){
  summ <- df %>% 
    group_by(speciesID) %>% 
    summarise(Rsq = cor(obs, sim)^2,
              RMSE = rmse(obs, sim)) %>% 
    mutate_if(is.numeric, round, digits=2)
  df.annotations <- data.frame()
  df.annotations <- rbind(df.annotations,
                            cbind(as.character(summ$speciesID),
                                  paste0("Rsq = ", summ$Rsq, ", RMSE = ", summ$RMSE)))
  colnames(df.annotations) <- c("speciesID", "label")
  return(df.annotations)
}
r2_rmse_flux <- function(df){
  summ <- df %>% 
    group_by(siteID) %>% 
    summarise(Rsq = cor(obs, sim)^2,
              RMSE = rmse(obs, sim)) %>% 
    mutate_if(is.numeric, round, digits=2)
  df.annotations <- data.frame()
  df.annotations <- rbind(df.annotations,
                          cbind(as.character(summ$siteID),
                                paste0("Rsq = ", summ$Rsq, ", RMSE = ", summ$RMSE)))
  colnames(df.annotations) <- c("siteID", "label")
  return(df.annotations)
}

calc_resid<- function(df){
  r <- df %>%
    group_by(speciesID) %>%
    do(data.frame(., augment(lm(sim ~ obs, data=.))[c(3:8)]))
  return(r)
}
calc_residF<- function(df){
  r <- df %>%
    group_by(siteID) %>%
    do(data.frame(., augment(lm(sim ~ obs, data=.))[c(3:8)]))
  return(r)
}

pdf(file="out/outupdate.pdf")

pV <- ggplot(allData[var=="V"],aes(x=obs,y=sim,col=speciesID)) +
  geom_point() + geom_abline() + ggtitle("V")+
  geom_text(data=r2_rmse(allData[var=="V"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4,6), size=4)
# pV_res<-ggplot(calc_resid(allData[var=="V"]),aes(x=obs,y=sim,col=speciesID)) + 
#   geom_smooth(method = "lm", se = FALSE) +
#   geom_segment(aes(xend = obs, yend = .fitted)) +
#   geom_point() 
pV_res <- ggplot(calc_resid(allData[var=="V"]), aes(x = .fitted, y = .resid, col=speciesID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='V: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pB <- ggplot(allData[var=="B"],aes(x=obs,y=sim,col=speciesID)) +
  geom_point() + geom_abline() + ggtitle("B")+ 
  geom_text(data=r2_rmse(allData[var=="B"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4,6), size=4)
pB_res <- ggplot(calc_resid(allData[var=="B"]), aes(x = .fitted, y = .resid, col=speciesID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='B: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pD <- ggplot(allData[var=="D"],aes(x=obs,y=sim,col=speciesID)) +
  geom_point() + geom_abline() + ggtitle("D")+ 
  geom_text(data=r2_rmse(allData[var=="D"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4,6), size=4)
pD_res <- ggplot(calc_resid(allData[var=="D"]), aes(x = .fitted, y = .resid, col=speciesID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='D: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pH <- ggplot(allData[var=="H"],aes(x=obs,y=sim,col=speciesID)) +
  geom_point() + geom_abline() + ggtitle("H")+ 
  geom_text(data=r2_rmse(allData[var=="H"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4,6), size=4)
pH_res <- ggplot(calc_resid(allData[var=="H"]), aes(x = .fitted, y = .resid, col=speciesID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='H: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pHc <- ggplot(allData[var=="Hc"],aes(x=obs,y=sim,col=speciesID)) +
  geom_point() + geom_abline() + ggtitle("Hc")+ 
  geom_text(data=r2_rmse(allData[var=="Hc"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4,6), size=4)
pHc_res <- ggplot(calc_resid(allData[var=="Hc"]), aes(x = .fitted, y = .resid, col=speciesID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='Hc: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pGPP <- ggplot(allData[var=="GPP"],aes(x=obs,y=sim,col=siteID)) +
  geom_point() + geom_abline() + ggtitle("GPP")+ 
  geom_text(data=r2_rmse_flux(allData[var=="GPP"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4), size=4)
pGPP_res <- ggplot(calc_residF(allData[var=="GPP"]), aes(x = .fitted, y = .resid, col=siteID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='GPP: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')

pET <- ggplot(allData[var=="ET"],aes(x=obs,y=sim,col=siteID)) +
  geom_point() + geom_abline() + ggtitle("ET")+ 
  geom_text(data=r2_rmse_flux(allData[var=="ET"]), aes(x=-Inf,y=+Inf,label=label),
            hjust = -0.3, vjust = c(2,4), size=4)
pET_res <- ggplot(calc_residF(allData[var=="ET"]), aes(x = .fitted, y = .resid, col=siteID)) + 
  geom_point() + geom_hline(yintercept = 0) + 
  labs(title='ET: Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
print(pV)
print(pB)
print(pD)
print(pH)
print(pHc)
print(pGPP)
print(pET)
dev.off()

pdf(file="out/outupdate_residuals.pdf")
print(pV_res)
print(pB_res)
print(pD_res)
print(pH_res)
print(pHc_res)
print(pGPP_res)
print(pET_res)
dev.off()


### error decomposition

resultV1d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='pine')],method = 1)
resultH1d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='pine')],method = 1)
resultD1d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='pine')],method = 1)
resultB1d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='pine')],method = 1)
resultHc1d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='pine')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='pine')],method = 1)

resultV2d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='spruce')],method = 1)
resultH2d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='spruce')],method = 1)
resultD2d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='spruce')],method = 1)
resultB2d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='spruce')],method = 1)
resultHc2d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='spruce')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='spruce')],method = 1)

resultV3d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='birch')],method = 1)
resultH3d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='birch')],method = 1)
resultD3d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='birch')],method = 1)
resultB3d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='birch')],method = 1)
resultHc3d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='birch')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='birch')],method = 1)

resultGd1<- MSEdec("GPP",allData$obs[which(allData$var=="GPP")],
                   allData$sim[which(allData$var=="GPP")],method = 1)
resultEd1<- MSEdec("ET",allData$obs[which(allData$var=="ET")],
                   allData$sim[which(allData$var=="ET")],method = 1)
method_Kob <- Map(c,resultV1d1,resultH1d1,resultD1d1,resultB1d1,resultHc1d1,
                  resultV2d1,resultH2d1,resultD2d1,resultB2d1,resultHc2d1,
                  resultV3d1,resultH3d1,resultD3d1,resultB3d1,resultHc3d1,
                  resultGd1,resultEd1)

method_Kob <- data.frame(matrix(unlist(method_Kob), nrow=length(method_Kob), byrow=T))
row.names(method_Kob)<-c("Var","sb","sdsd","lc","mse")
write.csv(method_Kob,'Z:/PREBAS_calibration/newCal/out/errordecomp_kobsal.csv')


resultV1d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='pine')],method = 2)
resultH1d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='pine')],method = 2)
resultD1d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='pine')],method = 2)
resultB1d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='pine')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='pine')],method = 2)
resultHc1d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='pine')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='pine')],method = 2)

resultV2d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='spruce')],method = 2)
resultH2d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='spruce')],method = 2)
resultD2d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='spruce')],method = 2)
resultB2d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='spruce')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='spruce')],method = 2)
resultHc2d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='spruce')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='spruce')],method = 2)

resultV3d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="V" & allData$speciesID=='birch')],method = 2)
resultH3d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="H" & allData$speciesID=='birch')],method = 2)
resultD3d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="D" & allData$speciesID=='birch')],method = 2)
resultB3d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='birch')],
                    allData$sim[which(allData$var=="B" & allData$speciesID=='birch')],method = 2)
resultHc3d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='birch')],
                     allData$sim[which(allData$var=="Hc" & allData$speciesID=='birch')],method = 2)

resultGd2<- MSEdec("GPP",allData$obs[which(allData$var=="GPP")],
                    allData$sim[which(allData$var=="GPP")],method = 2)
resultEd2<- MSEdec("ET",allData$obs[which(allData$var=="ET")],
                    allData$sim[which(allData$var=="ET")],method = 2)

method_Gauch <- Map(c,resultV1d2,resultH1d2,resultD1d2,resultB1d2,resultHc1d2,
                  resultV2d2,resultH2d2,resultD2d2,resultB2d2,resultHc2d2,
                  resultV3d2,resultH3d2,resultD3d2,resultB3d2,resultHc3d2,
                  resultGd2,resultEd2)

method_Gauch <- data.frame(matrix(unlist(method_Gauch), nrow=length(method_Gauch), byrow=T))
row.names(method_Gauch)<-c("Var","sb","sdsd","lc","mse")
write.csv(method_Gauch,'Z:/PREBAS_calibration/newCal/out/errordecomp_gauch.csv')


