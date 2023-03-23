#1# specify packages folder if you are on puhti
vLocal <- FALSE
newV=F
if(!vLocal){
  .libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
}
#1#

library(ggplot2)
# library(Metrics)
library(data.table)
library(BayesianTools)
library(ggpubr)
library(dplyr)
library(ggpmisc)
library(broom)

#2# install package if needed
if(newV){
  devtools::install_github("ForModLabUHel/Rprebasso", ref="newVersion")
}else{
  devtools::install_github("ForModLabUHel/Rprebasso", ref="master")
}
library(Rprebasso)

#3# setting working directory
if(vLocal){
  #you need to change this if you want to work on your local machine
  setwd("Z:/PREBAS_calibration/newCal")
}else{
  setwd("/scratch/project_2000994/calibrations/srinet/newCal")
}

###loading functions and settings
source("Rsrc/functions.r")
source("Rsrc/settings.r")
source("https://raw.github.com/ForModLabUHel/utilStuff/master/ErrorDecomposition/ErrorDecomposition.R")
## Warning: unsupported URL scheme
###load data for initialization
if(newV){
  load('inputs/init_set1_newVersion.rdata')
  load('inputs/init_set2_newVersion.rdata')
  load('inputs/init_set3_newVersion.rdata')
  load('inputs/init_set4_newVersion.rdata')
  load('inputs/init_set5Flux_newVersion.rdata')
}else{
  load('inputs/init_set1.rdata')
  load('inputs/init_set2.rdata')
  load('inputs/init_set3.rdata')
  load('inputs/init_set4.rdata')
  load('inputs/init_set5Flux.RData')
}

vapu_S<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_spruce.csv'))
nData_S <- length(vapu_S$plotNo)
vapu_P<-read.csv(url('https://raw.githubusercontent.com/ForModLabUHel/newCal/master/inputs/VAPU_pine.csv'))
nData_P <- length(vapu_P$plot)

#load("outCal/pMAP.rdata")
par<-read.csv('inputs/par_prebas_newCal.csv')
parmod<-par$parmod

# init_set1$ECMmod=0; init_set1$pCN_alfar=NULL
# init_set2$ECMmod=0; init_set2$pCN_alfar=NULL
# init_set3$ECMmod=0; init_set3$pCN_alfar=NULL
# init_set4$ECMmod=0; init_set4$pCN_alfar=NULL
# init_set5Flux$ECMmod=0; init_set5Flux$pCN_alfar=NULL

startX <- Sys.time()
modOut1 <- likelihood1(parmod,cal=F)
modOut2 <- likelihood2(parmod,cal=F)
modOut3 <- likelihood3(parmod,cal=F)
modOut4 <- likelihood4(parmod,cal=F)
modOut5 <- likelihood5Flux(parmod,cal=F)
endX <- Sys.time()
timeX = endX- startX
print(timeX)

gpp_annual<-read.csv("inputs/gpp_annual.csv")

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

allData <- rbind(allData, data.table(obs = modOut4$wf_p_obs,sim=modOut4$simWf1_p,cal_set='set 4',var="Wf1","pine", vapu_P$plot), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$wf_p_obs,sim=modOut4$simWf2_p,cal_set='set 4',var="Wf2","pine", vapu_P$plot), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$As_p_obs,sim=modOut4$simAs_p,cal_set='set 4',var="As","pine", vapu_P$plot), use.names=FALSE)

allData <- rbind(allData, data.table(obs = modOut4$wf_s_obs,sim=modOut4$simWf1_s,cal_set='set 4',var="Wf1","spruce", vapu_S$plot), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$wf_s_obs,sim=modOut4$simWf2_s,cal_set='set 4',var="Wf2","spruce", vapu_S$plot), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut4$As_s_obs,sim=modOut4$simAs_s,cal_set='set 4',var="As","spruce", vapu_S$plot), use.names=FALSE)


allData <- rbind(allData, data.table(obs = modOut5$obsGPP_day,sim=modOut5$simGPP_day,cal_set='set 5',var="GPP","all",fsiteNam[GPPdata_s5$outData[,1]]), use.names=FALSE)
allData <- rbind(allData, data.table(obs = modOut5$obsET,sim=modOut5$simET,cal_set='set 5',var="ET","all",fsiteNam[ETdata_s5$outData[,1]]), use.names=FALSE)

allData <- rbind(allData, data.table(obs = gpp_annual$GPPobs,sim=modOut5$simGPP_yr,cal_set='set 5',var="GPP_yr","all",fsiteNam[gpp_annual$siteID]), use.names=FALSE)

allData<-na.omit(allData)
# allData_master<-allData
# allData_newPar<-allData
# allData_newV<-allData
# allData_newV_Par<-allData
#
# save(allData_master, file="out/output/allData_master.RData")
# save(allData_newPar, file="out/output/allData_newPar.RData")
# save(allData_newV, file="out/output/allData_newV.RData")
# save(allData_newV_Par, file="out/output/allData_newV_Par.RData")
# load("out/output/allData_master.RData")


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

# pdf(file="out/output_all_compare1.pdf")
# 
# ### MASTER VERSION -- OBSERVED VS. PREDICTED
# pV_o <- ggplot(allData_master[var=="V"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("V - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="V"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 800)
# # pV_res<-ggplot(calc_resid(allData_master[var=="V"]),aes(x=obs,y=sim,col=speciesID)) +
# #   geom_smooth(method = "lm", se = FALSE) +
# #   geom_segment(aes(xend = obs, yend = .fitted)) +
# #   geom_point()
# pV_res_o <- ggplot(calc_resid(allData_master[var=="V"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='V - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pB_o <- ggplot(allData_master[var=="B"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("B - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="B"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 60)
# pB_res_o <- ggplot(calc_resid(allData_master[var=="B"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='B - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pD_o <- ggplot(allData_master[var=="D"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("D - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="D"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 75)
# pD_res_o <- ggplot(calc_resid(allData_master[var=="D"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='D - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pH_o <- ggplot(allData_master[var=="H"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("H - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="H"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 45)
# pH_res_o <- ggplot(calc_resid(allData_master[var=="H"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='H - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pHc_o <- ggplot(allData_master[var=="Hc"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Hc - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="Hc"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 25)
# pHc_res_o <- ggplot(calc_resid(allData_master[var=="Hc"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Hc - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf1_o <- ggplot(allData_master[var=="Wf1"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf1 - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="Wf1"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 40)
# pWf1_res_o <- ggplot(calc_resid(allData_master[var=="Wf1"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf1 - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf2_o <- ggplot(allData_master[var=="Wf2"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf2 - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="Wf2"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 25)
# pWf2_res_o <- ggplot(calc_resid(allData_master[var=="Wf2"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf2 - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pAs_o <- ggplot(allData_master[var=="As"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("As - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_master[var=="As"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 0.1)
# pAs_res_o <- ggplot(calc_resid(allData_master[var=="As"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   labs(title='As - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pGPP_o <- ggplot(allData_master[var=="GPP"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("GPP - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_master[var=="GPP"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pGPP_res_o <- ggplot(calc_residF(allData_master[var=="GPP"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP -  - current version\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pET_o <- ggplot(allData_master[var=="ET"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("ET - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_master[var=="ET"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pET_res_o <- ggplot(calc_residF(allData_master[var=="ET"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='ET -  - current version\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# dt_o<-allData_master[var=="GPP_yr"]
# dt_o$year<-c(1,2,3,1,2,3,4,5,6)
# pGPPyr_o <- ggplot(dt_o, aes(x= year)) + geom_point(aes(y=obs, col=siteID)) + 
#   geom_line(aes(y=sim, col=siteID)) +ylim(0,1500) + ggtitle("GPP annual - current version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10))
# pGPPyr_res_o <- ggplot(calc_residF(allData_master[var=="GPP_yr"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP annual - current version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# ### MASTER VERSION - NEW PARAMETERS -- OBSERVED VS. PREDICTED
# pV_newP <- ggplot(allData_newPar[var=="V"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("V - current version\nsuggested pValues\n") +
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="V"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 800)
# pV_res_newP <- ggplot(calc_resid(allData_newPar[var=="V"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='V - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pB_newP <- ggplot(allData_newPar[var=="B"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("B - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="B"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 60)
# pB_res_newP <- ggplot(calc_resid(allData_newPar[var=="B"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='B - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pD_newP <- ggplot(allData_newPar[var=="D"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("D - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="D"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 75)
# pD_res_newP <- ggplot(calc_resid(allData_newPar[var=="D"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='D - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pH_newP <- ggplot(allData_newPar[var=="H"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("H - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="H"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 45)
# pH_res_newP <- ggplot(calc_resid(allData_newPar[var=="H"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='H - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pHc_newP <- ggplot(allData_newPar[var=="Hc"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Hc - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="Hc"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 25)
# pHc_res_newP <- ggplot(calc_resid(allData_newPar[var=="Hc"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Hc - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf1_newP <- ggplot(allData_newPar[var=="Wf1"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf1 - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="Wf1"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 40)
# pWf1_res_newP <- ggplot(calc_resid(allData_newPar[var=="Wf1"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf1 - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf2_newP <- ggplot(allData_newPar[var=="Wf2"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf2 - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="Wf2"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 25)
# pWf2_res_newP <- ggplot(calc_resid(allData_newPar[var=="Wf2"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf2 - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pAs_newP <- ggplot(allData_newPar[var=="As"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("As - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newPar[var=="As"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 0.1)
# pAs_res_newP <- ggplot(calc_resid(allData_newPar[var=="As"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='As - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pGPP_newP <- ggplot(allData_newPar[var=="GPP"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("GPP - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newPar[var=="GPP"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pGPP_res_newP <- ggplot(calc_residF(allData_newPar[var=="GPP"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pET_newP <- ggplot(allData_newPar[var=="ET"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("ET - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newPar[var=="ET"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pET_res_newP <- ggplot(calc_residF(allData_newPar[var=="ET"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='ET - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# dt_newP<-allData_newPar[var=="GPP_yr"]
# dt_newP$year<-c(1,2,3,1,2,3,4,5,6)
# pGPPyr_newP <- ggplot(dt_newP, aes(x= year)) + geom_point(aes(y=obs, col=siteID)) + 
#   geom_line(aes(y=sim, col=siteID)) +ylim(0,1500) + ggtitle("GPP annual - current version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10))
# pGPPyr_res_newP <- ggplot(calc_residF(allData_newPar[var=="GPP_yr"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP annual - current version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# ###NEW VERSION -- OBSERVED VS. PREDICTED
# pV <- ggplot(allData_newV[var=="V"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("V - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="V"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 800)
# pV_res <- ggplot(calc_resid(allData_newV[var=="V"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='V - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pB <- ggplot(allData_newV[var=="B"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("B - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="B"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 60)
# pB_res <- ggplot(calc_resid(allData_newV[var=="B"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='B - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pD <- ggplot(allData_newV[var=="D"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("D - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="D"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 75)
# pD_res <- ggplot(calc_resid(allData_newV[var=="D"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='D - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pH <- ggplot(allData_newV[var=="H"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("H - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="H"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 45)
# pH_res <- ggplot(calc_resid(allData_newV[var=="H"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='H - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pHc <- ggplot(allData_newV[var=="Hc"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Hc - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="Hc"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 25)
# pHc_res <- ggplot(calc_resid(allData_newV[var=="Hc"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   labs(title='Hc - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf1 <- ggplot(allData_newV[var=="Wf1"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf1 - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="Wf1"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 40)
# pWf1_res <- ggplot(calc_resid(allData_newV[var=="Wf1"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf1 - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf2 <- ggplot(allData_newV[var=="Wf2"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf2 - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="Wf2"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 25)
# pWf2_res <- ggplot(calc_resid(allData_newV[var=="Wf2"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf2 - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pAs <- ggplot(allData_newV[var=="As"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("As - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV[var=="As"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 0.1)
# pAs_res <- ggplot(calc_resid(allData_newV[var=="As"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='As - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pGPP <- ggplot(allData_newV[var=="GPP"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("GPP - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newV[var=="GPP"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pGPP_res <- ggplot(calc_residF(allData_newV[var=="GPP"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pET <- ggplot(allData_newV[var=="ET"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("ET - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newV[var=="ET"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pET_res <- ggplot(calc_residF(allData_newV[var=="ET"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='ET - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# dt<-allData_newV[var=="GPP_yr"]
# dt$year<-c(1,2,3,1,2,3,4,5,6)
# pGPPyr <- ggplot(dt, aes(x= year)) + geom_point(aes(y=obs, col=siteID)) + 
#   geom_line(aes(y=sim, col=siteID)) +ylim(0,1500) + ggtitle("GPP annual - new version\nactual pValues\n")+
#   theme(plot.title = element_text(size=10))
# pGPPyr_res <- ggplot(calc_residF(allData_newV[var=="GPP_yr"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP annual - new version:\nactual pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# ### NEW VERSION - NEW PARAMETERS -- OBSERVED VS. PREDICTED
# pV_newVP <- ggplot(allData_newV_Par[var=="V"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("V - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="V"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 800)
# pV_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="V"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='V - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pB_newVP <- ggplot(allData_newV_Par[var=="B"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("B - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="B"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 60)
# pB_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="B"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='B - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pD_newVP <- ggplot(allData_newV_Par[var=="D"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("D - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="D"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 75)
# pD_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="D"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='D - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pH_newVP <- ggplot(allData_newV_Par[var=="H"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("H - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="H"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 45)
# pH_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="H"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='H - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pHc_newVP <- ggplot(allData_newV_Par[var=="Hc"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Hc - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="Hc"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4,6), size=3) + ylim(0, 25)
# pHc_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="Hc"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Hc - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf1_newVP <- ggplot(allData_newV_Par[var=="Wf1"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf1 - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="Wf1"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 40)
# pWf1_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="Wf1"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf1 - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pWf2_newVP <- ggplot(allData_newV_Par[var=="Wf2"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("Wf2 - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="Wf2"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 25)
# pWf2_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="Wf2"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='Wf2 - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pAs_newVP <- ggplot(allData_newV_Par[var=="As"],aes(x=obs,y=sim,col=speciesID)) +
#   geom_point() + geom_abline() + ggtitle("As - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse(allData_newV_Par[var=="As"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3) + ylim(0, 0.1)
# pAs_res_newVP <- ggplot(calc_resid(allData_newV_Par[var=="As"]), aes(x = .fitted, y = .resid, col=speciesID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='As - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pGPP_newVP <- ggplot(allData_newV_Par[var=="GPP"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("GPP - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newV_Par[var=="GPP"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pGPP_res_newVP <- ggplot(calc_residF(allData_newV_Par[var=="GPP"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# pET_newVP <- ggplot(allData_newV_Par[var=="ET"],aes(x=obs,y=sim,col=siteID)) +
#   geom_point() + geom_abline() + ggtitle("ET - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10)) +
#   geom_text(data=r2_rmse_flux(allData_newV_Par[var=="ET"]), aes(x=-Inf,y=+Inf,label=label),
#             hjust = 0, vjust = c(2,4), size=3)
# pET_res_newVP <- ggplot(calc_residF(allData_newV_Par[var=="ET"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='ET - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# dt_newVP<-allData_newV_Par[var=="GPP_yr"]
# dt_newVP$year<-c(1,2,3,1,2,3,4,5,6)
# pGPPyr_newVP <- ggplot(dt_newVP, aes(x= year)) + geom_point(aes(y=obs, col=siteID)) + 
#   geom_line(aes(y=sim, col=siteID)) +ylim(0,1500) + ggtitle("GPP annual - new version\nsuggested pValues\n")+
#   theme(plot.title = element_text(size=10))
# pGPPyr_res_newVP <- ggplot(calc_residF(allData_newV_Par[var=="GPP_yr"]), aes(x = .fitted, y = .resid, col=siteID)) +
#   geom_point() + geom_hline(yintercept = 0) +
#   theme(plot.title = element_text(size=10)) +
#   labs(title='GPP annual - new version:\nsuggested pValues\nResidual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
# 
# ggarrange(pV_o,pV_newP,pV,pV_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pB_o,pB_newP,pB,pB_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pD_o,pD_newP,pD,pD_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pH_o,pH_newP,pH,pH_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pHc_o,pHc_newP,pHc,pHc_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pWf1_o,pWf1_newP,pWf1,pWf1_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pWf2_o,pWf2_newP,pWf2,pWf2_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pAs_o,pAs_newP,pAs,pAs_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pGPPyr_o,pGPPyr_newP,pGPPyr,pGPPyr_newVP,  ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pGPP_o,pGPP_newP,pGPP,pGPP_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# ggarrange(pET_o,pET_newP,pET,pET_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# # ggarrange(pV_res_o,pV_res_newP,pV_res,pV_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pB_res_o,pB_res_newP,pB_res,pB_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pD_res_o,pD_res_newP,pD_res,pD_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pH_res_o,pH_res_newP,pH_res,pH_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_res_o,pHc_res_newP,pHc_res,pHc_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_res_o,pWf1_res_newP,pWf1_res,pWf1_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf2_res_o,pWf2_res_newP,pWf2_res,pWf2_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pAs_res_o,pAs_res_newP,pAs_res,pAs_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_res_o,pGPPyr_res_newP,pGPPyr_res,pGPPyr_res_newVP,  ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPP_res_o,pGPP_res_newP,pGPP_res,pGPP_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pET_res_o,pET_res_newP,pET_res,pET_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# 
# # ggarrange(pV_o,pB_o,pD_o,pH_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_o,pWf2_o,pAs_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_o,pGPP_o,pET_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # 
# # ggarrange(pV_res_o,pB_res_o,pD_res_o,pH_res_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_res_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_res_o,pWf2_res_o,pAs_res_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_res_o,pGPP_res_o,pET_res_o, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# # ggarrange(pV_newP,pB_newP,pD_newP,pH_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_newP,pWf2_newP,pAs_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_newP,pGPP_newP,pET_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # 
# # ggarrange(pV_res_newP,pB_res_newP,pD_res_newP,pH_res_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_res_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_res_newP,pWf2_res_newP,pAs_res_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_res_newP,pGPP_res_newP,pET_res_newP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# # ggarrange(pV,pB,pD,pH, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1,pWf2,pAs, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr,pGPP,pET, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # 
# # ggarrange(pV_res,pB_res,pD_res,pH_res, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_res, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_res,pWf2_res,pAs_res, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_res,pGPP_res,pET_res, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# # ggarrange(pV_newVP,pB_newVP,pD_newVP,pH_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_newVP,pWf2_newVP,pAs_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_newVP,pGPP_newVP,pET_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # 
# # ggarrange(pV_res_newVP,pB_res_newVP,pD_res_newVP,pH_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pHc_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pWf1_res_newVP,pWf2_res_newVP,pAs_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# # ggarrange(pGPPyr_res_newVP,pGPP_res_newVP,pET_res_newVP, ncol = 2, nrow = 2, legend = "right", common.legend = T)
# 
# 
# # print(pV)
# # print(pB)
# # print(pD)
# # print(pH)
# # print(pHc)
# # print(pGPP)
# #print(pET)
# #dev.off()
# # 
# # #pdf(file="out/outupdate_residuals.pdf")
# # print(pV_res)
# # print(pB_res)
# # print(pD_res)
# # print(pH_res)
# # print(pHc_res)
# # print(pGPP_res)
# # print(pET_res)
# dev.off()
# 
# 
# ### error decomposition
# 
# resultV1d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='pine')],method = 1)
# resultH1d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='pine')],method = 1)
# resultD1d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='pine')],method = 1)
# resultB1d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='pine')],method = 1)
# resultHc1d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='pine')],method = 1)
# resultWf11d1<- MSEdec("Wf1",allData$obs[which(allData$var=="Wf1" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="Wf1" & allData$speciesID=='pine')],method = 1)
# resultWf21d1<- MSEdec("Wf2",allData$obs[which(allData$var=="Wf2" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="Wf2" & allData$speciesID=='pine')],method = 1)
# resultAs1d1<- MSEdec("As",allData$obs[which(allData$var=="As" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="As" & allData$speciesID=='pine')],method = 1)
# 
# resultV2d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='spruce')],method = 1)
# resultH2d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='spruce')],method = 1)
# resultD2d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='spruce')],method = 1)
# resultB2d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='spruce')],method = 1)
# resultHc2d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='spruce')],method = 1)
# resultWf12d1<- MSEdec("Wf1",allData$obs[which(allData$var=="Wf1" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="Wf1" & allData$speciesID=='spruce')],method = 1)
# resultWf22d1<- MSEdec("Wf2",allData$obs[which(allData$var=="Wf2" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="Wf2" & allData$speciesID=='spruce')],method = 1)
# resultAs2d1<- MSEdec("As",allData$obs[which(allData$var=="As" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="As" & allData$speciesID=='spruce')],method = 1)
# 
# resultV3d1<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='birch')],method = 1)
# resultH3d1<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='birch')],method = 1)
# resultD3d1<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='birch')],method = 1)
# resultB3d1<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='birch')],method = 1)
# resultHc3d1<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='birch')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='birch')],method = 1)
# 
# resultGd1<- MSEdec("GPP",allData$obs[which(allData$var=="GPP")],
#                    allData$sim[which(allData$var=="GPP")],method = 1)
# resultEd1<- MSEdec("ET",allData$obs[which(allData$var=="ET")],
#                    allData$sim[which(allData$var=="ET")],method = 1)
# method_Kob <- Map(c,resultV1d1,resultH1d1,resultD1d1,resultB1d1,resultHc1d1,resultWf11d1,resultWf21d1,resultAs1d1,
#                   resultV2d1,resultH2d1,resultD2d1,resultB2d1,resultHc2d1,resultWf12d1,resultWf22d1,resultAs2d1,
#                   resultV3d1,resultH3d1,resultD3d1,resultB3d1,resultHc3d1,
#                   resultGd1,resultEd1)
# 
# method_Kob <- data.frame(matrix(unlist(method_Kob), nrow=length(method_Kob), byrow=T))
# row.names(method_Kob)<-c("Var","sb","sdsd","lc","mse")
# write.csv(method_Kob,'Z:/PREBAS_calibration/newCal/out/errordecomp_kobsal_newV.csv')
# 
# 
# resultV1d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='pine')],method = 2)
# resultH1d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='pine')],method = 2)
# resultD1d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='pine')],method = 2)
# resultB1d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='pine')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='pine')],method = 2)
# resultHc1d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='pine')],method = 2)
# resultWf11d2<- MSEdec("Wf1",allData$obs[which(allData$var=="Wf1" & allData$speciesID=='pine')],
#                       allData$sim[which(allData$var=="Wf1" & allData$speciesID=='pine')],method = 2)
# resultWf21d2<- MSEdec("Wf2",allData$obs[which(allData$var=="Wf2" & allData$speciesID=='pine')],
#                       allData$sim[which(allData$var=="Wf2" & allData$speciesID=='pine')],method = 2)
# resultAs1d2<- MSEdec("As",allData$obs[which(allData$var=="As" & allData$speciesID=='pine')],
#                      allData$sim[which(allData$var=="As" & allData$speciesID=='pine')],method = 2)
# 
# resultV2d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='spruce')],method = 2)
# resultH2d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='spruce')],method = 2)
# resultD2d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='spruce')],method = 2)
# resultB2d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='spruce')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='spruce')],method = 2)
# resultHc2d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='spruce')],method = 2)
# resultWf12d2<- MSEdec("Wf1",allData$obs[which(allData$var=="Wf1" & allData$speciesID=='spruce')],
#                       allData$sim[which(allData$var=="Wf1" & allData$speciesID=='spruce')],method = 2)
# resultWf22d2<- MSEdec("Wf2",allData$obs[which(allData$var=="Wf2" & allData$speciesID=='spruce')],
#                       allData$sim[which(allData$var=="Wf2" & allData$speciesID=='spruce')],method = 2)
# resultAs2d2<- MSEdec("As",allData$obs[which(allData$var=="As" & allData$speciesID=='spruce')],
#                      allData$sim[which(allData$var=="As" & allData$speciesID=='spruce')],method = 2)
# 
# resultV3d2<- MSEdec("V",allData$obs[which(allData$var=="V" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="V" & allData$speciesID=='birch')],method = 2)
# resultH3d2<- MSEdec("H",allData$obs[which(allData$var=="H" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="H" & allData$speciesID=='birch')],method = 2)
# resultD3d2<- MSEdec("D",allData$obs[which(allData$var=="D" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="D" & allData$speciesID=='birch')],method = 2)
# resultB3d2<- MSEdec("B",allData$obs[which(allData$var=="B" & allData$speciesID=='birch')],
#                     allData$sim[which(allData$var=="B" & allData$speciesID=='birch')],method = 2)
# resultHc3d2<- MSEdec("Hc",allData$obs[which(allData$var=="Hc" & allData$speciesID=='birch')],
#                      allData$sim[which(allData$var=="Hc" & allData$speciesID=='birch')],method = 2)
# 
# resultGd2<- MSEdec("GPP",allData$obs[which(allData$var=="GPP")],
#                    allData$sim[which(allData$var=="GPP")],method = 2)
# resultEd2<- MSEdec("ET",allData$obs[which(allData$var=="ET")],
#                    allData$sim[which(allData$var=="ET")],method = 2)
# 
# method_Gauch <- Map(c,resultV1d2,resultH1d2,resultD1d2,resultB1d2,resultHc1d2,resultWf11d2,resultWf21d2,resultAs1d2,
#                     resultV2d2,resultH2d2,resultD2d2,resultB2d2,resultHc2d2,resultWf12d2,resultWf22d2,resultAs2d2,
#                     resultV3d2,resultH3d2,resultD3d2,resultB3d2,resultHc3d2,
#                     resultGd2,resultEd2)
# 
# method_Gauch <- data.frame(matrix(unlist(method_Gauch), nrow=length(method_Gauch), byrow=T))
# row.names(method_Gauch)<-c("Var","sb","sdsd","lc","mse")
# write.csv(method_Gauch,'Z:/PREBAS_calibration/newCal/out/errordecomp_gauch_newV.csv')
# 
