####settings
if(!exists("vLocal")) vLocal <- FALSE

if(!exists("fX")) fX <- 2.38  ###parameter of DEzs algorithm
if(!exists("calAlg")) calAlg <- "DEzs"  ###MC algorithm
if(!exists("nworkers")) nworkers <- 20  ###number of cores

###siteIDs used in the different calibration sets
set1 <- 1:300
set2 <- 301:600
set3 <- 601:922
set4 <- 923:948
set5 <- 949:950

###process parameters
parSel_PREL <- c(5:11,14:18,21,31:34)
nparPREL <- length(parSel_PREL)
#CROBAS parameters
#parSel_CROB <- c(1:18,31:34,38,41) #original
#parSel_CROB <- c(1:15,17,18,38,41) #modified 
parSel_CROB <- c(1:15,17,18,21:25,38,41) #added alphar
nparCROB <- length(parSel_CROB)
param_all <- read.csv('inputs/par_prebas_newCal.csv')
parmod <- param_all[,2] 
parmin <- param_all[,3] 
parmax <- param_all[,4] 
parnam <- as.character(param_all[,1])
npar <- length(parmod)

###error parameters
a_GPPind <- which(parnam=="GPP_NT_a")
b_GPPind <- which(parnam=="GPP_NT_b")
a_ETind <- which(parnam=="ET_a")
b_ETind <- which(parnam=="ET_b")
a_Hind <- which(parnam=="a_H")
b_Hind <- which(parnam=="b_H")
a_Bind <- which(parnam=="a_B")
b_Bind <- which(parnam=="b_B")
a_Dind <- which(parnam=="a_D")
b_Dind <- which(parnam=="b_D")
a_Vind <- which(parnam=="a_V")
b_Vind <- which(parnam=="b_V")
a_Hcind <- which(parnam=="a_Hc")
b_Hcind <- which(parnam=="b_Hc")
a_Wf1ind <- which(parnam=="a_Wf1")
b_Wf1ind <- which(parnam=="b_Wf1")
a_Wf2ind <- which(parnam=="a_Wf2")
b_Wf2ind <- which(parnam=="b_Wf2")
a_Acind <- which(parnam=="a_Ac")
b_Acind <- which(parnam=="b_Ac")
a_WfDataind <- which(parnam=="a_WfData")
b_WfDataind <- which(parnam=="b_WfData")

nCores <- ifelse(vLocal,1,nworkers)

sets <- 1:5
likelihoods <- list()
likelihoods[[1]] <- likelihood1
likelihoods[[2]] <- likelihood2
likelihoods[[3]] <- likelihood3
likelihoods[[4]] <- likelihood4
likelihoods[[5]] <- likelihood5Flux
