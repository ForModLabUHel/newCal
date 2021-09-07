##I have intergrated the old growth data into initprebas and made it data set4
##questions:how to set pvalues in the likelihood function? a serial of numbers for each parameter


# initPrebas$pCROBAS <- pCROB
# initPrebas$pCROBAS[parSel,1] <- pValues[1:nparCROB]
# initPrebas$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
# initPrebas$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
# ###init_set2
# init_set2$pCROBAS <- pCROB
# init_set2$pCROBAS[parSel,1] <- pValues[1:nparCROB]
# init_set2$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
# init_set2$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
# ###init_set3
# init_set3$pCROBAS <- pCROB
# init_set3$pCROBAS[parSel,1] <- pValues[1:nparCROB]
# init_set3$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
# init_set3$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
# ###init_set4
# init_set4$pCROBAS <- pCROB
# init_set4$pCROBAS[parSel,1] <- pValues[1:nparCROB]
# init_set4$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
# init_set4$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]

# function to subset Data
subSetData <- function(outdata,setX,obs){
  selX <- which(outdata[,1] %in% setX)
  if(length(selX)>0){
    obsSub <- obs[selX]
    outdataSub <- outdata[selX,]
    outdataSub[,1] <- outdataSub[,1] - (min(setX)-1)
    return(list(obs=obsSub,outData= outdataSub))
  }
}

####Function to Subset initPrebas 
subInit <- function(initPrebas,setX){
  nYears <- initPrebas$nYears[setX]
  siteInfo=initPrebas$siteInfo[setX,]
  siteInfo[,1] <- siteInfo[,1] - (min(setX)-1)
  defaultThin=initPrebas$defaultThin[setX]
  ClCut = initPrebas$ClCut[setX]
  multiInitVar = initPrebas$multiInitVar[setX,,]
  yassoRun <- initPrebas$yassoRun[setX]
  multiThin = initPrebas$thinning[setX,,]
  multiThin[,,2][which(multiThin[,,2]>0)] <- 
    multiThin[,,2][which(multiThin[,,2]>0)] - (min(setX)-1)
  multiNthin = initPrebas$nThinning[setX]
  multiThin <- multiThin[,1:max(multiNthin),]
  
  climIDx <- sort(unique(siteInfo[,2]))
  siteInfo[,2] <- match(siteInfo[,2],climIDx)
  maxYears <- max(initPrebas$nYears[setX])
  CO2 <- TAir <- VPD <- Precip <- PAR <- matrix(-999,length(climIDx),maxYears*365)
  
  for(i in 1:length(climIDx)){
    ij <- climIDx[i]
    PAR[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,1]))
    TAir[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,2]))
    Precip[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,4]))
    VPD[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,3]))
    CO2[i,] <- as.vector(t(initPrebas$weather[ij,1:maxYears,,5]))
  }
  
  init_setX <- InitMultiSite(nYearsMS = nYears,
                             siteInfo=siteInfo,
                             # pCROBAS = pCROBAS,
                             # litterSize = litterSize,#pAWEN = parsAWEN,
                             defaultThin=defaultThin,
                             ClCut = ClCut,
                             multiInitVar = multiInitVar,
                             PAR = PAR,
                             TAir= TAir,
                             VPD= VPD,
                             Precip= Precip,
                             CO2= CO2,
                             yassoRun = yassoRun,#lukeRuns = initPrebas$lukeRuns,
                             # initCLcutRatio = initCLcutRatio
                             multiThin = as.array(multiThin),
                             multiNthin = as.vector(multiNthin)
  )
  return(init_setX)
}

likelihood1 <- function(pValues,cal=T){
  ###init_set1
  init_set1$pCROBAS <- pCROB
  init_set1$pCROBAS[parSel,1] <- pValues[1:nparCROB]
  init_set1$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
  init_set1$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set1)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s1$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s1$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s1$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s1$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s1$outData]-Hdata_s1$obs
  diff_Hc <- output[Hcdata_s1$outData]-Hcdata_s1$obs
  diff_D <- output[Ddata_s1$outData]-Ddata_s1$obs
  diff_B <- out_B-Bdata_s1$obs
  diff_V <- out_V-Vdata_s1$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[a_Hind]+pValues[b_Hind]*output[Hdata_s1$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[a_Dind]+pValues[b_Dind]*output[Ddata_s1$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[a_Bind]+pValues[b_Bind]*output[Bdata_s1$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[a_Hcind]+pValues[b_Hcind]*output[Hcdata_s1$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[a_Vind]+pValues[b_Vind]*output[Vdata_s1$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  if(cal==T){
    return(loglikelihood)  
  }else{
    return(list(simV=outV,simB=outB,simH=output[Hdata_s1$outData],
                obsH = Hdata_s1$obs,simHc = output[Hcdata_s1$outData],
                obsHc = Hcdata_s1$obs,simD=output[Ddata_s1$outData],
                obsD=Ddata_s1$obs,obsB =Bdata_s1$obs,
                obsV =Vdata_s1$obs,
                output=output))
  }
  
}


likelihood2 <- function(pValues,cal=T){
  ###init_set2
  init_set2$pCROBAS <- pCROB
  init_set2$pCROBAS[parSel,1] <- pValues[1:nparCROB]
  init_set2$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
  init_set2$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set2)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s2$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s2$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s2$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s2$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s2$outData]-Hdata_s2$obs
  diff_Hc <- output[Hcdata_s2$outData]-Hcdata_s2$obs
  diff_D <- output[Ddata_s2$outData]-Ddata_s2$obs
  diff_B <- out_B-Bdata_s2$obs
  diff_V <- out_V-Vdata_s2$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[a_Hind]+pValues[b_Hind]*output[Hdata_s2$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[a_Dind]+pValues[b_Dind]*output[Ddata_s2$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[a_Bind]+pValues[b_Bind]*output[Bdata_s2$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[a_Hcind]+pValues[b_Hcind]*output[Hcdata_s2$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[a_Vind]+pValues[b_Vind]*output[Vdata_s2$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  
  if(cal==T){
    return(loglikelihood)  
  }else{
    return(list(simV=outV,simB=outB,simH=output[Hdata_s2$outData],
                obsH = Hdata_s2$obs,simHc = output[Hcdata_s2$outData],
                obsHc = Hcdata_s2$obs,simD=output[Ddata_s2$outData],
                obsD=Ddata_s2$obs,obsB =Bdata_s2$obs,
                obsV =Vdata_s2$obs,
                output=output))
  }
}


likelihood3 <- function(pValues,cal=T){
  ###init_set3
  init_set3$pCROBAS <- pCROB
  init_set3$pCROBAS[parSel,1] <- pValues[1:nparCROB]
  init_set3$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
  init_set3$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set3)$multiOut
  # if (output==-999){
  #   loglikelihood= -Inf
  # } else {
  outdata_B2 <- Bdata_s3$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s3$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s3$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s3$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s3$outData]-Hdata_s3$obs
  diff_Hc <- output[Hcdata_s3$outData]-Hcdata_s3$obs
  diff_D <- output[Ddata_s3$outData]-Ddata_s3$obs
  diff_B <- out_B-Bdata_s3$obs
  diff_V <- out_V-Vdata_s3$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[a_Hind]+pValues[b_Hind]*output[Hdata_s3$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[a_Dind]+pValues[b_Dind]*output[Ddata_s3$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[a_Bind]+pValues[b_Bind]*output[Bdata_s3$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[a_Hcind]+pValues[b_Hcind]*output[Hcdata_s3$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[a_Vind]+pValues[b_Vind]*output[Vdata_s3$outData]))
  
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V)
  # }
  
  if(cal==T){
    return(loglikelihood)  
  }else{
    return(list(simV=outV,simB=outB,simH=output[Hdata_s3$outData],
                obsH = Hdata_s3$obs,simHc = output[Hcdata_s3$outData],
                obsHc = Hcdata_s3$obs,simD=output[Ddata_s3$outData],
                obsD=Ddata_s3$obs,obsB =Bdata_s3$obs,
                obsV =Vdata_s3$obs,
                output=output))
  }
}

##########old growth likelihood
likelihood4 <- function(pValues,cal=T){
  ###init_set4
  init_set4$pCROBAS <- pCROB
  init_set4$pCROBAS[parSel,1] <- pValues[1:nparCROB]
  init_set4$pCROBAS[parSel,2] <- pValues[(nparCROB + 1):(nparCROB*2)]
  init_set4$pCROBAS[parSel,3] <- pValues[(nparCROB*2 + 1):(nparCROB*3)]
  
  output <- multiPrebas(init_set4)$multiOut
  
  outdata_B2 <- Bdata_s4$outData; outdata_B2[,5] <- 2
  outdata_V2 <- Vdata_s4$outData; outdata_V2[,5] <- 2
  
  out_B <-  output[Bdata_s4$outData] + output[outdata_B2]
  out_V <-  output[Vdata_s4$outData] + output[outdata_V2]
  diff_H <- output[Hdata_s4$outData] - Hdata_s4$obs
  diff_Hc <- output[Hcdata_s4$outData] - Hcdata_s4$obs
  diff_D <- output[Ddata_s4$outData] - Ddata_s4$obs
  diff_B <- out_B-Bdata_s4$obs
  diff_V <- out_V-Vdata_s4$obs
  
  ##Sivia likelihood
  ll_H <- sum(Sivia_log(diff_H,sd = pValues[a_Hind]+pValues[b_Hind]*output[Hdata_s4$outData]))
  ll_D <- sum(Sivia_log(diff_D,sd = pValues[a_Dind]+pValues[b_Dind]*output[Ddata_s4$outData]))
  ll_B <- sum(Sivia_log(diff_B,sd = pValues[a_Bind]+pValues[b_Bind]*output[Bdata_s4$outData]))
  ll_Hc <- sum(Sivia_log(diff_Hc,sd = pValues[a_Hcind]+pValues[b_Hcind]*output[Hcdata_s4$outData]))
  ll_V <- sum(Sivia_log(diff_V,sd = pValues[a_Vind]+pValues[b_Vind]*output[Vdata_s4$outData]))
  
  ##vapu data
  As_p_obs <- vapu_P$Ac
  wf_p_obs <- vapu_P$Wf.1 ## wf.1 is the carbon, wf is the kg.
  Lc_p <- vapu_P$Hc
  As_s_obs <- vapu_S$Ac
  wf_s_obs <- vapu_S$Wf.1
  Lc_s<-vapu_S$Hc.m

  ksi_p <- init_set4$pCROBAS[38,1]
  rhof_p <- init_set4$pCROBAS[15,1]
  z_p <- init_set4$pCROBAS[11,1]
  Wf1_p <- rhof_p*As_p_obs
  Wf2_p <- ksi_p*Lc_p^z_p
  As_p <- ksi_p/rhof_p*Lc_p^z_p

  diff_wf1_p <- Wf1_p-wf_p_obs
  diff_wf2_p <- Wf2_p-wf_p_obs
  diff_As_p <- As_p-As_p_obs

  Wfsd1 = pValues[a_WfDataind]+pValues[b_WfDataind]*wf_p_obs +pValues[a_Wf1ind]+pValues[b_Wf1ind]*wf_p_obs
  Wfsd2 = pValues[a_WfDataind]+pValues[b_WfDataind]*wf_p_obs +pValues[a_Wf2ind]+pValues[b_Wf2ind]*wf_p_obs
  ll_wf1_p <- sum(Sivia_log(diff_wf1_p,sd = Wfsd1))
  ll_wf2_p <- sum(Sivia_log(diff_wf2_p,sd = Wfsd2))
  ll_As_p <- sum(Sivia_log(diff_As_p,sd = pValues[a_Acind]+pValues[b_Acind]*As_p_obs))
  
  ksi_s <- init_set4$pCROBAS[38,2]
  rhof_s <- init_set4$pCROBAS[15,2]
  z_s <- init_set4$pCROBAS[11,2]
  Wf1_s <- rhof_s*As_s_obs
  Wf2_s <- ksi_s*Lc_s^z_s
  As_s <- ksi_s/rhof_s*Lc_s^z_s

  diff_wf1_s <- Wf1_s-wf_s_obs
  diff_wf2_s <- Wf2_s-wf_s_obs
  diff_As_s <- As_s-As_s_obs
      
  Wfsd1 = pValues[a_WfDataind]+pValues[b_WfDataind]*wf_s_obs + pValues[a_Wf1ind]+pValues[b_Wf1ind]*Wf1_s
  Wfsd2 = pValues[a_WfDataind]+pValues[b_WfDataind]*wf_s_obs + pValues[a_Wf2ind]+pValues[b_Wf2ind]*Wf2_s
  ll_wf1_s <- sum(Sivia_log(diff_wf1_s,sd = Wfsd1))
  ll_wf2_s <- sum(Sivia_log(diff_wf2_s,sd = Wfsd2))
  ll_As_s <- sum(Sivia_log(diff_As_s,sd = pValues[a_Acind]+pValues[b_Acind]*As_s_obs))
             
  ###Normal distribution
  # ll_H <- sum(dnorm(diff_H,sd = pValues[73]+pValues[74]*output[outdata_H],log=T))
  # ll_D <- sum(dnorm(diff_D,sd = pValues[75]+pValues[76]*output[outdata_D],log=T))
  # ll_B <- sum(dnorm(diff_B,sd = pValues[77]+pValues[78]*output[outdata_B],log=T))
  # ll_Hc <- sum(dnorm(diff_Hc,sd = pValues[79]+pValues[80]*output[outdata_Hc],log=T))
  # ll_V <- sum(dnorm(diff_V,sd = pValues[81]+pValues[82]*output[outdata_V],log=T))
  # 
  loglikelihood <-  sum(ll_H,ll_D,ll_B,ll_Hc,ll_V,ll_wf1_p,ll_wf2_p,ll_As_p,ll_wf1_s,ll_wf2_s,ll_As_s)
  
  if(cal==T){
    return(loglikelihood)  
  }else{
    return(list(simV=outV,simB=outB,simH=output[Hdata_s4$outData],
                obsH = Hdata_s4$obs,simHc = output[Hcdata_s4$outData],
                obsHc = Hcdata_s4$obs,simD=output[Ddata_s4$outData],
                obsD=Ddata_s4$obs,obsB =Bdata_s4$obs,
                obsV =Vdata_s4$obs,
                simWf1_p=Wf1_p,wf_p_obs=wf_p_obs,
                simWf2_p= Wf2_p,wf_p_obs=wf_p_obs,
                simAs_p =As_p,As_p_obs=As_p_obs,
                simWf1_s=Wf1_s,wf_s_obs=wf_s_obs,
                simWf2_s= Wf2_s,wf_s_obs=wf_s_obs,
                simAs_s =As_s,As_s_obs=As_s_obs,
                output=output
                ))
  }
}

##Sivia likelihood
Sivia_log<-function(diff,sd){
  # sd[which(sd<=0)]<-11e-6e-6
  diff[which(abs(diff)<=1e-6)]<-1e-6
  # diff[which(abs(diff)<=1e-6)]<-1e-4
  R2<-(diff/sd)^2
  prob<-1/(sd*(pi*2)^0.5)*(1-exp(-R2/2))/R2
  log(prob)
}


likelihood <- function(pValues){
  logLike <- mclapply(sets, function(jx) {
    likelihoods[[jx]](pValues)  
  }, mc.cores = nCores)  
  llP <- sum(unlist(logLike ))
  return(llP)
}

