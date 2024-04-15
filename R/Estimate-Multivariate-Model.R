#  Initialise  ----

# Refactor this file (and associated saved objects) to work with latest Pkg.  Ver 0.9.4.53 ####

rm(list=ls())
library(MTVGARCH)  # ver. 0.9.4.6
library(knitr)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

setwd("C:/Source/Repos/tv_betas")
setwd("G:/My Drive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project8_tvBetas_mean")   #Anna's Laptop - GOOGLE DRIVE
setwd("G:/My Drive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project8_tvBetas_mean")   #Anna's Surface book Laptop - GOOGLE DRIVE
setwd("C:/Users/silvenno/Google Drive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project8_tvBetas_mean") # Anna's work PC


# Get Data  ----
# Get all the final estimated univariate models & compile them into a multivariate ntvgarch Object

tvg_list = list()
tvg_list[[1]] = readRDS("R/Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("R/Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("R/Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("R/Results/WBC_Final_TVG_model.RDS")
tvg_list[[5]] = readRDS("R/Results/FourB_Final_TVG_model.RDS")
tvg_list[[6]] = readRDS("R/Results/STW_Final_TVG_model.RDS")
tvg_list[[7]] = readRDS("R/Results/SPY_Final_TVG_model.RDS")
tvg_list[[8]] = readRDS("R/Results/PR_Final_TVG_model.RDS")
tvg_list[[9]] = readRDS("R/Results/CRAU_Final_TVG_model.RDS")
tvg_list[[10]] = readRDS("R/Results/CRUS_Final_TVG_model.RDS")

# Make a multivariate model object ----

tvg_names = c("ANZ","CBA","NAB","WBC","FourB","STW","SPY","PR","CRAU","CRUS")
tvgModels <- ntvgarch(tvg_list,tvg_names)

#saveRDS(tvgModels,"R/Results/multivar-spec.RDS")

#tvgModels = readRDS("R/Results/multivar-spec.RDS")

# Refactor: tvgModels now has $e ####
e <- tvgModels$e
colnames(e)<- tvg_names

# Refactor: tvgModels$N has changed to tvgModels$N  ####

# all returns into e
e1 <- matrix()
for (n in 1:tvgModels$N){
  if(n==1) e1 <- tvgModels[[1]]@e
  else e1 = cbind(e1,tvgModels[[n]]@e)
}
# Or reading from Returns_USlagged_4B.RDS, col1=date
if (FALSE){
  e1<-readRDS("R/data/Returns_USlagged_4B.RDS")
}
colnames(e1)<- tvg_names

# Refactor - Proof  ####
identical(e1,e)


# Bivariate CCC: ####
## Run Tests ----
Tobs = NROW(e)
st = seq(0,1,length.out=Tobs)
cccTestMat1 <- matrix(0,nrow=10,ncol=10)
cccTestMat2 <- matrix(0,nrow=10,ncol=10)
for (i in 1:(tvgModels$N-1)){
  for (j in (i+1):tvgModels$N){
    tvgAssetNames <- c(tvg_names[i],tvg_names[j])
    tvgAssetList <- list(tvg_list[[i]],tvg_list[[j]])
    tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
    NullHyp = ccc(tvgAssets)
    NullHyp = estimateCCC(NullHyp)
    testOrd = 1
    ccctest <- test.CCCvSTCC1(NullHyp,st,testOrd)
    ccctest <- (round(ccctest[2],5))
    cccTestMat1[j,i]<- ccctest
    testOrd = 2
    ccctest <- test.CCCvSTCC1(NullHyp,st,testOrd)
    ccctest <- (round(ccctest[2],5))
    cccTestMat2[j,i]<- ccctest
  }
}
colnames(cccTestMat1)<-colnames(cccTestMat2)<-rownames(cccTestMat1)<-rownames(cccTestMat2)<-tvg_names

#saveRDS(cccTestMat1,"R/Results/cccTestMat1.RDS")
#saveRDS(cccTestMat2,"R/Results/cccTestMat2.RDS")

cccTestMat1.old = readRDS("R/Results/cccTestMat1.RDS")
cccTestMat2.old = readRDS("R/Results/cccTestMat2.RDS")
# Refactor Proof
identical(cccTestMat1,cccTestMat1.old)
identical(cccTestMat2,cccTestMat2.old)
# Drop Refactor proof objects
rm(e1,cccTestMat1.old,cccTestMat2.old)


## Test Results ----

# testorder=1
round(cccTestMat1,4)
# > round(cccTestMat1,4)
#          ANZ    CBA    NAB    WBC  FourB    STW    SPY     PR CRAU CRUS
# ANZ   0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# CBA   0.0014 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# NAB   0.0874 0.6928 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# WBC   0.0845 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# FourB 0.0001 0.0026 0.0001 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# STW   0.0000 0.0241 0.0084 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# SPY   0.0259 0.0850 0.1104 0.0177 0.0374 0.0000 0.0000 0.0000    0    0
# PR    0.2038 0.0174 0.0144 0.1839 0.0331 0.0000 0.0571 0.0000    0    0
# CRAU  0.2003 0.0026 0.1669 0.2678 0.0386 0.1518 0.0003 0.4708    0    0
# CRUS  0.0008 0.0015 0.0008 0.0126 0.0008 0.0000 0.0000 0.7382    0    0

# time-varying correlations: 1:ANZ 2:CBA 3:NAB 4:WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
# 1:ANZ  & 2:CBA 5:FourB 6:STW 10:CRUS
# 2:CBA  & 4:WBC 5:FourB 9:CRAU 10:CRUS
# 3:NAB  & 4:WBC 5:FourB 6:STW 10:CRUS 
# 4:WBC  & 5:FourB 6:STW 10:CRUS
# 5:FourB& 6:STW 10:CRUS
# 6:STW  & 7:SPY 8:PR 10:CRUS
# 7:SPY  & 9:CRAU 10:CRUS
# 8:PR   & -
# 9:CRAU & 10:CRUS

# testorder=2
round(cccTestMat2,4)

# > round(cccTestMat2,4)
#          ANZ    CBA    NAB    WBC  FourB    STW    SPY     PR CRAU CRUS
# ANZ   0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# CBA   0.0049 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# NAB   0.0012 0.8220 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# WBC   0.0050 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# FourB 0.0000 0.0022 0.0001 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# STW   0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    0    0
# SPY   0.0012 0.0008 0.0003 0.0352 0.0004 0.0000 0.0000 0.0000    0    0
# PR    0.0155 0.0039 0.0010 0.0219 0.0029 0.0000 0.0576 0.0000    0    0
# CRAU  0.0455 0.0005 0.2069 0.0320 0.0087 0.0011 0.0013 0.2972    0    0
# CRUS  0.0012 0.0038 0.0005 0.0168 0.0012 0.0000 0.0000 0.1638    0    0

# time-varying correlations: 1:ANZ 2:CBA 3:NAB 4:WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
# 1:ANZ  & 2:CBA 3:NAB 4:WBC 5:FourB 6:STW 7:SPY 10:CRUS
# 2:CBA  & 4:WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
# 3:NAB  & 4:WBC 5:FourB 6:STW 7:SPY 8:PR 10:CRUS 
# 4:WBC  & 5:FourB 6:STW
# 5:FourB& 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
# 6:STW  & 7:SPY 8:PR 9:CRAU 10:CRUS
# 7:SPY  & 9:CRAU 10:CRUS
# 8:PR   & -
# 9:CRAU & 10:CRUS




## TVCC estimation ----
cccTestMat1=readRDS("R/Results/cccTestMat1.RDS")
cccTestMat2=readRDS("R/Results/cccTestMat2.RDS")

### Extract only Correlated-Pairs Indices ####

# Set the PValues <= 1% as TRUE
cccTest1_TF = cccTestMat1 <= 0.010
# Set the Upper Tri to FALSE
cccTest1_TF[upper.tri(cccTest1_TF)] <- FALSE
# Get the indexes into a 2 column matrix
cccTest1_idx <- which(cccTest1_TF, arr.ind = TRUE)
#
# cccTest1_idx is a 2-col matrix, holding the indices for all correlated pairs
# Note: for convenience we left the self-correlation in this matrix and use "if(i != j)" below

### Estimate stcc1 for each pair

stcc1_biv <- list()
estControl <- list(calcSE = FALSE,verbose = TRUE)
for (n in 1:NROW(cccTest1_idx)){
  i <- cccTest1_idx[n,2]
  j <- cccTest1_idx[n,1]
  if(i != j){
    tvgAssetNames <- c(tvg_names[i],tvg_names[j])
    tvgAssetList <- list(tvg_list[[i]],tvg_list[[j]])
    tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
    stccObj <- stcc1(tvgAssets)
    stccObj$pars["speed"] <- 2.5
    stccObj$optimcontrol$parscale <- c(1,1,21,3)
    stccObj <- estimateSTCC1(stccObj,estControl)
    itemName <- paste0(tvgAssetNames[1], " ", tvgAssetNames[2])
    stcc1_biv[[itemName]] <- stccObj
    cat("\n",NROW(cccTest1_idx)-n, " estimations left to go...\n")
  }
}


# Set the PValues <= 1% as TRUE
cccTest2_TF = cccTestMat2 <= 0.010
# Set the Upper Tri to FALSE
cccTest2_TF[upper.tri(cccTest2_TF)] <- FALSE
# Get the indexes into a 2 column matrix
cccTest2_idx <- which(cccTest2_TF, arr.ind = TRUE)
  
stcc2_biv <- list()
estControl <- list(calcSE = FALSE,verbose = TRUE)
  for (n in 1:NROW(cccTest2_idx)){
    i <- cccTest2_idx[n,2]
    j <- cccTest2_idx[n,1]
    if(i != j){
      tvgAssetNames <- c(tvg_names[i],tvg_names[j])
      tvgAssetList <- list(tvg_list[[i]],tvg_list[[j]])
      tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
      stccObj <- stcc2(tvgAssets)
      stccObj$pars["speed1"] <- 2.5
      stccObj$pars["speed2"] <- 2.5
      stccObj$optimcontrol$parscale <- c(1,1,1,21,3,21,5)
      stccObj <- estimateSTCC2(stccObj,estControl)
      itemName <- paste0(tvgAssetNames[1], " ", tvgAssetNames[2])
      stcc2_biv[[itemName]] <- stccObj
      cat("\n",NROW(cccTest2_idx)-n, " estimations left to go...\n")
    }
  }

# Refactor Proof:
stcc1_biv.old = readRDS("R/Results/stcc1_biv.RDS")
stcc2_biv.old = readRDS("R/Results/stcc2_biv.RDS")
# Refactor Proof
X = list()
for(n in 1:length(stcc1_biv)){
    # Confirm the Pt matrix is very similar:
    #X[n] <- all.equal(stcc1_biv[[n]]$Estimated$Pt,stcc1_biv.old[[n]]$Estimated$Pt)
    X[n] <- all.equal(stcc1_biv[[n]]$Estimated$pars[[1]],stcc1_biv.old[[n]]$Estimated$pars[[1]])
    #X[n] <- all.equal(stcc1_biv[[n]]$Estimated$pars[[2]],stcc1_biv.old[[n]]$Estimated$pars[[2]])
    #X[n] <- all.equal(stcc2_biv.old[[n]]$Estimated$Pt,stcc2_biv[[n]]$Estimated$Pt)    
    X[n] <- stcc1_biv[[n]]$Estimated$value - stcc1_biv.old[[n]]$Estimated$value  # Should be +ve
   
}
## Conclusion - Estimated models are a bit different.  Qs. Which is better?
## Probably need to examine the plots, so save old objects for later
#

saveRDS(stcc1_biv.old,"R/Results/stcc1_biv.old.RDS")
saveRDS(stcc2_biv.old,"R/Results/stcc2_biv.old.RDS")
  
#saveRDS(stcc1_biv,"R/Results/stcc1_biv.RDS")
#saveRDS(stcc2_biv,"R/Results/stcc2_biv.RDS")
stcc1_biv <- readRDS("R/Results/stcc1_biv.RDS")
stcc2_biv <- readRDS("R/Results/stcc2_biv.RDS")

# Drop Refactor proof objects
rm(stcc1_biv.old,stcc2_biv.old)

  
## Plots ----
##> update ----
# testorder=1
for (i in 1:length(stcc1_biv)){
  plot(stcc1_biv[[i]]$Estimated$Pt,type='l',main=names(stcc1_biv[i]))  
# Refactor proof by plots:
  plot(stcc1_biv.old[[i]]$Estimated$Pt,type='l',main=paste0(names(stcc1_biv.old[i]), " old"))  
}
for (i in 1:length(stcc2_biv)){
  plot(stcc2_biv[[i]]$Estimated$Pt,type='l',main=names(stcc2_biv[i]))  
# Refactor proof by plots:
  plot(stcc2_biv.old[[i]]$Estimated$Pt,type='l',main=paste0(names(stcc2_biv.old[i]), " old"))  
}

# Refactor: change tvgModels@Tobs to tvgModels$Tobs ####

# bivariate cor plots
tmp <- matrix(0,nrow=tvgModels$Tobs,ncol=1)
x <- 100
idx <- c(6,8) # select two: 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
for(t in x:(tvgModels$Tobs-x)){
  tmp[t,]<-vecL(cor(e[(t-x+1):(t+x),idx]))
}

plot(tmp[,1],type='l')




# Multivariate TVCC ####

# View Bi-variate Plots: ####
which(names(stcc2_biv) == "CRAU CRUS")
indices <- c(5,6,7,28,30,31,32,33,34)
indices <- c(30)
for (n in seq_along(indices)){
    i = indices[n]
    plot(stcc2_biv[[i]]$Estimated$Pt,type='l',main=names(stcc2_biv[i]))  
}

## model 1: ANZ - STW - PR (Done) ####
#indices <- c(5,29)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[1]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0,-0.2))
stccObj$P2 = unVecL(c(0.8,NA,0))
stccObj$P3 = unVecL(c(0.6,NA,0.3))
stccObj$pars[1:4] = c(4,0.2,2,0.67)
estCtrl = list(verbose=TRUE,calcSE=TRUE)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)


# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)
#
plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-PR")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-PR")

#saveRDS(stccObj,"R/Results/ANZ_STW_PR.RDS")

## model 1: CBA - STW - PR (Done) ####
# indices <- c(10,12,29)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[2],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[2]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,-0.1,-0.2))
stccObj$P2 = unVecL(c(0.8,-0.2,0))
stccObj$P3 = unVecL(c(0.6,NA,0.3))
stccObj$pars[1:4] = c(4,0.25,4,0.65)
estCtrl = list(verbose=TRUE,calcSE=TRUE)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="CBA-PR")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-PR")

#saveRDS(stccObj,"R/Results/CBA_STW_PR.RDS")


## model 1: NAB - STW - PR (Done) ####
# indices <- c(17,19,29)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[3],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[3]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0,-0.15))
stccObj$P2 = unVecL(c(0.85,-0.3,0))
stccObj$P3 = unVecL(c(0.7,-0.15,0.3))
stccObj$pars[1:4] = c(3,0.2,4,0.75)
estCtrl = list(verbose=TRUE,calcSE=TRUE)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="NAB-PR")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-PR")

#saveRDS(stccObj,"R/Results/NAB_STW_PR.RDS")


## model 1: WBC - STW - PR (Done) ####
# indices <- c(22,29)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[4],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[4]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,-0.1,-0.2))
stccObj$P2 = unVecL(c(0.8,NA,0))
stccObj$P3 = unVecL(c(0.6,NA,0.3))
stccObj$pars[1:4] = c(4,0.2,4,0.75)
estCtrl = list(verbose=TRUE,calcSE=TRUE)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)
plot(stccObj$Estimated$Pt[,1],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="WBC-PR")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-PR")

#saveRDS(stccObj,"R/Results/WBC_STW_PR.RDS")


## (skip) model 1: ANZ~CBA~NAB~WBC - STW - PR ####

# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[2],tvg_names[3],tvg_names[4],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[1]],tvg_list[[2]],tvg_list[[3]],tvg_list[[4]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0.6,0.8,0.6,0.1,   0.8,0.6,0.6,0.1,  0.7,0.6,0,    0.6,-0.1, -0.2))
stccObj$P2 = unVecL(c(0.8,0.8,NA,0.8,-0.2,   NA,0.8,0.8,-0.2,  0.9,0.8,-0.2, 0.8,NA,   0))
stccObj$P3 = unVecL(c(0.6,0.6,NA,0.6,NA,     NA,0.6,0.6,NA,    0.7,0.6,NA,   0.6,NA,   0.3))
stccObj$pars = c(2.5,0.25,NA,2.5,0.75,NA)
estCtrl = list(verbose=TRUE,calcSE=TRUE)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-CBA")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-NAB")
plot(stccObj$Estimated$Pt[,3],type='l',main="ANZ-WBC")
plot(stccObj$Estimated$Pt[,4],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,5],type='l',main="ANZ-PR")
plot(stccObj$Estimated$Pt[,6],type='l',main="CBA-NAB")
plot(stccObj$Estimated$Pt[,7],type='l',main="CBA-WBC")
plot(stccObj$Estimated$Pt[,8],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,9],type='l',main="CBA-PR")
plot(stccObj$Estimated$Pt[,10],type='l',main="NAB-WBC")
plot(stccObj$Estimated$Pt[,11],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,12],type='l',main="NAB-PR")
plot(stccObj$Estimated$Pt[,13],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,14],type='l',main="WBC-PR")
plot(stccObj$Estimated$Pt[,15],type='l',main="STW-PR")


#saveRDS(stccObj,"R/Results/ANZ_CBA_NAB_WBC_STW_PR.RDS")


## model 1: FourB portfolio - STW - PR (Done-ish) ####
#indices <- c(23,25,29)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[5],tvg_names[6],tvg_names[8]) 
tvgAssetList <- list(tvg_list[[5]],tvg_list[[6]],tvg_list[[8]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.7,0.2,-0.1))
stccObj$P2 = unVecL(c(0.9,-0.2,0))
stccObj$P3 = unVecL(c(0.75,-0.1,0.2))
stccObj$pars[1:4] = c(4,0.2,4,0.8)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)
stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="FourB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="FourB-PR")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-PR")

#saveRDS(stccObj,"R/Results/FourB_STW_PR_2.RDS")


## model 2: ANZ - STW - CRAU (Done) ####
# indices <- c(5,30)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[6],tvg_names[9])
tvgAssetList <- list(tvg_list[[1]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.7,0.15,-0.05))
stccObj$P2 = unVecL(c(0.85,NA,0.2))
stccObj$P3 = unVecL(c(0.65,NA,-0.05))
stccObj$pars[1:4] = c(5,0.2,3,0.7)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,21,5)
#stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-CRAU")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/ANZ_STW_CRAU.RDS")


## model 2: CBA - STW - CRAU (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[2],tvg_names[6],tvg_names[9])
tvgAssetList <- list(tvg_list[[2]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.65,0.3,-0.05))
stccObj$P2 = unVecL(c(0.85,-0.05,0.15))
stccObj$P3 = unVecL(c(0.7,0.05,0.05))
stccObj$pars[1:4] = c(4,0.2,3,0.6)
stccObj$optimcontrol$reltol = 1e-05
stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,18,5)
#stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="CBA-CRAU")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/CBA_STW_CRAU.RDS")


## model 2: NAB - STW - CRAU (Done) ####
# indices <- c(17,30)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[3],tvg_names[6],tvg_names[9])
tvgAssetList <- list(tvg_list[[3]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0.1,-0.1))
stccObj$P2 = unVecL(c(0.85,NA,0.2))
stccObj$P3 = unVecL(c(0.7,NA,0.05))
stccObj$pars[1:4] = c(2.5,0.3,2,0.6)
stccObj$optimcontrol$reltol = 1e-05
#stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,18,5)
#stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="NAB-CRAU")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/NAB_STW_CRAU.RDS")


## model 2: WBC - STW - CRAU (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[4],tvg_names[6],tvg_names[9])
tvgAssetList <- list(tvg_list[[4]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0.15,-0.1))
stccObj$P2 = unVecL(c(0.8,NA,0.2))
stccObj$P3 = unVecL(c(0.6,NA,-0.05))
stccObj$pars[1:4] = c(3.5,0.2,2,0.6)
stccObj$optimcontrol$reltol = 1e-05
#stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,18,5)
#stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="WBC-CRAU")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/WBC_STW_CRAU.RDS")


## (skip) model 2: ANZ~CBA~NAB~WBC - STW - CRAU ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[2],tvg_names[3],tvg_names[4],tvg_names[6],tvg_names[9]) 
tvgAssetList <- list(tvg_list[[1]],tvg_list[[2]],tvg_list[[3]],tvg_list[[4]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0.6,0.8,0.6,0.15,   0.8,0.6,0.6,0.4,  0.7,0.6,0.15,    0.6,0.15, -0.1))
stccObj$P2 = unVecL(c(0.8,0.8,NA,0.8,NA,      NA,0.8,0.8,0,     0.9,0.8,NA,      0.8,NA,   0.2))
stccObj$P3 = unVecL(c(0.6,0.6,NA,0.6,NA,     NA,0.6,0.6,NA,    0.7,0.6,NA,   0.6,NA,   -0.1))
stccObj$pars = c(2.5,0.25,NA,2.5,0.75,NA)
names(stccObj$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")
stccObj$sel.vec2 = as.numeric(!is.na(vecL(stccObj$P2)))
stccObj$sel.vec2 = as.numeric(!is.na(vecL(stccObj$P3)))
estCtrl = list(verbose=TRUE,calcSE=FALSE)
stccObj<-myestimateSTCC2.R(stccObj,estCtrl)
plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-CBA")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-NAB")
plot(stccObj$Estimated$Pt[,3],type='l',main="ANZ-WBC")
plot(stccObj$Estimated$Pt[,4],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,5],type='l',main="ANZ-CRAU")
plot(stccObj$Estimated$Pt[,6],type='l',main="CBA-NAB")
plot(stccObj$Estimated$Pt[,7],type='l',main="CBA-WBC")
plot(stccObj$Estimated$Pt[,8],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,9],type='l',main="CBA-CRAU")
plot(stccObj$Estimated$Pt[,10],type='l',main="NAB-WBC")
plot(stccObj$Estimated$Pt[,11],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,12],type='l',main="NAB-CRAU")
plot(stccObj$Estimated$Pt[,13],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,14],type='l',main="WBC-CRAU")
plot(stccObj$Estimated$Pt[,15],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/ANZ_CBA_NAB_WBC_STW_CRAU.RDS")


## model 2: FourB portfolio - STW - CRAU (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[5],tvg_names[6],tvg_names[9]) 
tvgAssetList <- list(tvg_list[[5]],tvg_list[[6]],tvg_list[[9]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.7,0.3,-0.05))
stccObj$P2 = unVecL(c(0.9,0.0,0.15))
stccObj$P3 = unVecL(c(0.75,0.1,0.05))
stccObj$pars[1:4] = c(4,0.35,4,0.75)
stccObj$optimcontrol$reltol = 1e-05
#stccObj$optimcontrol$parscale <- c(1,1,1,1,1,1,1,1,1,21,3,18,5)
#stccObj$optimcontrol$ndeps = c(rep(1e-05,9),1e-06,1e-06,1e-06,1e-09)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="FourB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="FourB-CRAU")
plot(stccObj$Estimated$Pt[,3],type='l',main="STW-CRAU")

#saveRDS(stccObj,"R/Results/FourB_STW_CRAU.RDS")


## model 3: ANZ - STW - SPY - CRAU - CRUS (Done) ####
# indices <- c(5,6,7,28,30,31,32,33,34)
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[6],tvg_names[7],tvg_names[10])
tvgAssetList <- list(tvg_list[[1]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)

stccObj$P1 = unVecL(c(0.6,0.25,  0.2,0.35, 0.2, 0.5))
stccObj$P2 = unVecL(c(0.8,  NA, -0.2,  NA,-0.2, 0.1))
stccObj$P3 = unVecL(c(0.6,  NA,  0.0,  NA, 0.0,-0.2))

stccObj$pars[1:4] = c(4,0.35,4,0.75)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-SPY")
plot(stccObj$Estimated$Pt[,3],type='l',main="ANZ-CRUS")
plot(stccObj$Estimated$Pt[,4],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,5],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,6],type='l',main="SPY-CRUS")

#saveRDS(stccObj,"R/Results/ANZ_STW_SPY_CRUS.RDS")


## model 3: CBA - STW - SPY - CRAU - CRUS (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[2],tvg_names[6],tvg_names[7],tvg_names[10])
tvgAssetList <- list(tvg_list[[2]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)

stccObj$P1 = unVecL(c(0.6, 0.2, 0.2,0.35, 0.2, 0.5))
stccObj$P2 = unVecL(c(0.85, NA, 0.0,  NA,-0.2, 0.1))
stccObj$P3 = unVecL(c(0.7,  NA,-0.1,  NA, 0.0,-0.2))

stccObj$pars[1:4] = c(4,0.15,4,0.8)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="CBA-SPY")
plot(stccObj$Estimated$Pt[,3],type='l',main="CBA-CRUS")
plot(stccObj$Estimated$Pt[,4],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,5],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,6],type='l',main="SPY-CRUS")

#saveRDS(stccObj,"R/Results/CBA_STW_SPY_CRUS.RDS")


## model 3: NAB - STW - SPY - CRAU - CRUS (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[3],tvg_names[6],tvg_names[7],tvg_names[10])
tvgAssetList <- list(tvg_list[[3]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6, 0.27, 0.2,0.15, 0.2, 0.5))
stccObj$P2 = unVecL(c(0.85,  NA, 0.0,  NA, 0.0, 0.1))
stccObj$P3 = unVecL(c(0.7,   NA, 0.0,  NA,  NA,-0.1))

stccObj$pars[1:4] = c(3,0.2,3,0.6)
estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="NAB-SPY")
plot(stccObj$Estimated$Pt[,3],type='l',main="NAB-CRUS")
plot(stccObj$Estimated$Pt[,4],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,5],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,6],type='l',main="SPY-CRUS")

#saveRDS(stccObj,"R/Results/NAB_STW_SPY_CRUS.RDS")


## model 3: WBC - STW - SPY - CRAU - CRUS (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[4],tvg_names[6],tvg_names[7],tvg_names[10])
tvgAssetList <- list(tvg_list[[4]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6, 0.22,0.2,0.35,0.2, 0.5))
stccObj$P2 = unVecL(c(0.85,  NA,0.0,  NA, NA, 0.1))
stccObj$P3 = unVecL(c(0.6,   NA, NA,  NA,0.0,-0.1))
stccObj$pars[1:4] = c(3,0.2,3,0.6)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="WBC-SPY")
plot(stccObj$Estimated$Pt[,3],type='l',main="WBC-CRUS")
plot(stccObj$Estimated$Pt[,4],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,5],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,6],type='l',main="SPY-CRUS")

#saveRDS(stccObj,"R/Results/WBC_STW_SPY_CRUS.RDS")



## (skip) model 3: ANZ~CBA~NAB~WBC - STW - SPY - CRUS ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[1],tvg_names[2],tvg_names[3],tvg_names[4],tvg_names[6],tvg_names[7],tvg_names[10]) 
tvgAssetList <- list(tvg_list[[1]],tvg_list[[2]],tvg_list[[3]],tvg_list[[4]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)
stccObj$P1 = unVecL(c(0.6,0.6,0.8,0.6,0.25,0.2,   0.8,0.6,0.6,0.2,0.2,  0.7,0.6,0.27,0.15,    0.6,0.25,0.2, 0.35,0.2,  0.5))
stccObj$P2 = unVecL(c(0.8,0.8,NA,0.8,NA,-0.2,      NA,0.8,0.8,NA,0,     0.9,0.8,NA,0,         0.8,NA,0,   NA,-0.2,  0.1))
stccObj$P3 = unVecL(c(0.6,0.6,NA,0.6,NA,0,         NA,0.6,0.6,NA,0,     0.7,0.6,NA,NA,        0.6,NA,NA,  NA,0,   -0.2))
stccObj$pars = c(2.5,0.25,NA,2.5,0.75,NA)
names(stccObj$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")
stccObj$sel.vec2 = as.numeric(!is.na(vecL(stccObj$P2)))
stccObj$sel.vec2 = as.numeric(!is.na(vecL(stccObj$P3)))
estCtrl = list(verbose=TRUE,calcSE=FALSE)
stccObj<-myestimateSTCC2.R(stccObj,estCtrl)
plot(stccObj$Estimated$Pt[,1],type='l',main="ANZ-CBA")
plot(stccObj$Estimated$Pt[,2],type='l',main="ANZ-NAB")
plot(stccObj$Estimated$Pt[,3],type='l',main="ANZ-WBC")
plot(stccObj$Estimated$Pt[,4],type='l',main="ANZ-STW")
plot(stccObj$Estimated$Pt[,5],type='l',main="ANZ-SPY")
plot(stccObj$Estimated$Pt[,6],type='l',main="ANZ-CRUS")
plot(stccObj$Estimated$Pt[,7],type='l',main="CBA-NAB")
plot(stccObj$Estimated$Pt[,8],type='l',main="CBA-WBC")
plot(stccObj$Estimated$Pt[,9],type='l',main="CBA-STW")
plot(stccObj$Estimated$Pt[,10],type='l',main="CBA-SPY")
plot(stccObj$Estimated$Pt[,11],type='l',main="CBA-CRUS")
plot(stccObj$Estimated$Pt[,12],type='l',main="NAB-WBC")
plot(stccObj$Estimated$Pt[,13],type='l',main="NAB-STW")
plot(stccObj$Estimated$Pt[,14],type='l',main="NAB-SPY")
plot(stccObj$Estimated$Pt[,15],type='l',main="NAB-CRUS")
plot(stccObj$Estimated$Pt[,16],type='l',main="WBC-STW")
plot(stccObj$Estimated$Pt[,17],type='l',main="WBC-SPY")
plot(stccObj$Estimated$Pt[,18],type='l',main="WBC-CRUS")
plot(stccObj$Estimated$Pt[,19],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,20],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,21],type='l',main="SPY-CRUS")

#saveRDS(stccObj,"R/Results/ANZ_CBA_NAB_WBC_STW_SPY_CRUS.RDS")


## model 3: FourB - STW - SPY - CRUS (Done) ####
# 1:ANZ 2: CBA 3:NAB 4: WBC 5:FourB 6:STW 7:SPY 8:PR 9:CRAU 10:CRUS
tvgAssetNames <- c(tvg_names[5],tvg_names[6],tvg_names[7],tvg_names[10])
tvgAssetList <- list(tvg_list[[5]],tvg_list[[6]],tvg_list[[7]],tvg_list[[10]])
tvgAssets <- ntvgarch(tvgAssetList,tvgAssetNames)
stccObj <- stcc2(tvgAssets)

stccObj$P1 = unVecL(c(0.6, 0.5, 0.2, 0.5, 0.3, 0.6))
stccObj$P2 = unVecL(c(0.85,0  , NA ,   0,  NA,   0))
stccObj$P3 = unVecL(c(0.7, 0.4,-0.1, 0.6,   0,  NA))
stccObj$pars[1:4] = c(4,0.35,4,0.65)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

# Estimate using Pkg:
stccObj <- estimateSTCC2(stccObj,estCtrl)

plot(stccObj$Estimated$Pt[,1],type='l',main="FourB-STW")
plot(stccObj$Estimated$Pt[,2],type='l',main="FourB-SPY")
plot(stccObj$Estimated$Pt[,3],type='l',main="FourB-CRUS")
plot(stccObj$Estimated$Pt[,4],type='l',main="STW-SPY")
plot(stccObj$Estimated$Pt[,5],type='l',main="STW-CRUS")
plot(stccObj$Estimated$Pt[,6],type='l',main="SPY-CRUS")


#saveRDS(stccObj,"R/Results/FourB_STW_SPY_CRUS.RDS")


tmp <- matrix(0,nrow=stccObj@Tobs,ncol=6)
x<-100
for(t in x:(stccObj@Tobs-x)){
  tmp[t,]<-vecL(cor(e[(t-x+1):(t+x),c(10,5,6,9)]))
}
plot(tmp[,1],type='l')
plot(tmp[,2],type='l')
plot(tmp[,3],type='l')
plot(tmp[,4],type='l')
plot(tmp[,5],type='l')
plot(tmp[,6],type='l')





