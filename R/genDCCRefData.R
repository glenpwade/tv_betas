
library(MTVGARCH)
library(doParallel)
library(knitr)

setwd("C:/Source/Repos/tv_betas") 

SimRuns = 1500
Tobs=1000
estCtrl = list(verbose=FALSE,calcSE=FALSE)

## Pilot_N2_T1000_Tplz05_a100_b800.RDS ####
## Pilot_N2_T1000_Tplz05_a050_b900.RDS ####
## Pilot_N5_T1000_Tplz09_a050_b900.RDS ####

## Pilot_N5_T1000_Tplz05_a025_b900.RDS

# Select N = 2 , 5 , 10  
N = 2
# Select rho = 0.5 , 0.9 , 0.0 AND rhoname = 05 , 09 , 00
rho = 0.5
rhoname = "05"
Qbar = toeplitz(rho^seq.int(0, N-1))
# Select a = 0.1 , 0.05 , 0.025 AND aname = 100 , 050 , 025
a=0.100
aname = "100"
# Select b = 0.8 , 0.9 , 0.95 AND bname = 800 , 900 , 950
b=0.800
bname="800"

numCores = 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

timestamp()
dccTest= foreach(n = 1:SimRuns, .combine = rbind, .inorder = FALSE, .packages = "MTVGARCH")%dopar%{

    e <- generateDCCRefData(N,Tobs,Qbar,a,b)
    st = seq(-0.5,0.5,length.out=Tobs)
   
    TV = list()
    MTV = list()
    Names = vector()
    for(i in 1:N){
        TV[[i]] <- tv(st,tvshape$delta0only)
        TV[[i]] <- estimateTV(e[,i],TV[[i]],estCtrl)
        MTV[[i]] = tvgarch(TV[[i]],garchtype$noGarch)
        MTV[[i]] <- estimateTVGARCH(e[,i],MTV[[i]],estCtrl)
        Names <- c(Names,i)
    }
    NTV = ntvgarch(MTV, as.character(Names))
    cccObj = ccc(NTV)
    cccObj = estimateCCC(cccObj,estimationCtrl = estCtrl)
    
    #testOrd=1
    #ccctest <- test.CCCParsim(cccObj,st,testOrd)
    testOrd=1
    stcctest1 = test.CCCvSTCC1(cccObj,st,testOrd)
    testOrd=2
    stcctest2 = test.CCCvSTCC1(cccObj,st,testOrd)
    # returns statistic and pvalue
    
    c(as.integer(n),stcctest1["Statistic",1],stcctest1["P_value",1],stcctest2["Statistic",1],stcctest2["P_value",1])  # stat1 ~ pval ~ stat2 ~ pval
    
    
}
stopCluster(cl)

saveName = paste0("Pilot_N",N,"_T",Tobs,"_Tplz",rhoname,"_a",aname,"_b",bname,".RDS")
saveRDS(dccTest,file=saveName)
timestamp()
# 



## Review the output:   ####

dccTest = readRDS("Pilot_N2_T1000_Tplz05_a100_b800.RDS")
dccTest = readRDS("Pilot_N2_T1000_Tplz05_a050_b900.RDS")
dccTest = readRDS("Pilot_N5_T1000_Tplz09_a050_b900.RDS")
dccTest = readRDS("Pilot_N5_T1000_Tplz05_a025_b900.RDS")

hist(dccTest[,2],breaks = 25)
sum(as.integer( dccTest[,3]<0.01) )/nrow(dccTest)
sum(as.integer( dccTest[,3]<0.05) )/nrow(dccTest)
sum(as.integer( dccTest[,3]<0.10) )/nrow(dccTest)
summary(dccTest[,3])


##  Gen CCC Data - Unit Test ####

# Generate Noise Data:
e <- matrix(rnorm(Tobs * N),nrow=Tobs, ncol=N)
corrObj <- ccc(N)
# - - - CCC - - -
if (isTRUE(class(corrObj) == "ccc_class")){
    if(is.null(corrObj$Estimated)) {
        P <- corrObj$P
    } else P <- corrObj$Estimated$P
    P.sqrt <-  sqrt_mat1(P)
    e <- t(P.sqrt %*% t(u))
}
#End: Generate Correlated Data















