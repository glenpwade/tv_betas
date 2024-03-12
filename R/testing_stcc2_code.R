##  Test stcc2 ####

library(MTVGARCH)
library(doParallel)
library(knitr)

## -- calc.Gt -- ####
calc.Gt = function(stcc1Obj){
  this <- stcc1Obj
  
  speed <- this$Estimated$pars["speed"]
  loc1 <- this$Estimated$pars["loc1"]
  loc2 <- this$Estimated$pars["loc2"]
  
  st_c <- 0
  if(this$shape == corrshape$single) { st_c <- this@st - loc1 }
  if(this$shape == corrshape$double) { st_c <- (this@st - loc1)*(this@st - loc2) }
  if(this$shape == corrshape$double1loc) { st_c <- (this@st - loc1)^2 }
  
  G <- 0
  if(this$speedopt == corrspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
  if(this$speedopt == corrspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(this@st)*st_c)) }
  if(this$speedopt == corrspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }
  
  return(matrix(G,nrow = this@Tobs,ncol = 1))
}

calc.Pt =   function(stcc1Obj){
  this <- stcc1Obj
  
  vP1 <- vecL(this$Estimated$P1)
  vP2 <- vecL(this$Estimated$P2)
  
  Gt <- calc.Gt(this)
  Pt <- apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2)
  
  if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1) else Pt <- t(Pt)
  
  return(Pt)
} 

calc.Gt2 = function(stcc2Obj){
    this <- stcc2Obj
    
    speed1 <- this$Estimated$pars["speed1"]
    loc11 <- this$Estimated$pars["loc11"]
    loc12 <- this$Estimated$pars["loc12"]
    speed2 <- this$Estimated$pars["speed2"]
    loc21 <- this$Estimated$pars["loc21"]
    loc22 <- this$Estimated$pars["loc22"]
    
    st_c_1 <- 0
    if(this$shape == corrshape$single) { st_c_1 <- this@st - loc11 }
    if(this$shape == corrshape$double) { st_c_1 <- (this@st - loc11)*(this@st - loc12) }
    if(this$shape == corrshape$double1loc) { st_c_1 <- (this@st - loc11)^2 }
    st_c_2 <- 0
    if(this$shape == corrshape$single) { st_c_2 <- this@st - loc21 }
    if(this$shape == corrshape$double) { st_c_2 <- (this@st - loc21)*(this@st - loc22) }
    if(this$shape == corrshape$double1loc) { st_c_2 <- (this@st - loc21)^2 }
    
    G <- matrix(0,nrow = this@Tobs, ncol = 2)
    if(this$speedopt == corrspeedopt$gamma) { G[,1] <- 1/(1+exp(-speed1*st_c_1)) }
    if(this$speedopt == corrspeedopt$gamma_std) { G[,1] <- 1/(1+exp(-speed1/sd(this@st)*st_c_1)) }
    if(this$speedopt == corrspeedopt$eta) { G[,1] <- 1/(1+exp(-exp(speed1)*st_c_1)) }
    if(this$speedopt == corrspeedopt$gamma) { G[,2] <- 1/(1+exp(-speed2*st_c_2)) }
    if(this$speedopt == corrspeedopt$gamma_std) { G[,2] <- 1/(1+exp(-speed2/sd(this@st)*st_c_2)) }
    if(this$speedopt == corrspeedopt$eta) { G[,2] <- 1/(1+exp(-exp(speed2)*st_c_2)) }
    
    return(matrix(G,nrow = this@Tobs,ncol = 2))
}

calc.Pt2 =   function(stcc2Obj){
    this <- stcc2Obj
    # debug
    # this = corrObj
    
    vP1 <- matrix(vecL(this$Estimated$P1),nrow=1)
    vP2 <- matrix(vecL(this$Estimated$P2),nrow=1)
    vP3 <- matrix(vecL(this$Estimated$P3),nrow=1)
    
    G <- calc.Gt2(this) # T x 2
    Pt <- ((1-G[,2])*(1-G[,1]))%*%vP1 + ((1-G[,2])*G[,1])%*%vP2 + G[,2]%*%vP3 # T x N(N-1)/2
    
    if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1) 
    
    return(Pt)
    
}


N = 3
Tobs = 2000
e <- matrix(rnorm(N*Tobs),Tobs,N)
st = 1:Tobs/Tobs
estCtrl = list(verbose=FALSE,calcSE=FALSE)

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
corrObj <- stcc2(NTV)
# Adjust P1 P2 P3:
corrObj$P1 = toeplitz(0.1^seq.int(0, N-1))
corrObj$P2 = toeplitz(0.5^seq.int(0, N-1))
corrObj$P3 = toeplitz(0.9^seq.int(0, N-1))
corrObj$P1 = unVecL(c(0.2,0,0.5))
corrObj$P2 = unVecL(c(0.8,0,0.2))
corrObj$P3 = unVecL(c(0.2,0.5,0.2))

corrObj$pars = c(5,0.33,NA,5,0.67,NA)
names(corrObj$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")

# Fake Estimation:
corrObj$Estimated <- list()
corrObj$Estimated$error = FALSE
corrObj$Estimated$value = 42
corrObj$Estimated$P1 = corrObj$P1
corrObj$Estimated$P2 = corrObj$P2
corrObj$Estimated$P3 = corrObj$P3
corrObj$Estimated$pars = corrObj$pars

if(FALSE){
    # DGP:
    u <- e
    Pt <- calc.Pt2(corrObj)
    for (t in 1:corrObj@Tobs){
        mPt <- unVecL(Pt[t,,drop=FALSE])
        mPt.sqrt <- sqrt_mat1(mPt)
        e[t,] <- t( mPt.sqrt %*% t(u[t,,drop=FALSE]) )
    }
    
    #Save e!
    saveRDS(e,"SomeCorrDATA.RDS")
    
}else e = readRDS("SomeCorrDATA.RDS")

estCtrl = list(verbose=TRUE,calcSE=FALSE)
estCorr = stcc2(NTV)
# Set Starting Pars:
estCorr$P1 = toeplitz(0.15^seq.int(0, N-1))
estCorr$P2 = toeplitz(0.55^seq.int(0, N-1))
estCorr$P3 = toeplitz(0.95^seq.int(0, N-1))
estCorr$P1 = unVecL(c(0.25,0,0.55))
estCorr$P2 = unVecL(c(0.85,0,0.25))
estCorr$P3 = unVecL(c(0.25,0.55,0.25))
estCorr$pars = c(4,0.25,NA,4,0.75,NA)
names(estCorr$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")

## PKG Estimator:  ####
estCorr = estimateSTCC2(estCorr,estCtrl)

plot(estCorr$Estimated$Pt[,1],type='l')
lines(estCorr$Estimated$Pt[,2],type='l',col="grey75")
lines(estCorr$Estimated$Pt[,3],type='l',col="green")




###############################


tmp <- matrix(0,nrow=Tobs,ncol=N*(N-1)/2)
for (t in 150:(Tobs-150)){
  tmp[t,] <- vecL(cor(e[(t-148):t+150,]) )
}
# N=5:
plot(tmp[150:(Tobs-150),1],type='l')
plot(tmp[150:(Tobs-150),5],type='l')
plot(tmp[150:(Tobs-150),8],type='l')
plot(tmp[150:(Tobs-150),10],type='l')

plot(tmp[150:(Tobs-150),2],type='l')
plot(tmp[150:(Tobs-150),3],type='l')
plot(tmp[150:(Tobs-150),9],type='l')

plot(tmp[150:(Tobs-150),3],type='l')
plot(tmp[150:(Tobs-150),7],type='l')

plot(tmp[150:(Tobs-150),4],type='l')

#N=3
plot(tmp[150:(Tobs-150),1],type='l')
plot(tmp[150:(Tobs-150),2],type='l')
plot(tmp[150:(Tobs-150),3],type='l')


