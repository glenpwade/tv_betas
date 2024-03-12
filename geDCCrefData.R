## -- sqrt_mat1 -- ####
setGeneric(name="sqrt_mat1",
           valueClass = "matrix",
           signature = c("m"),
           def = function(m){
             ## Note: this only works for positive semi-definite (or definite) matrices that are diagonalizable
             ## (no normal Jordan forms, etc.)
             
             # Handle the special case where we have a 1 x 1 matrix:
             if(NROW(m) == 1){
               m.sqrt <- matrix(sqrt(m[1,1]),nrow=1,ncol=1)
             } else {
               m.eig <- eigen(m)
               m.sqrt <- m.eig$vectors %*% diag(sqrt(m.eig$values)) %*% solve(m.eig$vectors)
             }
             return(m.sqrt)
           }
)

## -- sqrt_mat2 -- ####
setGeneric(name="sqrt_mat2",
           valueClass = "list",
           signature = c("mat"),
           def = function(mat){
             # Using the Denman-Beavers algorithm:
             maxit <- 50
             stopifnot(nrow(mat) == ncol(mat))
             niter <- 0
             y <- mat
             z <- diag(rep(1,nrow(mat)))
             for (niter in 1:maxit) {
               y.temp <- 0.5*(y+solve(z))
               z <- 0.5*(z+solve(y))
               y <- y.temp
             }
             return(list(sqrt=y,sqrt.inv=z))
           }
)



## -- generateDCCRefData -- ####
generateDCCRefData=function(nr.series,nr.obs,Qbar,a,b)
{
  
  refData <- matrix(NA,nrow = nr.obs, ncol = nr.series)
  
  # Step1: Generate iid Data & create Correlation
  # u = un-correlated data
  u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)
  
  # e = correlated data (defaults to 'u' if there is no correlation object)
  e <- u
  
  # - - - DCC - - -
  # pass in Qbar matrix and "a" alpha & "b" beta for DCC
  # need to create more than just nr.obs, add discard amount then drop the discard part before returning
  discard <- 2000
  discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
  u <- rbind(discardData,u)
  e <- u
  endRow <- discard + nr.obs
  # starting from t=2! (set Qt[1]=Qbar)
  Qt_1 <- Qbar
  for (t in 2:endRow){	
    Qt <- (1-a-b)*Qbar + a*t(e[t-1,,drop=FALSE])%*%e[t-1,,drop=FALSE] + b*Qt_1 # N x N
    #scale Qt by inverse of its sqrt diagonals (from front and back) to make it correlation matrix
    Pt <- diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)%*%Qt%*%diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)  # N x N
    # create DCC correlated data
    Pt.sqrt <- sqrt_mat1(Pt)
    e[t,] <- t( Pt.sqrt %*% t(u[t,,drop=FALSE]) )
    # for the next loop round, store the "previous" Qt
    Qt_1 <- Qt 
  }
  #return e
  startRow <- discard + 1
  refData <- e[(startRow:endRow), ]
  
  #Return:
  refData
  
}

## -- calc.Gt -- ####

 calc.Gt = function(speed,loc,Tobs){

   st = seq(0,1,length.out=Tobs-1)
   st_c <- st - loc
   G <- 1/(1+exp(-exp(speed)*st_c)) 
   
   return(matrix(G,nrow = Tobs,ncol = 1))
 }



## -- generateDCCRefData -- ####
generateDynDCCRefData=function(nr.series,nr.obs,Q1,Q2,speed,loc,a,b)
{
  
  # Step1: Generate iid Data & create Correlation
  # u = un-correlated data
  u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)
  
  # e = correlated data (defaults to 'u' if there is no correlation object)
  e <- u
  
  # - - - dynDCC - - -
  # pass in Q1 and Q2 matrices, speed and loc, and "a" alpha & "b" beta for DCC
  # need to create more than just nr.obs, add discard amount then drop the discard part before returning
  discard <- 2000
  discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
  u <- rbind(discardData,u)
  e <- u
  endRow <- discard + nr.obs
  # starting from t=2! (set Qt[1]=Qbar)
  Qt_1 <- Q1
  for (t in 2:endRow){	
    if (t>discard){
      Gt <- .calc.Gt(speed,loc,Tobs) # st=timetrend, loc=0.5, gamma s.t. transition takes place between 0.25 and 0.75 of the sample
      Qbar <- (1-Gt)%*%Q1+Gt%*%Q2
    } else Qbar<-Q1
    Qt <- (1-a-b)*Qbar + a*t(e[t-1,,drop=FALSE])%*%e[t-1,,drop=FALSE] + b*Qt_1 # N x N
    #scale Qt by inverse of its sqrt diagonals (from front and back) to make it correlation matrix
    Pt <- diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)%*%Qt%*%diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)  # N x N
    # create DCC correlated data
    Pt.sqrt <- sqrt_mat1(Pt)
    e[t,] <- t( Pt.sqrt %*% t(u[t,,drop=FALSE]) )
    # for the next loop round, store the "previous" Qt
    Qt_1 <- Qt 
  }
  #return e
  startRow <- discard + 1
  e[(startRow:endRow), ]
  
}

# Make some data DCC:
N = 2
Tobs=1000
rho = 0.5
Qbar = toeplitz(rho^seq.int(0, N-1))
a=0.1
b=0.8


# Make some data dynamic DCC:
N = 2 or 5
Tobs=1000
rho = 0.5 or 0.9
Q1 = diag(N)
Q2 = toeplitz(rho^seq.int(0, N-1))
a= 0.05
b= 0.9
speed= 2.5, 5
loc = 0.5



dccData = generateDCCRefData(N,Tobs,Qbar,a,b)
# Simulation code

# N=2,5,10 series, Tobs=1000 length, R=500 replications

# rho = 0.0, 0.5, 0.9 used for Qbar matrix (constant matrix for DCC intercept)

# a = 0.100 b = 0.80 used for DCC dynamics parameters
# a = 0.050 b = 0.90 used for DCC dynamics parameters
# a = 0.025 b = 0.95 used for DCC dynamics parameters

# generate refData with Qbar = toeplitz(rho^seq.int(0, N-1)) and a and b.

# estimate CCC model

# test CCC vs TVCC using st=(1:Tobs)/Tobs


library(doParallel)
library(knitr)

numCores = 7
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

SimRuns = 500
Tobs=1000
# Select N = 2 , 5 , 10  
  N = 2
# Select rho = 0.5 , 0.9 , 0.0 AND rhoname = 5 , 9 , 0
  rho = 0.5
  rhoname =5
  Qbar = toeplitz(rho^seq.int(0, N-1))
# Select a = 0.1 , 0.05 , 0.025 AND aname = 100 , 050 , 025
  a=0.100
  aname = 100
# Select b = 0.8 , 0.9 , 0.95 AND bname = 800 , 900 , 950
  b=0.800
  bname=800
  
dccTest= foreach(n = 1:SimRuns, .combine = cbind, .inorder = FALSE)%dopar%{
  
  dccData = generateDCCRefData(N,Tobs,Qbar,a,b)
  
  # How do these get created??
  cccObj = ccc(N,ntv)
  cccObj = estimateCCC(e,cccObj,estCtrl)

  st = seq(-0.5,0.5,length.out=Tobs)
  #testOrd=1
  #ccctest <- test.CCCParsim(e,cccObj,st,testOrd)
  testOrd=1
  stcctest1 = test.CCCvSTCC1(e,cccObj,st,testOrd)
  testOrd=2
  stcctest2 = test.CCCvSTCC1(e,cccObj,st,testOrd)
  # returns statistic and pvalue
  
  c(stcctest1,stcctest2)  # stat1 ~ pval ~ stat2 ~ pval
  

}
stopCluster(cl)

saveName = paste0("Pilot_N",N,"_T",Tobs,"_Tplz",rhoname,"_a",aname,"_b",bname,".RDS")
saveRDS(dccTest,file=saveName)
