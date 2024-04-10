##  Test stcc2 ####

library(MTVGARCH)
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

## -- loglik.stcc1.R() --####
loglik.stcc1.R = function(optimpars,z,stcc1Obj){
             
             err_output <- -1e10
             this <- stcc1Obj
             
             this$Estimated$pars <- tail(optimpars,this@nr.trPars)
             tmp.par <- optimpars
             
             # P1 always full of parameters:
             vP1 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P1 <- mP
             
             #Remove the P1 corPars, then extract the P2 corPars
             tmp.par <- tail(tmp.par,-this@nr.corPars)
             nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
             freepars <- tmp.par[1:nr.freepars2]
             vP2<-vP1 # copy of P1
             # replace vP2 elements that have 1's in selvec by free pars
             vP2[which(as.logical(this$sel.vec2))]<-freepars
             mP <- unVecL(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P2 <- mP
             
             #### ======== constraint checks ======== ####
             
             # Check 2: Check the boundary values for speed params:
             speed <- this$Estimated$pars[1]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)
             
             # Check 3: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars[2]
             if(this$shape == corrshape$double) loc2 <- this$Estimated$pars[3] else loc2 <- NA
             if (loc1 < min(this@st)) return(err_output)
             if (loc1 > max(this@st)) return(err_output)
             if (!is.na(loc2)) {
               if (loc2 < min(this@st)) return(err_output)
               if (loc2 > max(this@st)) return(err_output)
             }
             
             
             #### ======== calculate loglikelihood ======== ####
             Pt <- .calc.Pt(this)  # T x nr.corPars
             
             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))
             
}


## -- loglik.stcc2() --####
loglik.stcc2 = function(optimpars,z,stcc2Obj){
             
             err_output <- -1e10
             this <- stcc2Obj
             
             this$Estimated$pars <- tail(optimpars,2*this@nr.trPars)
             tmp.par <- optimpars
             
             vP1 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P1 <- mP
             
             #Remove the P1 corPars, then extract the P2 corPars
             tmp.par <- tail(tmp.par,-this@nr.corPars)
             vP2 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P2 <- mP
             
             #Remove the P2 corPars, then extract the P3 corPars
             tmp.par <- tail(tmp.par,-this@nr.corPars)
             vP3 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP3)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P3 <- mP
             
             #### ======== constraint checks ======== ####
             
             # Check 2.1: Check the boundary values for speed params:
             pos <- 1
             speed <- this$Estimated$pars[pos]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)
             pos<- pos+1
             
             # Check 3.1: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars[pos]
             if(this$shape == corrshape$double){
               pos<-pos+1
               loc2 <- this$Estimated$pars[pos]
             } else {
               loc2 <- NA
             }
             if (loc1 < min(this@st)) return(err_output)
             if (loc1 > max(this@st)) return(err_output)
             if (!is.na(loc2)) {
               if (loc2 < min(this@st)) return(err_output)
               if (loc2 > max(this@st)) return(err_output)
             }
             loc11<-loc1 # keep for later comparison with the second loc1
             pos<-pos+1
             
             # Check 2.2: Check the boundary values for speed params:
             speed <- this$Estimated$pars[pos]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)
             pos<- pos+1
             
             # Check 3.2: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars[pos]
             if(this$shape == corrshape$double){
               pos<-pos+1
               loc2 <- this$Estimated$pars[pos]
             } else {
               loc2 <- NA
             }
             if (loc1 < min(this@st)) return(err_output)
             if (loc1 > max(this@st)) return(err_output)
             if (!is.na(loc2)) {
               if (loc2 < min(this@st)) return(err_output)
               if (loc2 > max(this@st)) return(err_output)
             }
             
             loc21<-loc1 # keep for comparison with the first loc1
             if (loc21<loc11) return(err_output)
             #if (loc21-loc11<0.2) return(err_output)
             
             
             #### ======== calculate loglikelihood ======== ####
             Pt <- calc.Pt2(this)  # T x nr.corPars
             
             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))
             
}

## -- loglik.stcc2.R() --####
loglik.stcc2.R = function(optimpars,z,stcc2Obj){
  
  err_output <- -1e10
  this <- stcc2Obj
  
  this$Estimated$pars <- tail(optimpars,2*this@nr.trPars)
  tmp.par <- optimpars
  
  # P1 always full of parameters:
  vP1 <- tmp.par[1:this@nr.corPars]
  mP <- unVecL(vP1)
  eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
  # Check for SPD - positive-definite check:
  if (min(eig$values) <= 0) return(err_output)
  this$Estimated$P1 <- mP
  
  #Remove the P1 corPars, then extract the P2 corPars
  tmp.par <- tail(tmp.par,-this@nr.corPars)
  nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
  freepars <- tmp.par[1:nr.freepars2]
  vP2<-vP1 # copy of P1
  # replace vP2 elements that have 1's in selvec by free pars
  vP2[which(as.logical(this$sel.vec2))]<-freepars
  mP <- unVecL(vP2)
  eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
  # Check for SPD - positive-definite check:
  if (min(eig$values) <= 0) return(err_output)
  this$Estimated$P2 <- mP
  
  #Remove the P2 corPars, then extract the P3 corPars
  tmp.par <- tail(tmp.par,-nr.freepars2)
  nr.freepars3 <- sum(this$sel.vec3,na.rm=TRUE) # count of number of free pars
  freepars <- tmp.par[1:nr.freepars3]
  vP3<-vP2 # copy of P2
  # replace vP3 elements that have 1's in selvec by free pars
  vP3[which(as.logical(this$sel.vec3))]<-freepars
  mP <- unVecL(vP3)
  eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
  # Check for SPD - positive-definite check:
  if (min(eig$values) <= 0) return(err_output)
  this$Estimated$P3 <- mP
  

  #### ======== constraint checks ======== ####
  
  # Check 2.1: Check the boundary values for speed params:
  pos <- 1
  speed <- this$Estimated$pars[pos]
  maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
  if (speed > maxSpeed) return(err_output)
  if (speed < 0) return(err_output)
  pos<- pos+1
  
  # Check 3.1: Check the locations fall within min-max values of st
  loc1 <- this$Estimated$pars[pos]
  if(this$shape == corrshape$double){
    pos<-pos+1
    loc2 <- this$Estimated$pars[pos]
  } else {
    loc2 <- NA
  }
  if (loc1 < min(this@st)) return(err_output)
  if (loc1 > max(this@st)) return(err_output)
  if (!is.na(loc2)) {
    if (loc2 < min(this@st)) return(err_output)
    if (loc2 > max(this@st)) return(err_output)
  }
  loc11<-loc1 # keep for later comparison with the second loc1
  pos<-pos+1
  
  # Check 2.2: Check the boundary values for speed params:
  speed <- this$Estimated$pars[pos]
  maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
  if (speed > maxSpeed) return(err_output)
  if (speed < 0) return(err_output)
  pos<- pos+1
  
  # Check 3.2: Check the locations fall within min-max values of st
  loc1 <- this$Estimated$pars[pos]
  if(this$shape == corrshape$double){
    pos<-pos+1
    loc2 <- this$Estimated$pars[pos]
  } else {
    loc2 <- NA
  }
  if (loc1 < min(this@st)) return(err_output)
  if (loc1 > max(this@st)) return(err_output)
  if (!is.na(loc2)) {
    if (loc2 < min(this@st)) return(err_output)
    if (loc2 > max(this@st)) return(err_output)
  }
  
  loc21<-loc1 # keep for comparison with the first loc1
  if (loc21<loc11) return(err_output)
  #if (loc21-loc11<0.2) return(err_output)
  
  
  #### ======== calculate loglikelihood ======== ####
  Pt <- calc.Pt2(this)  # T x nr.corPars
  
  llt <- vector("numeric")
  for(t in 1:this@Tobs) {
    mPt <- unVecL(Pt[t,,drop=FALSE])
    llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
  }
  # Return:
  return(sum(llt))
  
}

## --- estimateSTCC1 --- ####
myestimateSTCC1 = function(stcc1Obj,estimationCtrl){
             this <- stcc1Obj
             e <- this@e
             
             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose
             
             this$Estimated <- list()
             
             optimpars <- c( vecL(this$P1), vecL(this$P2), this$pars )
             optimpars <- optimpars[!is.na(optimpars)]
             
             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10
             
             # Filter the data:
             z <- w <- e <- this@e
             for(n in 1:this@N){
               w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv$g)
               z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch$h)
             }
             
             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.stcc1,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
             
             ### ---  Interpret the response from optim --- ###
             # An unhandled error could result in a NULL being returned by optim()
             if (is.null(tmp)) {
               this$Estimated$value <- -1e10
               this$Estimated$error <- TRUE
               return(this)
             }
             
             #Optim converged successfully => we expect tmp$par to have good estimates!
             if (tmp$convergence==0) {
               this$Estimated$error <- FALSE
               this$Estimated$value <- tmp$value
               
               tmp.par <- tmp$par
               this$Estimated$P1 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$P2 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$pars <- tmp.par
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)
               
               if (calcSE) {
                 cat("\nCalculating STCC standard errors...\n")
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
                 
                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   this$Estimated$P2.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   
                   this$Estimated$pars.se <- vecSE
                   if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- .calc.Pt(this)
               
             } else {
               #Failed to converge
               this$Estimated$error <- TRUE
               this$Estimated$value <- -1e10
               this$Estimated$optimoutput <- tmp
             }
             
             if (verbose) this$Estimated$optimoutput <- tmp
             #Return:
             return(this)
}


myestimateSTCC1.R = function(stcc1Obj,estimationCtrl){
  this <- stcc1Obj
  e <- this@e
  
  calcSE <- estimationCtrl$calcSE
  verbose <- estimationCtrl$verbose
  
  this$Estimated <- list()
  
  optimpars <- c( vecL(this$P1), vecL(this$P2), this$pars )
  optimpars <- optimpars[!is.na(optimpars)]
  
  ### ---  Call optim to calculate the estimate --- ###
  if (verbose) this$optimcontrol$trace <- 10
  
  # Filter the data:
  z <- w <- e <- this@e
  for(n in 1:this@N){
    w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv$g)
    z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch$h)
  }
  
  tmp <- NULL
  try(tmp <- optim(optimpars,.loglik.stcc1.R,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    this$Estimated$value <- -1e10
    this$Estimated$error <- TRUE
    return(this)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    
    this$Estimated$error <- FALSE
    this$Estimated$value <- tmp$value
    
    tmp.par <- tmp$par
    vP1<- tmp.par[1:this@nr.corPars]
    this$Estimated$P1 <- unVecL(vP1)
    tmp.par <- tail(tmp.par,-this@nr.corPars)
    
    nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
    freepars <- tmp.par[1:nr.freepars2]
    vP2<-vP1 # copy of P1
    # replace vP2 elements that have 1's in selvec by free pars
    vP2[which(as.logical(this$sel.vec2))]<-freepars
    this$Estimated$P2 <- unVecL(vP2)
    tmp.par <- tail(tmp.par,-nr.freepars2)
    
    this$Estimated$pars <- tmp.par
    if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
    names(this$Estimated$pars) <- names(this$pars)
    
    if (calcSE) {
      cat("\nCalculating STCC standard errors...\n")
      this$Estimated$hessian <- tmp$hessian
      vecSE <- vector("numeric")
      try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
      
      if(length(vecSE) > 0) {
        vecSE1 <- vecSE[1:this@nr.corPars]
        this$Estimated$P1.se <- unVecL(vecSE1)
        vecSE <- tail(vecSE,-this@nr.corPars)
        
        vecSE2<-vecSE1
        vecSE2[which(as.logical(this$sel.vec2))]<- vecSE[1:nr.freepars2]
        this$Estimated$P2.se <- unVecL(vecSE2)
        vecSE <- tail(vecSE,-nr.freepars2)
        
        this$Estimated$pars.se <- vecSE
        if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
        names(this$Estimated$pars.se) <- names(this$pars)
      }
    }
    
    this$Estimated$Pt <- .calc.Pt(this)
    
  } else {
    #Failed to converge
    this$Estimated$error <- TRUE
    this$Estimated$value <- -1e10
    this$Estimated$optimoutput <- tmp
  }
  
  if (verbose) this$Estimated$optimoutput <- tmp
  #Return:
  return(this)
}

## --- estimateSTCC2 --- ####
myestimateSTCC2 = function(stcc2Obj,estimationCtrl){
             this <- stcc2Obj
             e <- this@e
             
             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose
             
             this$Estimated <- list()
             
             optimpars <- c( vecL(this$P1), vecL(this$P2), vecL(this$P3), this$pars )
             optimpars <- optimpars[!is.na(optimpars)]
             
             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10
             
             # Filter the data:
             z <- w <- e <- this@e
             for(n in 1:this@N){
               w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv$g)
               z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch$h)
             }
             
             tmp <- NULL
             try(tmp <- optim(optimpars,loglik.stcc2,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
             
             ### ---  Interpret the response from optim --- ###
             # An unhandled error could result in a NULL being returned by optim()
             if (is.null(tmp)) {
               this$Estimated$value <- -1e10
               this$Estimated$error <- TRUE
               return(this)
             }
             
             #Optim converged successfully => we expect tmp$par to have good estimates!
             if (tmp$convergence==0) {
               this$Estimated$error <- FALSE
               this$Estimated$value <- tmp$value
               
               tmp.par <- tmp$par
               this$Estimated$P1 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$P2 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$P3 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               # TO DO : probaly going to fall over but likely never to be generlised to double
               pars1 <- tmp.par[1:this@nr.trPars]
               pars2 <- tail(tmp.par,-this@nr.trPars)
               this$Estimated$pars <- pars1
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               this$Estimated$pars <- c(this$Estimated$pars,pars2)
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)
               
               if (calcSE) {
                 cat("\nCalculating STCC standard errors...\n")
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
                 
                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   this$Estimated$P2.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   this$Estimated$P3.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   vecSE1 <- vecSE[1:this@nr.trPars]
                   vecSE2 <- tail(vecSE,-this@nr.trPars)
                   this$Estimated$pars.se <- vecSE1
                   if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   this$Estimated$pars.se <- c(this$Estimated$pars.se,vecSE2)
                   if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   #names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- calc.Pt2(this)
               
             } else {
               #Failed to converge
               this$Estimated$error <- TRUE
               this$Estimated$value <- -1e10
               this$Estimated$optimoutput <- tmp
             }
             
             if (verbose) this$Estimated$optimoutput <- tmp
             #Return:
             return(this)
}

myestimateSTCC2.R = function(stcc2Obj,estimationCtrl){
  this <- stcc2Obj
  e <- this@e
  
  calcSE <- estimationCtrl$calcSE
  verbose <- estimationCtrl$verbose
  
  this$Estimated <- list()
  
  optimpars <- c( vecL(this$P1), vecL(this$P2), vecL(this$P3), this$pars )
  optimpars <- optimpars[!is.na(optimpars)]
  
  ### ---  Call optim to calculate the estimate --- ###
  if (verbose) this$optimcontrol$trace <- 10
  
  # Filter the data:
  z <- w <- e <- this@e
  for(n in 1:this@N){
    w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv$g)
    z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch$h)
  }
  
  tmp <- NULL
  try(tmp <- optim(optimpars,loglik.stcc2.R,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    this$Estimated$value <- -1e10
    this$Estimated$error <- TRUE
    return(this)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    this$Estimated$error <- FALSE
    this$Estimated$value <- tmp$value
    
    tmp.par <- tmp$par
    vP1<- tmp.par[1:this@nr.corPars]
    this$Estimated$P1 <- unVecL(vP1)
    tmp.par <- tail(tmp.par,-this@nr.corPars)
    
    nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
    freepars <- tmp.par[1:nr.freepars2]
    vP2<-vP1 # copy of P1
    # replace vP2 elements that have 1's in selvec by free pars
    vP2[which(as.logical(this$sel.vec2))]<-freepars
    this$Estimated$P2 <- unVecL(vP2)
    tmp.par <- tail(tmp.par,-nr.freepars2)
    
    nr.freepars3 <- sum(this$sel.vec3,na.rm=TRUE) # count of number of free pars
    freepars <- tmp.par[1:nr.freepars3]
    vP3<-vP2 # copy of P2
    # replace vP3 elements that have 1's in selvec by free pars
    vP3[which(as.logical(this$sel.vec3))]<-freepars
    this$Estimated$P3 <- unVecL(vP3)
    tmp.par <- tail(tmp.par,-nr.freepars3)
    
    # TO DO : probably going to fall over but likely never to be generlised to double
    trPars <- tail(tmp$par,(2*this@nr.trPars))
    pars1 <- trPars[1:this@nr.trPars]
    pars2 <- tail(trPars,-this@nr.trPars)
    # First map the starting pars onto the Estimated pars, so we get the correct size & names
    #this$Estimated$pars <- this$pars
    # Then overwrite the non-NA pars with the Estimated values
    #this$Estimated$pars[which(!is.na(this$Estimated$pars))] <- pars1
    #
    if(this$shape != corrshape$double) this$Estimated$pars <- c(pars1,NA)
    if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,pars2,NA)
    names(this$Estimated$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")
    
    if (calcSE) {
      cat("\nCalculating STCC standard errors...\n")
      this$Estimated$hessian <- tmp$hessian
      vecSE <- vector("numeric")
      try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
      
      if(length(vecSE) > 0) {
        vecSE1 <- vecSE[1:this@nr.corPars]
        this$Estimated$P1.se <- unVecL(vecSE1)
        vecSE <- tail(vecSE,-this@nr.corPars)
        
        vecSE2<-vecSE1
        vecSE2[which(as.logical(this$sel.vec2))]<- vecSE[1:nr.freepars2]
        this$Estimated$P2.se <- unVecL(vecSE2)
        vecSE <- tail(vecSE,-nr.freepars2)
        
        vecSE3<-vecSE2
        vecSE3[which(as.logical(this$sel.vec3))]<- vecSE[1:nr.freepars3]
        this$Estimated$P3.se <- unVecL(vecSE3)
        vecSE <- tail(vecSE,-nr.freepars3)
        
        vecSE1 <- vecSE[1:this@nr.trPars]
        vecSE2 <- tail(vecSE,-this@nr.trPars)
        this$Estimated$pars.se <- vecSE1
        if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
        this$Estimated$pars.se <- c(this$Estimated$pars.se,vecSE2)
        if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
        names(this$Estimated$pars.se) <- names(this$pars)
      }
    }
    this$Estimated$Pt <- calc.Pt2(this)
    
  } else {
    #Failed to converge
    this$Estimated$error <- TRUE
    this$Estimated$value <- -1e10
    this$Estimated$optimoutput <- tmp
  }
  
  if (verbose) this$Estimated$optimoutput <- tmp
  #Return:
  return(this)
}

# Debug:
pars = startpars
DCClist = GCH

## DCC ####
myestimateDCC = function(z,pars,DCClist){
  # e = T x N returns
  # z: standardised returns:

  optimpars <- pars # (DCC a and b)
  Qbar <- cor(z) # N x N unoconditional correlation matrix
  Tobs <- nrow(z)
  N <- ncol(z)
  
  tmp <- optim(optimpars,myloglik.DCC,z,Qbar,gr=NULL,method="BFGS",control=list(fnscale=-1))
  
  #tmp <- NULL
  #try(tmp <- optim(optimpars,myloglik.DCC,z,Qbar,gr=NULL,method="BFGS"))

  ### ---  Interpret the response from optim --- ###

    #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    tmp$P <- matrix(0,nrow=Tobs,ncol=N*(N-1)/2)
    Qt_1 = Qbar
    a <- tmp$par[1]
    b <- tmp$par[2]
    for (t in 2:Tobs){
      Qt <- (1-a-b)*Qbar + a*t(z[t-1,,drop=FALSE])%*%z[t-1,,drop=FALSE] + b*Qt_1
      Pt <- diag((diag(Qt))^(-0.5))%*%Qt%*%diag((diag(Qt))^(-0.5))
      tmp$P[t,] <- vecL(Pt)
      Qt_1<-Qt
    }
    cat("Loglik Value: ", tmp$value)
    cat("Params: ", tmp$par)
  }
  return(tmp)
  
}

#Debug:

t=2

myloglik.DCC = function(optimpars,z,Qbar){
  err_output <- -1e10
  a <- optimpars[1]
  b <- optimpars[2]
  Tobs <- nrow(z)
  if (a<0) return(err_output)
  if (b<0) return(err_output)
  if ((a+b)>1) return(err_output)
  ll <- vector("numeric",Tobs)
  Qt_1 = Qbar
  for (t in 2:Tobs){
    Qt <- (1-a-b)*Qbar + a*t(z[t-1,,drop=FALSE])%*%z[t-1,,drop=FALSE] + b*Qt_1
    Pt <- diag((diag(Qt))^(-0.5))%*%Qt%*%diag((diag(Qt))^(-0.5))
    ll[t] <- -0.5*log(det(Pt)) -0.5*( z[t,,drop=FALSE] %*% (solve(Pt)) %*% t(z[t,,drop=FALSE]) )
    Qt_1<-Qt
  }
  
  return(sum(ll))
  
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
corrObj$P2 = unVecL(c(0.6,0,0.5))
corrObj$P3 = unVecL(c(0.2,0,0.2))

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

if(TRUE){
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

# check e has data with the expected correlation pattern
x<-100
tmp <- matrix(0,nrow=Tobs,ncol=N*(N-1)/2)
for(t in x:(Tobs-x)){
  tmp[t,]<-vecL(cor(e[(t-x+1):(t+x),]))
}
plot(tmp[,1],type='l')
plot(tmp[,2],type='l')
plot(tmp[,3],type='l')

# create stcc2 object to be estimated:
estCtrl = list(verbose=TRUE,calcSE=TRUE)
# how does the correlated e data get into NTV?? NTV was created using e=noise 
estCorr = stcc2(NTV)
estCorr@e <- e


# check extCorr has data with the expected correlation pattern
x<-100
tmp <- matrix(0,nrow=Tobs,ncol=N*(N-1)/2)
for(t in x:(Tobs-x)){
  tmp[t,]<-vecL(cor(estCorr@e[(t-x+1):(t+x),]))
}
plot(tmp[,1],type='l')
plot(tmp[,2],type='l')
plot(tmp[,3],type='l')

# Set Starting Pars:
estCorr$P1 = toeplitz(0.15^seq.int(0, N-1))
estCorr$P2 = toeplitz(0.55^seq.int(0, N-1))
estCorr$P3 = toeplitz(0.95^seq.int(0, N-1))
estCorr$P1 = unVecL(c(0.15,0,0.55))
estCorr$P2 = unVecL(c(0.65,NA,NA))
estCorr$P3 = unVecL(c(0.25,NA,0.25))
estCorr$pars = c(4,0.25,NA,4,0.75,NA)
names(estCorr$pars) <- c("speed1","loc11","loc12","speed2","loc21","loc22")
estCorr$sel.vec2 <- c(1,0,0)
estCorr$sel.vec3 <- c(1,0,1)

## PKG Estimator:  ####
estCorr = myestimateSTCC2(estCorr,estCtrl)
estCorr = myestimateSTCC2.R(estCorr,estCtrl)

plot(estCorr$Estimated$Pt[,1],type='l')
plot(estCorr$Estimated$Pt[,2],type='l',col="grey75")
plot(estCorr$Estimated$Pt[,3],type='l',col="green")

# estCorr_stcc2 <- estCorr
# estCorr_stcc2.R <- estCorr

saveRDS(estCorr_stcc2,"Results/estCorr_stcc2.RDS")
saveRDS(estCorr_stcc2.R,"Results/estCorr_stcc2.R.RDS")

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



## -- plot(stcc1) ####
setMethod("plot",signature = c(x="stcc1_class",y="missing"),
          function(x, y,...){
              #TODO: Allow override of main/title by user
              for(n in 1:x@N)
              plot(x$Estimated$Pt[,n],type='l')
              
          })

## -- plot(stcc1) ####
setMethod("plot",signature = c(x="stcc2_class",y="missing"),
          function(x, y,...){
              #TODO: Allow override of main/title by user
              for(n in 1:x@N)
                  plot(x$Estimated$Pt[,n],type='l')
              
          })



