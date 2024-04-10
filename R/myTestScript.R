
this <- stccObj
this$Estimated <- list()
z <- this@z

calcSE <- TRUE
verbose <- TRUE

veclP1 <- vecL(this$P1)
veclP2 <- vecL(this$P2)
veclP3 <- vecL(this$P3)


allPars <- c( veclP1, veclP2, veclP3, this$pars )
this@fixedPar.idx <- which(is.na(allPars))

if(any(is.na(allPars))){
    this@fixedPar.idx <- which(is.na(allPars))
    optimpars <- allPars[-this@fixedPar.idx]
    # Manage the optimControl params to deal with any NA's above
    if(length(this$optimcontrol$parscale) == length(optimpars)){
        # Have to assume the user made correct changes
    } else {
        this$optimcontrol$parscale <- this$optimcontrol$parscale[-this@fixedPar.idx]
    }
    #
    if(length(this$optimcontrol$ndeps) == length(optimpars)){
        # Have to assume the user made correct changes
    } else {
        this$optimcontrol$ndeps <- this$optimcontrol$ndeps[-this@fixedPar.idx]
    }
    
} else optimpars <- allPars
#Note: NA's can be sent due to restricted (constant cor) parameters


### ---  Call optim to calculate the estimate --- ###
if (verbose) this$optimcontrol$trace <- 10


# Run Loglik:

err_output <- -1e10

# Rebuild 'this' with the values from optimpars
tmp.par <- optimpars
if(length(this@fixedPar.idx) > 0 ){
    # First reset the optimpars to the original length with NA's
    tmp.par <- vector.insert(optimpars,this@fixedPar.idx,rep(NaN,length(this@fixedPar.idx)))
    # Add the previous correlation param in place of the NA.  This will be put into the $Estimated$Pn matrices
    tmp.par <- .fixCorPars(tmp.Par,this@fixedPar.idx,this@nr.corPars)
}

this$Estimated$pars <- tail(tmp.par,2*this@nr.trPars)

vecP <- tmp.par[1:this@nr.corPars]
mP <- unVecL(vecP)
eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
# Check for SPD - positive-definite check:
if (min(eig$values) <= 0) return(err_output)
this$Estimated$P1 <- mP
# Now drop these values from tmp.par
tmp.par <- tail(tmp.par,-this@nr.corPars)
#
vecP <- tmp.par[1:this@nr.corPars]
mP <- unVecL(vecP)
eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
# Check for SPD - positive-definite check:
if (min(eig$values) <= 0) return(err_output)
this$Estimated$P2 <- mP
# Now drop these values from tmp.par
tmp.par <- tail(tmp.par,-this@nr.corPars)
#
vecP <- tmp.par[1:this@nr.corPars]
mP <- unVecL(vecP)
eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
# Check for SPD - positive-definite check:
if (min(eig$values) <= 0) return(err_output)
this$Estimated$P3 <- mP




#### ======== constraint checks ======== ####

# Check 1: Check the boundary values for speed1 param:
speed <- this$Estimated$pars[1]
maxSpeed <- switch(this$speedopt,1000,(1000/sd(this$st)),7.0,0.30)
if (speed > maxSpeed) return(err_output)
if (speed < 0) return(err_output)

# Check 2: Check the first location falls within min-max values of st
loc1 <- this$Estimated$pars[2]
if (loc1 < min(this$st)) return(err_output)
if (loc1 > max(this$st)) return(err_output)

# Check 3: Check the boundary values for speed2 param:
speed <- this$Estimated$pars[3]
maxSpeed <- switch(this$speedopt,1000,(1000/sd(this$st)),7.0,0.30)
if (speed > maxSpeed) return(err_output)
if (speed < 0) return(err_output)

# Check 4: Check the second location falls within min-max values of st
loc2 <- this$Estimated$pars[4]
if (loc2 < min(this$st)) return(err_output)
if (loc2 > max(this$st)) return(err_output)

# Check 5: Check that loc1 is before loc 2
if (loc2 < loc1) return(err_output)


#### ======== calculate loglikelihood ======== ####
#Pt <- .calc.Pt2(this)  # T x nr.corPars
vP1 <- matrix(vecL(this$Estimated$P1),nrow=1)
vP2 <- matrix(vecL(this$Estimated$P2),nrow=1)
vP3 <- matrix(vecL(this$Estimated$P3),nrow=1)

#Gt <- .calc.Gt2(this) # T x 2

speed1 <- this$Estimated$pars[1]
loc1 <- this$Estimated$pars[2]
speed2 <- this$Estimated$pars[3]
loc2 <- this$Estimated$pars[4]

st_c_1 <- this$st - loc1
st_c_2 <- this$st - loc2

G <- matrix(0,nrow = this$Tobs, ncol = 2)
if(this$speedopt == corrspeedopt$gamma) { G[,1] <- 1/(1+exp(-speed1*st_c_1)) }
if(this$speedopt == corrspeedopt$gamma_std) { G[,1] <- 1/(1+exp(-speed1/sd(this$st)*st_c_1)) }
if(this$speedopt == corrspeedopt$eta) { G[,1] <- 1/(1+exp(-exp(speed1)*st_c_1)) }
#
if(this$speedopt == corrspeedopt$gamma) { G[,2] <- 1/(1+exp(-speed2*st_c_2)) }
if(this$speedopt == corrspeedopt$gamma_std) { G[,2] <- 1/(1+exp(-speed2/sd(this$st)*st_c_2)) }
if(this$speedopt == corrspeedopt$eta) { G[,2] <- 1/(1+exp(-exp(speed2)*st_c_2)) }

Gt <- (matrix(G,nrow = this$Tobs,ncol = 2))


Pt <- ((1-Gt[,2])*(1-Gt[,1]))%*%vP1 + ((1-Gt[,2])*Gt[,1])%*%vP2 + Gt[,2]%*%vP3 # T x N(N-1)/2
if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1)
 

optimpars <- optimpars * 0.86

## -- .fixCorPars -- ####
.fixCorPars = function(optimpars,fixed.idx,nr.corPars){
    # Function Logic:
    for(n in 1:length(fixed.idx)){ optimpars[fixed.idx[n]] <- optimpars[fixed.idx[n]-nr.corPars] }
    return(optimpars)
}

x = c(1,2,3,4,NA,6,7,NA,9,10,11)
napos = which(is.na(x))
tmp.par <- .fixCorPars(x,napos,3)

pars <- tail(tmp.par,2)

vecP <- tmp.par[1:3]
P1 <- unVecL(vecP)
# Now drop these values from tmp.par
tmp.par <- tail(tmp.par,-3)
#
vecP <- tmp.par[1:3]
P2 <- unVecL(vecP)
# Now drop these values from tmp.par
tmp.par <- tail(tmp.par,-3)
#
vecP <- tmp.par[1:3]
P3 <- unVecL(vecP)




