estimateGARCH_old <- function(e,garchObj,estimationControl,tvObj){
    this <- garchObj
    
    if(this$type == garchtype$noGarch) {
        message("Cannot estimateGARCH for type: NoGarch")
        return(this)
    }
    
    # Attach results of estimation to the object
    this$Estimated <- list()
    this$Estimated$method <- "MLE"
    
    if (isTRUE(estimationControl$verbose)) {
        this$optimcontrol$trace <- 10
        cat("\nEstimating GARCH object...\n")
    } else this$optimcontrol$trace <- 0
    
    # Get Optimpars from garch$pars
    optimpars <- as.vector(this$pars)
    names(optimpars) <- rownames(this$pars)
    
    # Now call optim:
    tmp <- NULL
    try(tmp <- optim(optimpars,loglik.garch.univar,gr=NULL,e,this,tvObj, method="BFGS",control=this$optimcontrol,hessian = estimationControl$calcSE))
    
    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
        this$Estimated$value <- -Inf
        this$Estimated$error <- TRUE
        warning("estimateGARCH() - optim failed unexpectedly and returned NULL. Check the optim controls & starting params")
        return(this)
    }
    if (tmp$convergence!=0) {
        this$Estimated$value <- -Inf
        this$Estimated$error <- TRUE
        this$Estimated$optimoutput <- tmp
        warning("estimateGARCH() - failed to converge. Check the optim controls & starting params")
        return(this)
    }
    
    this$Estimated$value <- tmp$value
    this$Estimated$error <- FALSE
    
    this$Estimated$pars <- .parsVecToMatrix(this,tmp$par)
    # Get conditional variance
    
    this@h <- .calculate_h(this,e/sqrt(tvObj@g))
    
    # Calc Std Errors
    if (isTRUE(estimationControl$calcSE)) {
        cat("\nCalculating GARCH standard errors...\n")
        #this$Estimated$hessian <- optimHess(tmp$par,loglik.garch.univar,gr=NULL,e,this,tvObj,control=this$optimcontrol)
        this$Estimated$hessian <- tmp$hessian
        StdErrors <- NULL
        try(StdErrors <- sqrt(-diag(invertMatrix_NPD(tmp$hessian))))
        if(is.null(StdErrors)) {
            this$Estimated$se <- matrix(NA,nrow=this@nr.pars)
        }else {
            this$Estimated$se <- matrix(StdErrors,nrow=this@nr.pars)
        }
        rownames(this$Estimated$se) <- rownames(this$pars)
        colnames(this$Estimated$se) <- "se"
    }
    if (isTRUE(estimationControl$verbose)) this$Estimated$optimoutput <- tmp
    
    return(this)
}

estimateGARCH_new <- function(e,garchObj,estimationControl,tvObj){
    this <- garchObj
    
    if(this$type == garchtype$noGarch) {
        message("Cannot estimateGARCH for type: NoGarch")
        return(this)
    }
    
    # Attach results of estimation to the object
    this$Estimated <- list()
    this$Estimated$method <- "MLE"
    
    if (isTRUE(estimationControl$verbose)) {
        this$optimcontrol$trace <- 10
        cat("\nEstimating GARCH object...\n")
    } else this$optimcontrol$trace <- 0
    
    # Get Optimpars from garch$pars
    optimpars <- as.vector(this$pars)
    names(optimpars) <- rownames(this$pars)
    
    # Now call optim:
    tmp <- NULL
    try(tmp <- optim(optimpars,loglik.garch.univar,gr=NULL,e,this,tvObj, method="BFGS",control=this$optimcontrol))
    
    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
        this$Estimated$value <- -Inf
        this$Estimated$error <- TRUE
        warning("estimateGARCH() - optim failed unexpectedly and returned NULL. Check the optim controls & starting params")
        return(this)
    }
    if (tmp$convergence!=0) {
        this$Estimated$value <- -Inf
        this$Estimated$error <- TRUE
        this$Estimated$optimoutput <- tmp
        warning("estimateGARCH() - failed to converge. Check the optim controls & starting params")
        return(this)
    }
    
    this$Estimated$value <- tmp$value
    this$Estimated$error <- FALSE
    
    this$Estimated$pars <- .parsVecToMatrix(this,tmp$par)
    # Get conditional variance
    
    this@h <- .calculate_h(this,e/sqrt(tvObj@g))
    
    # Calc Std Errors
    if (isTRUE(estimationControl$calcSE)) {
        cat("\nCalculating GARCH standard errors...\n")
        this$Estimated$hessian <- optimHess(tmp$par,loglik.garch.univar,gr=NULL,e,this,tvObj,method="BFGS",control=this$optimcontrol)
        #this$Estimated$hessian <- tmp$hessian
        StdErrors <- NULL
        try(StdErrors <- sqrt(-diag(invertMatrix_NPD(tmp$hessian))))
        if(is.null(StdErrors)) {
            this$Estimated$se <- matrix(NA,nrow=this@nr.pars)
        }else {
            this$Estimated$se <- matrix(StdErrors,nrow=this@nr.pars)
        }
        rownames(this$Estimated$se) <- rownames(this$pars)
        colnames(this$Estimated$se) <- "se"
    }
    if (isTRUE(estimationControl$verbose)) this$Estimated$optimoutput <- tmp
    
    return(this)
}

.parsVecToMatrix <-   function(garchObj,pars){
   this <- garchObj
   
   if(this$type == garchtype$noGarch) {
       message("Cannot create Garch Params for type: NoGarch")
       return(this)
   }
   
   maxLag <- max(this@order)
   
   # Set the row names:
   garchparsRownames <- c("omega","alpha","beta","gamma")
   # Return the formatted matrix
   matrix(pars,nrow = this@nr.pars ,ncol = maxLag,dimnames = list(garchparsRownames[1:this@nr.pars],"Est"))
   
}

.calculate_h <- function(garchObj,e){
    this <- garchObj
    
    if(this$type == garchtype$noGarch) return(this@h)
    
    Tobs <- NROW(e)
    h <- rep(0,Tobs)
    h[1] <- sum(e*e)/Tobs
    
    # TODO: Extend the below to handle more lags (higher order Garch)
    for(t in 2:Tobs) {
        h[t] <- this$Estimated$pars["omega",1] + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1]
        if(this$type == garchtype$gjr) h[t] <- h[t] + this$Estimated$pars["gamma",1]*(min(e[t-1],0))^2
    }
    
    return(h)
    
}


Tobs = NROW(e)
tvObj <- tv(1,tvshape$delta0only)
tvObj@g <- rep(1,Tobs)

estCtrl = list(verbose=TRUE,calcSE=TRUE)

Garch1 <- garch(garchtype$general)
Garch1 <- estimateGARCH_old(e,Garch1,estCtrl,tvObj)
summary(Garch1)

Garch2 <- garch(garchtype$general)
Garch2 <- estimateGARCH_old(e,Garch2,estCtrl,tvObj)
summary(Garch2)

identical(Garch1$Estimated$hessian,Garch2$Estimated$hessian)

