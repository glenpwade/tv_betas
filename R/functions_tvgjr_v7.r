#### =========================================== File Notes =========================================== ####
###                                                                                                     ###
###   TVGJR:  All the code in this functions file uses custom-defined list structures for TV & GARCH    ###
###           To use this file effectively, please read the TV & GARCH object definitions carefully.    ###
###                                                                                                     ###
###           The main functions in this file are:                                                      ###
###           * EstimateTV()                                                                            ###
###           * EstimateGARCH()                                                                         ###
###           * CalcProbabilityDist()                                                                   ###
###           Think of these as 'public methods'.  They provide a reasonable level of error-handling.   ###
###                                                                                                     ###
###                                                                                                     ###
###           Latest Ver.7.0.0                                                                          ###
###                                                                                                     ###
###     Change Log:  Ver.2.x.x
###         
###     27-Feb-2017 - List Objects changed, to remove need for position vectors                        
###     03-Mar-2017 - $var_target added to TV & GARCH objects
###                 - EstimateGARCH() modified to use compiled LogLik Fn for better performance   
###     10-Mar-2017 - Cut-down version of this file saved as: 'functions_tvgjr_WIG20.r'                 ###
###     18-Mar-2017 - EstimateTVGARCH() function finished & tested                                      ###
###
###     Change Log:  Ver.4.x.x
###         
###     20-Feb-2018 - CalculateTest_TV() changed to include a 'useRcpp' parameter
###                 - CalcProbabilityDist() changed to fix errors with length of runSimRow
###                 - dg.dt, myG & calculate_g re-written to improve performance
###                 - TV$linparsMode removed as it is not needed.
###                 - TV$Tobs (Total number of observations in sample) added.
###                 - Tests added for CCC_STCC correlation...  Work in progress
###
###                                                                                                     ###
### =================================================================================================== ###


####  ========================== Main_Object_Definitions ============================ ####

  ### ============================================== ###
  ###        TV list object definition               ###
  ### ============================================== ###
  
  # TV <- list()
  # TV$delta0 <- vector("numeric",1L)   # (MANDATORY) delta0 - single element vector (not integer) since we need to use it in vector functions 
  # TV$pars <- vector() or NULL         # delta1, speed1, loc1.1, <loc1.2>=NA, delta2, speed2, loc2.1, <loc2.2>,... loc2.n MUST provide NA if no value
  # TV$linpars <- vector() or NULL      # [l1,l2,l3,l4] = constant, linear, quadratic, cubic terms. MUST provide NA if no value.  All code expects a 4-element vector!
  # TV$var_target <- logical()          # Variance Targeting switch. BOOL = FALSE => delta0 will NOT be used in estimation!!!  Always set = TRUE for this code.
  # TV$optimcontrol <- list()           # Control list passed to the optim() function.  Must be set for ALL parameters - Estimate.. functions will handle var targetting.
  # TV$speedoption <- integer()         # speedoptions are: 1=gamma, 2=gamma/sd(st), 3=exp(eta), 4=1/lambda^2
  # TV$shape <- vector() or NULL        # vector of integers describing the transitions (G's)
  #                                     #    E.g.  c(1,2,3) => 3 transitions  Must be NULL for TV Order0, as 'is.na(TV$shape[1])' is used to test for Order=0.
  #                                     #    1=Single location, 2=Double location, 3=Single location, with Squared transition
  # TV$st <- vector()                   # smooth transition variable.  Typically = seq(1,Tobs)/Tobs, i.e. linear transition
  ##
  ## -- Elements returned from the EstimateTV() function: -- ##
  ##
  # TV$Estimated$delta0 <- vector("numeric",1L)   # (MANDATORY) delta0 - single element vector (not integer) since we need to use it in vector functions
  # TV$Estimated$lastdelta0 <- vector("numeric",1L)   # keeps a copy of the last good estimate of delta0, to be used when variance targetting
  # TV$Estimated$pars <- vector() or NULL         # delta1, speed1, loc1.1, <loc1.2>=NA, delta2, speed2, loc2.1, <loc2.2>,... loc2.n MUST provide NA if no value
  # TV$Estimated$linpars <- vector() or NULL      # [l1,l2,l3,l4] = constant, linear, quadratic, cubic terms. MUST provide NA if no value.  All code expects a 4-element vector!
  # TV$Estimated$value <- numeric()               # numeric value returned from the optim() function for the log-liklihood value
  # TV$Estimated$hessian <- matrix()                 # Hessian matrix of TV$pars, returned from optim() or optimHess()
  # TV$Estimated$error <- logical()               # True/False.  Set to TRUE when optim fails, False otherwise.
  # TV$Estimated$stderr <- matrix()               # Residual Std Err matrix calculated from the returned Hessian matrix - Often causes a numerical error
  # TV$Estimated$condvars <- vector()             # vector of conditional variances, 'g'
  # TV$Estimated$optimoutput <- list()            # list returning the full output from the optim()


  
  
  ### ============================================== ###
  ###       GARCH list object definition             ###
  ### ============================================== ###
  
  # GARCH <- list()
  # GARCH$type <- numeric()          # 0: Constant GARCH, 1:GARCH, 2:GJR-GARCH, 3: GJR-GARCH with omega=0
  # GARCH$pars <- vector()           # omega, alpha, beta, delta  (Single-order Garch only implemented)
  #										               #    Note: All the code relies on all 4 parameters being provided, so don't forget to include NA where needed (Garch Types 0,1 & 3)!!
  # GARCH$var_target <- logical()    # Variance Targeting switch. BOOL = FALSE => All $pars will be used in estimation. (TRUE => omega will be calculated, not free)
  #										               #    Note: Variance targetting is not used in this code, but comes into play for TV with GARCH estimation.
  # GARCH$optimcontrol <- list()     # Control list passed to the optim() function.  Must be set for ALL parameters - Estimate.. functions will handle var targetting.
  ##
  ## -- Elements returned from the EstimateGARCH() function: -- #
  ##
  # GARCH$Estimated$value <- numeric()         # numeric value returned from the optim() function for the log-liklihood value
  # GARCH$Estimated$pars <- vector()           # omega, alpha, beta, delta  (Single-order Garch only implemented)
  # GARCH$Estimated$condvars <- vector()       # vector of conditional variances 'h'
  # GARCH$Estimated$error <- logical()         # True/False.  Set to TRUE when optim fails, False otherwise
  # GARCH$Estimated$hessian <- matrix()           # Hessian matrix of GARCH$pars, returned from optim() or optimHess()
  # GARCH$Estimated$stderr <- matrix()         # Residual Std Err matrix calculated from the returned Hessian matrix - Often causes a numerical error
  
  
  ### ============================================== ###
  ###         CORR list object definition            ###
  ### ============================================== ###
  
  # CORR <- list()                 # This object is just a wrapper container for all the specific correlation models
  # CORR$type <- vector(char)      # "CCC", "STCC", "STEC", ...
  # CORR$CCC <- list()
  # CORR$STCC <- list()
  # CORR$STEC <- list()
  #
  
  ### ============================================== ###
  ###         STCC list object definition            ###
  ### ============================================== ###
  
  # STCC <- list()
  # STCC$type <- as.numeric()             # 0:"CCC", 1:"STCC"
  # STCC$P1pars <- vector()
  # STCC$P2pars <- vector()
  # STCC$TRpars <- as.vector()            # A vector of the remaining pars: Speed, Loc1, <Loc2>  
  # STCC$speedoption <- as.integer()      # speedoptions are: 1=gamma 2=gamma/sd(st) 3=exp(eta) 4=...
  # STCC$st <- as.vector()                # transition variable.  Default = seq(1,Tobs)/Tobs, i.e. linear
  # STCC$shape <- as.vector()             # vector of integers describing the transitions (G's)
  #                                       # e.g. c(NA) => No transitions.  c(1,2,3) => 3 transitions
  #                                       # 1=Single location, 2=Double location, 3=Single location, with Squared transition
  # -- Elements returned from the Optim() function: -- #
  # STCC$Estimated$P1pars <- vector()
  # STCC$Estimated$P2pars <- vector()
  # STCC$Estimated$TRpars <- as.vector()      
  # STCC$Estimated$condcorrs <- as.matrix()         # time-series of conditional correlations, 'rho'  (Size = T x N(N-1)/2)
  # STCC$Estimated$stderr <- as.matrix()            # Residual Std Err matrix calculated from the returned Hessian matrix  
  # STCC$Estimated$value <- as.numeric()            # numeric value returned from the optim() function for the log-liklihood value
  # STCC$Estimated$hessian <- as.matrix()              # Hessian matrix returned from optim() (or optimHess())
  
  
  ### ============================================== ###
  ###         STEC list object definition            ###
  ### ============================================== ###
  
  # STEC <- list()
  # STEC$type <- as.numeric()             # 0:"CEC" (only 1 par = p1),  1:"STEC" (up to 5 pars)
  # STEC$pars <- as.vector()              # A vector containing the parameters, p1, p2, Speed, Loc1, <Loc2> (all scalars)
  # STEC$linpars <- as.vector() or NULL   # linpar1, linpar2,...  Note: NULL (or non-existing) => All functions will operate with linearised = FALSE
  # STEC$speedoption <- as.integer()      # speedoptions are: 1=gamma 2=gamma/sd(st) 3=exp(eta) 4=...
  # STEC$st <- as.vector()                # transition variable.  Default = seq(1,Tobs)/Tobs, i.e. linear
  # STEC$shape <- as.vector()             # vector of integers describing the transitions (G's)
  #                                       # e.g. c(NA) => No transitions.  c(1,2,3) => 3 transitions
  #                                       # 1=Single location, 2=Double location, 3=Single location, with Squared transition
  # STEC$condcorrs <- as.matrix()         # time-series of conditional correlations, 'rho'  (Size = T x N(N-1)/2)
  # STEC$stderr <- as.matrix()            # Residual Std Err matrix calculated from the returned Hessian matrix
  # -- Elements returned from the Optim() function: -- #
  # STEC$value <- as.numeric()            # numeric value returned from the optim() function for the log-liklihood value
  # STEC$hess <- as.matrix()              # Hessian matrix returned from optim() (or optimHess())
  
  
  
  ### ++++++++++++ MULTIVARIATE OBJECTS - START ++++++++++++ ###
  
  #nTV <- list(N)              # List of TV() list objects, where N = total number of objects
  #nTV[[n]]                    # Reference to one specific TV() list object
  #nTV$delta0 <- as.vector()   # Vector containing c(TV[[1]]$delta0, TV[[2]]$delta0,..., TV[[N]]$delta0)
  #nTV$pars <- as.vector()     # Vector containing c(TV[[1]]$pars, TV[[2]]$pars,..., TV[[N]]$pars)
  #nTV$hessian <- as.matrix       # Hessian matrix of all TV$pars, returned from optim() (or optimHess())
  #nTV$stderr <- as.matrix     # Residual Std Err matrix calculated from the returned Hessian matrix
  
  #nGARCH <- list(N)              # List of GARCH() list objects, where N = total number of objects
  #nGARCH[[n]]                    # Reference to one specific GARCH() list object
  #nGARCH$pars <- as.vector()     # Vector containing c(GARCH[[1]]$pars, GARCH[[2]]$pars,..., GARCH[[N]]$pars)
  #nGARCH$hess <- as.matrix       # Hessian matrix of all GARCH$pars (with Omega removed), returned from optim() (or optimHess())
  #nGARCH$stderr <- as.matrix     # Residual Std Err matrix calculated from the returned Hessian matrix
  
  
  # Note: STCC is by definition a multivariate object.  We cannot have a univariate STCC
  
  ### ++++++++++++  MULTIVARIATE OBJECTS - END  ++++++++++++ ###
  
####  ========================== Main_Object_Definitions ============================ ###


#### ======================================== Initialisation =========================================== ####


## --- Define Global Variables --- ##
err_output <- -1e10

## ---  Define Global Enum's  --- ##
# Transition shape
TRshape <- list()
TRshape$none = as.integer(0)
TRshape$single = as.integer(1)
TRshape$double = as.integer(2)
TRshape$triple = as.integer(3)

# Transition speed option
TRspeedopt <- list()
TRspeedopt$gamma <- as.integer(1)
TRspeedopt$gamma_std <- as.integer(2)
TRspeedopt$eta <- as.integer(3)
TRspeedopt$lamda2_inv <- as.integer(4)

GARCHtype <- list()
GARCHtype$none <- NA
GARCHtype$constant <- as.integer(0)
GARCHtype$general <- as.integer(1)
GARCHtype$GJR <- as.integer(2)
GARCHtype$GJR_alpha0 <- as.integer(3)

STCCtype <- list()
STCCtype$CCC <- as.integer(0)
STCCtype$STCC <- as.integer(1)

TESTtype <- list()
TESTtype$TR2 <- as.integer(1)
TESTtype$Robust <- as.integer(2)
TESTtype$ALL <- as.integer(3)

TESTorder <- list()
TESTorder$H0 <- as.integer(0)
TESTorder$H01 <- as.integer(1)
TESTorder$H02 <- as.integer(2)
TESTorder$H03 <- as.integer(3)

### ======================================== Initialisation =========================================== ###


#### ================ TV class ================ ####
# Define the TV object class:
newTV <- setClass("tv_class",contains = c("namedList"))

setMethod("initialize","tv_class",
          function(.Object,...){
            # Initialise fields with defaults:
            .Object$st <- NA
            .Object$Tobs <- NA
            .Object$delta0 <- 1.0
            .Object$pars <- NA
            .Object$linpars <- NA
            .Object$shape <- TRshape$none
            .Object$speedoption <- TRspeedopt$eta
            .Object$var_target <- TRUE
            .Object$optimcontrol <- list(fnscale = -1)
            # Return:
            .Object
          })

setGeneric(name="newTV",def=function(st,del0,vecpars,shape,speedopt,...){
  standardGeneric("newTV")
}, valueClass = "tv_class", signature = c("st","del0","vecpars","shape","speedopt"))

# Create default constructor:
setMethod(f="newTV",signature = c("numeric","numeric","numeric","numeric","numeric"),
          definition = function(st,del0,vecpars,shape,speedopt, ...){
            # Validate params
            if (is.null(st)) stop("Error: missing 'st'\nAll parameters are required to construct a TV object")
            if (is.null(del0)) stop("Error: missing 'del0'\nAll parameters are required to construct a TV object")
            if (is.null(vecpars)) stop("Error: missing 'vecpars'\nAll parameters are required to construct a TV object")
            if (is.null(shape)) stop("Error: missing 'shape'\nAll parameters are required to construct a TV object")
            if (is.null(speedopt)) stop("Error: missing 'speedopt'\nAll parameters are required to construct a TV object")
            #
            tv <- new("tv_class") 
            tv$st <- st
            tv$Tobs <- length(st)  
            tv$delta0 <- del0
            tv$shape <- shape
            tv$pars <- constructTVpars(vecpars,shape)
            tv$speedoption <- speedopt
            tv$linpars <- NA
            # Return:
            tv
          })

# Create special-case 'no transitions' constructor - just requires delta0 & st
setMethod(f="newTV",signature = c("numeric","numeric","missing","missing","missing"),
          definition = function(st,del0, ...){
            tv <- new("tv_class")
            if(!is.null(st)) {
              tv$Tobs <- length(st) 
              tv$st <- st
            } else stop("Error: missing 'st' - Cannot construct TV object")
            if(!is.null(del0)) tv$delta0 <- del0 else tv$delta0 <- 1.0
            # Return:
            tv
          })

## TODO: It would be nice to 'pass by ref' for these methods
setGeneric(name="addH03pars",def=function(tv){
  standardGeneric("addH03pars")
}, valueClass = "tv_class", signature = c("tv"))

setMethod(f="addH03pars",signature = c("tv_class"),
          definition = function(tv){
            if(!is.null(tv)) {
              tv$linpars <- c(NA,1,1)
            } else stop("Error: missing 'tv' - Cannot add H03 params to NULL")
            # Return:
            tv
          })

setGeneric(name="addH02pars",def=function(tv){
  
  standardGeneric("addH02pars")
}, valueClass = "tv_class", signature = c("tv"))

setMethod(f="addH02pars",signature = c("tv_class"),
          definition = function(tv){
            if(!is.null(tv)) {
              tv$linpars <- c(NA,1)
            } else stop("Error: missing 'tv' - Cannot add H02 params to NULL")
            # Return:
            tv
          })

# Create summary() method
setMethod(f="summary",signature = c("tv_class"),
          definition = function(object, ...){
            if(is.null(object)) stop("Error in summary(tv) :\n argument 'tv' is missing,  with no default")
            #
            tv <- object
            cat("\nCurrent estimation results for this TV object:\n")
            cat("\nError:",tv$Estimated$error)
            cat("\nLoglikelihood Value:",tv$Estimated$value)
            cat("\n\nModel Parameters:")
            cat("\ndelta0:",tv$Estimated$delta0,"\nPars:\n")
            print(tv$Estimated$pars)
            cat("\nlastdelta0",tv$Estimated$lastdelta0,"\n")
            cat("\nTo see the full object, including starting params and optim controls, use print(tv)")
            cat("\nTo view the conditional variances, or smooth transition variable, use plot(tv)")
            # Return:
            
          })

# Create print() method
setMethod(f="print",signature = c("tv_class"),
          definition = function(x=tv_class, ...){
            if(is.null(x)) stop("Error in print(tv) :\n argument 'tv' is missing,  with no default")
            #
            tv <- x
            cat("\nSample size:",tv$Tobs)
            cat("\nStarting values for this TV object:\n")
            cat("\ndelta0:",tv$delta0,"\nPars:\n")
            print(tv$pars)
            cat("\nshape:",as.character(tv$shape))
            cat("\nspeedoption:",tv$speedoption)
            cat("\n\noptim controls:\n")
            print(tv$optimcontrol)
            cat("\nVariance Targeting:",tv$var_target)
            cat("\n\nCurrent estimation results for this TV object:\n")
            cat("\nError:",tv$Estimated$error)
            cat("\nLoglikelihood Value:",tv$Estimated$value)
            cat("\n\nModel Parameters:")
            cat("\ndelta0:",tv$Estimated$delta0,"\nPars:\n")
            print(tv$Estimated$pars)
            cat("\nlastdelta0",tv$Estimated$lastdelta0,"\n")
            cat("\nTo see a summary, use summary(tv)")
            cat("\nTo view the conditional variances, or smooth transition variable, use plot(tv)")
            # Return:
            
          })

setGeneric(name="plot.tv",def=function(tv,type, ...){
  standardGeneric("plot.tv")
}, signature = c("tv","type"))

# Create plot() method
setMethod(f="plot.tv", signature = c("tv_class"),
          definition = function(tv="tv_class",...){
            if(is.null(tv)) stop("Error in plot(tv) :\n argument 'tv' is missing,  with no default")
              plot(tv$Estimated$condvars,type='l', ...)
            })


#### ================ GARCH class ================ ####

# Define the GARCH object class:
newGARCH <- setClass("garch_class",contains = c("namedList"))

setMethod("initialize","garch_class",
          function(.Object,...){
            # Initialise fields with defaults:
            .Object$type <- GARCHtype$general
            .Object$pars <- c(0.05,0.05,0.85)
            .Object$var_target <- FALSE
            .Object$optimcontrol <- list(fnscale = -1)
            # Return:
            .Object
          })

setGeneric(name="newGARCH",def=function(type,pars,...){
  standardGeneric("newGARCH")
}, valueClass = "garch_class", signature = c("type","pars"))

# Create default constructor:
sig1 <- c("numeric","numeric")
setMethod(f="newGARCH",signature = sig1,
          definition = function(type,pars, ...){
            # Validate params
            OutofRange <- FALSE 
            if (!is.integer(type)) OutofRange <- TRUE  else {
              if (type < 0 & type > 3) OutofRange <- TRUE   
            }
            if (OutofRange) stop("Error: Invalid 'type'\nAll parameters are required to construct a GARCH object.  Use GarchType$...")
            
            if (!is.vector(pars,mode = "numeric") || is.na(pars)) OutofRange <- TRUE else {
              if (length(pars) > 4) OutofRange <- TRUE 
            }
            if(OutofRange) stop("Error: Invalid 'pars'\nAll parameters are required to construct a GARCH object.  Please use a numeric vector, max length = 4")
            #
            garch <- new("garch_class") 
            garch$type <- type
            garch$pars <- constructGARCHpars(garch,pars)
            if (type == GARCHtype$general) {
              nDeps=c(1e-4,1e-4,1e-9)
              garch$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3),ndeps=nDeps)  
            } else if (type == GARCHtype$GJR) {
              nDeps=c(1e-4,1e-4,1e-9,1e-4)
              garch$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=nDeps)  
            } else if (type == GARCHtype$GJR_alpha0) {
              nDeps=c(1e-4,1e-4,1e-9,1e-4)
              garch$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=nDeps)  
            }

            # Return:
            garch
          })

# Create special-case 'General Garch' constructor - no params required.
sig2 <- c("missing","missing")
setMethod(f="newGARCH",signature = sig2,
          definition = function(...){
            garch <- new("garch_class") 
            # Return:
            garch
          })


#### ================ STCC class ================ ####

# Define the GARCH object class:
newSTCC <- setClass("stcc_class",contains = c("namedList"))

setMethod("initialize","stcc_class",
          function(.Object,...){
            # Initialise fields with defaults:
            .Object$type <- STCCtype$STCC
            .Object$Tobs <- 0
            .Object$st <- NA
            .Object$P1 <- diag(3)
            .Object$P2 <- matrix(1,3,3)
            .Object$TRpars <- c(2.0,0.5)  #speed = 2, location = 0.5
            .Object$shape <- TRshape$single
            .Object$speedoption <- TRspeedopt$eta    
            .Object$optimcontrol <- list(fnscale = -1)
            # Return:
            .Object
          })

setGeneric(name="newSTCC",def=function(Tobs,st,type,p1,p2,trpars,shape,speedopt,...){
  standardGeneric("newSTCC")
}, valueClass = "stcc_class", signature = c("Tobs","st","type","p1","p2","trpars","shape","speedopt"))


# Create default constructor:
sig1 <- c("missing","numeric","numeric","matrix","matrix","numeric","numeric","numeric")
setMethod(f="newSTCC",signature = sig1,
          definition = function(st,type,p1,p2,trpars,shape,speedopt, ...){
            # Validate params
            if (is.null(st)) stop("Error: missing 'st'\nAll parameters are required to construct an STCC object")
            if (is.null(type)) stop("Error: missing 'type'\nAll parameters are required to construct an STCC object")
            if (is.null(p1)) stop("Error: missing 'p1'\nAll parameters are required to construct an STCC object")
            if (is.null(p2)) stop("Error: missing 'p2'\nAll parameters are required to construct an STCC object")
            if (is.null(trpars)) stop("Error: missing 'trpars'\nAll parameters are required to construct an STCC object")
            if (is.null(shape)) stop("Error: missing 'shape'\nAll parameters are required to construct an STCC object")
            if (is.null(speedopt)) stop("Error: missing 'speedopt'\nAll parameters are required to construct an STCC object")
            #
            stcc <- new("stcc_class") 
            stcc$st <- st
            stcc$Tobs <- NROW(st)
            stcc$type <- type
            stcc$P1 <- p1
            stcc$P2 <- p2
            stcc$TRpars <- trpars
            stcc$shape <- shape
            stcc$speedoption <- speedopt
            # Return:
            stcc
          })

# Create convenience 'quick' constructor: Linear Transition from 0 to 1, using eta
sig2 <- c("numeric","missing","missing","missing","missing","missing","missing","missing")
setMethod(f="newSTCC",signature = sig2,
          definition = function(Tobs){
            # Validate params
            if (is.null(Tobs)) stop("Error: missing 'Tobs'\nPlease enter a valid number of observations to create an STCC object")
            #
            stcc <- new("stcc_class") 
            stcc$Tobs <- Tobs
            stcc$st <- 1:Tobs/Tobs
            # Return:
            stcc
          })



#### ===============  Utility Functions  ================== ####

beep <- function(n = 2){
  for(i in seq(n)){
    system("rundll32 user32.dll,MessageBeep -1")
    Sys.sleep(.5)
  }
}

constructTVpars <- function(pars,shape) {
## This function will construct the TV$pars matrix required by the code
## Inputs:  pars  - vector of TV pars in form: d1,s1,l1.1,<loc1.2>,d2,s2,l2.1,<loc2.2>...
##                  ***  Note 1: delta0 is not passed in here!! ***
  
##                  Note 2: <locn.2> is only included where the shape = 2.  No need to add NA
##          shape - vector describing the transition type for each transition, 
##                  e.g. c(TRshape$single, TRshape$double)
##                  Note: TV$shape should be set to TRShape$none if there are no transitions (i.e. delta0 only)
  
  if(shape[1]==TRshape$none) stop("Can't construct TV$pars with no shape!!")
  # 1. Check that the pars & shape match
  expectedNumPars <- sum(NROW(shape)*3) + sum(NROW(shape[which(shape==TRshape$double)])) 
  if (NROW(pars)!=expectedNumPars) stop("pars & shape don't match\n 'Make sure there are 2 locations for all quadratic transitions and only 1 location for all others\n  No NA's should be included in pars'")
  
  # 2. Now add NA's for all missing locn.2 pars:
  naPars <- NULL
  for (i in seq_along(shape)) {
    if (shape[i] == TRshape$double) {
      naPars <- c(naPars,pars[1:4])
      pars <- pars[-(1:4)]
    } else {
      naPars <- c(naPars,pars[1:3],NA)
      pars <- pars[-(1:3)]
    }
  }
  # Return:
  matrix(naPars,nrow=4,ncol=NROW(shape),dimnames=list(c("delta","speed","loc1","loc2"),NULL))
  
} # End: constructTVpars()

getEstimatedTVpars <- function(tv,pars) {
## This function will write the values contained in the 'pars' parameter
## into the tv$Estimated list, and return the updated tv object.
## Input:  tv   - a tv object
##         pars - vector of TV pars in form: <d0>,d1,s1,l1.1,<loc1.2>,d2,s2,l2.1,<loc2.2>...
##                Note: <locn.2> is only included where the shape = 2.  No need to add NA
  
  # Handle delta0:
  if (tv$var_target) tv$Estimated$delta0 <- pars[1] else tv$Estimated$delta0 <- tv$Estimated$lastdelta0
  names(tv$Estimated$delta0) <- "delta0"
  
  # Now extract the transition parameters from optimpars back into the tv() object
  if (tv$shape[1] > TRshape$none) {
    numTVpars <- length(tv$pars[!is.na(tv$pars)])
    startPos <- endPos <- 0
    
    if (tv$var_target) {
      startPos <- 2 
      endPos <- numTVpars + 1
      } else {
        startPos <- 1
        endPos <- numTVpars
    }
    tv$Estimated$pars <- constructTVpars(pars[startPos:endPos],tv$shape)
  }
  
  # Now extract any linear terms:
  if(NROW(tv$linpars) > 1) tv$Estimated$linpars <- tail(pars,length(tv$linpars[!is.na(tv$linpars)]))
  
  # Return:
  tv
  
} # End: getEstimatedTVpars()

getTVparsVector <- function(m_tvPars) {
##  Return the matrix of tv params as a column-based vector, with all NA's stripped out.
  if (is.na(m_tvPars[1])) return(NA)
  v_tvPars <- as.vector(m_tvPars)
  # Set the parameter names
   
  namesVec <- vector("character")
  for (i in 1:NCOL(m_tvPars)) {
    namesVec <- c(namesVec,paste0("delta",i),paste0("speed",i),paste0("loc",i,".1"),paste0("loc",i,".2"))
  }
  names(v_tvPars) <- namesVec
  
  # Return:
  round(v_tvPars[!is.na(v_tvPars)],6)
}

constructGARCHpars <- function(garch,pars) {
##  This function ensures that the Garch Type & parameters match
##  It also takes care of Variance Targeting
  
  garchPars <- NULL
  
  if (garch$type == GARCHtype$constant) garchPars <- c(pars[1],0,0,0)
  else if (garch$type == GARCHtype$general) {
    if (garch$var_target) garchPars <- c((1-pars[1]-pars[2]),pars[1],pars[2],0) 
    else garchPars <- c(pars[1],pars[2],pars[3],0)
  }
  else if (garch$type == GARCHtype$GJR) {
    if (garch$var_target) garchPars <- c((1-pars[1]-pars[2]-pars[3]/2),pars[1],pars[2],pars[3]) 
    else garchPars <- c(pars[1],pars[2],pars[3],pars[4])
  }
  else if (garch$type == GARCHtype$GJR_alpha0) {
    if (garch$var_target) garchPars <- c((1-pars[1]-pars[2]/2),0,pars[1],pars[2]) 
    else garchPars <- c(pars[1],0,pars[2],pars[3])
  } 
  else stop("Unrecognised Garch Type!")
  
  # Return as a single-column named matrix:
  garchPars <- as.matrix(garchPars) 
  names(garchPars) = c("omega","alpha","beta","gamma")
  # Return:
  garchPars
  
}

matrix.sortbyCol <- function (mat, sortbyColumnNr, increasing = TRUE) {
  # This function will sort a matrix by a column number
  m <- do.call("order", c(as.data.frame(mat[, sortbyColumnNr]), decreasing = !increasing))
  mat[m, ]
}

unregisterDoParallel <- function(cluster=NULL) {
  # This function will release all connections created by the cluster
  # Very useful when testing, as the environment builds up quickly!
  if(!is.null(cluster)) stopCluster(cluster)
  stopImplicitCluster()
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

calculateOmega <- function(alpha,beta,gamma=0) {
  # When VarienceTargetting is ON: omega is calculated as below:
  # Since the 'gamma' param is only used for GJR Garch types, we default it to 0
  #Return:
  omega <- 1-alpha-beta-0.5*gamma
}

calculate_local_variance <- function(e,localvarwindow=500) {
  # Local Varience Targetting is ON: We can estimate Garch parameters, while minimising the effect of shock-persistence on alpha & beta
  # omega can then be calculated as: localVar * (1-alpha-beta-0.5*delta)
  
  Tobs <- length(e)
  localVar <- vector("numeric",Tobs)
  
  for (t in 1:Tobs) {
    if (t < (localvarwindow/2)) {
      localVar[t] <- var(e[1:localvarwindow])
    } else if (t > (Tobs - localvarwindow/2)) {
      localVar[t] <- var(e[(Tobs - localvarwindow):Tobs])
    } else localVar[t] <- var(e[(t - localvarwindow/2):(t + localvarwindow/2)])
  }
    
  # Return:
  localVar
}

calculate_local_kurtosis <- function(e,localkurtwindow=500) {
  Tobs <- length(e)
  localKurt <- vector("numeric",Tobs)
  
  for (t in 1:Tobs) {
    if (t < (localkurtwindow/2)) {
      localKurt[t] <- kurtosis(e[1:localkurtwindow])
    } else if (t > (Tobs - localkurtwindow/2)) {
      localKurt[t] <- kurtosis(e[(Tobs - localkurtwindow):Tobs])
    } else localKurt[t] <- kurtosis(e[(t - localkurtwindow/2):(t + localkurtwindow/2)])
  }
  
  # Return:
  localKurt
}


#### ===============  Utility Functions  ================== ###


#### ===============  Anna's Tricky Matrix Functions  ================== ####


myVecl <- function(M){
  ## Returns the lower triangle of a square matrix in vector format.
  ## Note: This operation can be reversed using myUnVecl()
  idx = 0
  N <- ncol(M)
  vM <- matrix(0,N*(N-1)/2,1)
  for (idxCol in seq(1,N-1)){
    for (idxRow in seq(idxCol+1,N)){
      idx <- idx+1
      vM[idx,1] <- M[idxRow,idxCol]
    }
  }
  #Return:
  vM
}

myUnVecl <- function(vM){
  ## Returns a square matrix constructed using a vector of it's lower triangle.
  ## Note: This operation can be reversed using myVecl()
  k <- length(vM)
  N <- (1+sqrt(1+8*k))/2
  M <- matrix(0,N,N)
  M[lower.tri(M,diag=FALSE)] <- vM
  M <- M + t(M) + diag(N)
  #Return:
  M
}


myQ.EC <- function(N)
{
  # compute matrix of eigenvectors (columns of Q) for an EQUI-correlation model
  # N = number of series
  Q <- matrix(0,N,N)
  for (i in 1:N)
  {
    if (i==1) Q[,i] <- (1/sqrt(N))*matrix(1,N,1)
    else
    {
      tmp <- 1/(sqrt((N-i+2)*(N-i+1)))
      for (j in seq((i-1),N))
      {
        if (j==(i-1)) Q[j,i] <- tmp*(N-i+1)
        else Q[j,i] <- tmp*(-1)
      }
    }
  }
  #Return
  Q # NXN
}  #End: myQ.EC(N)

myL.EC <- function(N,rho){
  # compute matrix of eigenvvalues for an EQUI-correlation model
  # N = number of series
  # rho = equicorrelation parameter (scalar)
  L <- rep((1-rho),N)
  L[1] <- L[1]+rho*N
  #Return vector:
  L # Nx1
}

myFilter <- function(mX,vB){
  # mX -- T x s matrix
  # vB -- 1 x s vector of coefficients
  # does AR type filtering with lag 1 only
  # output -- mY -- T x s matrix
  # mY[1,] = 0...0
  # mY[t,s] = mX[t,s]+vB[s]*mY[t-1,s]
  s <- length(vB)
  Tobs <- NROW(mX)
  # GLEN: add error check, NCOL(mX)=s
  mY <- matrix(mX[1,],1,s)
  for (t in 2:Tobs){
    mY <- rbind(mY,mX[t,]+vB*mY[(t-1),])
  }
  # Return:
  mY
}

myVecd <- function(M){
  # vectorises main diagonal and off diagonal
  N <- NCOL(M)
  if (N==3){
    vM <- matrix(0,nrow=N*(N+1)/2,ncol=1)
    vM[1:3] <- diag(M)^2
    vM[4] <- sqrt(2)*M[3,2]
    vM[5] <- sqrt(2)*M[3,1]
    vM[6] <- sqrt(2)*M[2,1]
    return(vM)
  }
  else return(NULL)
  
}

#### ===============  Anna's Tricky Functions  ================== ###


#### ====== Derivative Functions & 'g' Calculators ======= ####

myG <- function(speed,loc1,loc2,st,shape,speedoption){
  # input:	speed,loc1,<loc2>
  #   			st   = transition variable
  #         shape = TRshape: none,linear,quadratic,cubic
  #         speedoption = 1: gamma, 2: gamma/sd(st), 3: exp(eta)
  # output:	Tx1 matrix
  
  st_c <- switch(shape,st-loc1,(st-loc1)*(st-loc2),(st-loc1)*(st-loc1))
  G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
  return(1/G_inv)
  
}

calculate_g <- function(tv) {
  ###  Calculate the Conditional Variances of a TV model
  ###  Input: TV object
  ###  Output: column-based vector 'g', of Length = number of observations
  
  # 1. Initialise g to a constant variance = delta0
  if(is.null(tv$Estimated$delta0)) tv$Estimated$delta0 <- tv$Estimated$lastdelta0 <- tv$delta0
  if (tv$var_target) g <- rep(tv$Estimated$delta0,tv$Tobs) else  g <- rep(tv$Estimated$lastdelta0,tv$Tobs)
  
  
  # 2. Update based on any transition parameters in the model
  if (tv$shape[1] > TRshape$none) {
    # initialise some variables
    st_c <- 0
    Gi <- 0
    # calulate 'g'
    for (i in seq_along(tv$shape)) {
        st_c <- switch(tv$shape[i],tv$st-tv$Estimated$pars["loc1",i],(tv$st-tv$Estimated$pars["loc1",i])*(tv$st-tv$Estimated$pars["loc2",i]),(tv$st-tv$Estimated$pars["loc1",i])^2)  
        switch (tv$speedoption,
          Gi <- 1/(1+exp(-tv$Estimated$pars["speed",i] * st_c)),
          Gi <- 1/(1+exp(-tv$Estimated$pars["speed",i] * st_c/sd(tv$st))),
          Gi <- 1/(1+exp(-exp(tv$Estimated$pars["speed",i]) * st_c)) 
        ) #end: switch()
        g <- g + tv$Estimated$pars["delta",i]*Gi
    }
  }
  
  # 3. Update based on any linearised parameters in the model
  if (NROW(tv$linpars) > 1) {
    #if (!is.na(tv$linpars[1])) g <- g + tv$Estimated$linpars[1]
    if (!is.na(tv$linpars[2])) g <- g + tv$Estimated$linpars[1]*tv$st
    if (!is.na(tv$linpars[3])) g <- g + tv$Estimated$linpars[2]*tv$st^2
    if (!is.na(tv$linpars[4])) g <- g + tv$Estimated$linpars[3]*tv$st^3
  }
  #Return:
  g
}  # End: calculate_g()

dg.dt <- function(tv){
  # derivative of g(t) w.r.t. TV-parameters (d0,d1,eta,c1,<c2>,...,<lin1>,<lin2>,...)
  # output: T x (#of TV pars (incl. delta0) in H0) matrix
  Tobs <- length(tv$condvars)
  # 1. Process TV$delta0 + TV$pars - We will process TV$linpars(if any) later
  if (tv$shape[1]==TRshape$none){
    rtn <- matrix(nrow=Tobs,ncol=1) 
  } else {
    numTVpars <- length(getTVparsVector(tv$pars))
    rtn <- matrix(nrow=Tobs,ncol=numTVpars+1)
  }
  col_idx <- 1
  rtn[,col_idx] <- 1  # derivative of delta0
  
  if (tv$shape[1] > TRshape$none) {
    # initialise some variables
    stdev_st <- sd(tv$st)
    st_c <- speed_transf <- Gi <- 0
    
      for (i in seq_along(tv$shape)) {
    
        # Define: Transition Variable-Location, st-c = 'st_c'
        st_c <- switch(tv$shape[i],tv$st-tv$Estimated$pars["loc1",i],(tv$st-tv$Estimated$pars["loc1",i])*(tv$st-tv$Estimated$pars["loc2",i]),(tv$st-tv$Estimated$pars["loc1",i])^2)
        
        switch (tv$speedoption,
          {speed_transf <- tv$Estimated$pars["speed",i]
           Gi <- 1/(1+exp(-tv$Estimated$pars["speed",i] * st_c))},
          
          {speed_transf <- tv$Estimated$pars["speed",i]/stdev_st
           Gi <- 1/(1+exp(-tv$Estimated$pars["speed",i]/stdev_st*st_c))},
          
          {speed_transf <- exp(tv$Estimated$pars["speed",i])
           Gi <- 1/(1+exp(-exp(tv$Estimated$pars["speed",i])*st_c))}
        ) #end: switch()
  
        deriv_const <- tv$Estimated$pars["delta",i]*speed_transf*Gi*(1-Gi)
        
        col_idx <- col_idx + 1
        rtn[,col_idx] <- Gi    # derivative of delta1..n
        col_idx <- col_idx + 1
        rtn[,col_idx] <- deriv_const*st_c    # derivative of speed1..n
        
        switch (tv$shape[i],
                {
                  col_idx <- col_idx + 1
                  rtn[,col_idx] <- -deriv_const    # derivative of loc1..n (shape=TRshape$single)        
                },
                {
                  col_idx <- col_idx + 1
                  rtn[,col_idx] <- -deriv_const*(tv$st-tv$Estimated$pars["loc2",i])  # derivative of loc1..n (shape=TRshape$double)
                  col_idx <- col_idx + 1
                  rtn[,col_idx] <- -deriv_const*(tv$st-tv$Estimated$pars["loc1",i])  # derivative of loc2..n (shape=TRshape$double)
                },
                {
                  col_idx <- col_idx + 1
                  rtn[,col_idx] <- -deriv_const*2*(tv$st-tv$Estimated$pars["loc1",i])    # derivative of loc1..n (shape=TRshape$triple)
                }
        )
      } # End: for (i in seq_along(shape))

  } # eND: if (is.na(tv$shape[1])) 
  
  # 2. Now calculate the derivatives of the linearised parameters (if any):
  if (NROW(tv$linpars) > 1) {
    #if (!is.na(tv$linpars[1])) rtn <- cbind(rtn,1)  # Skip this, we already have a constant
    if (!is.na(tv$linpars[2])) rtn <- cbind(rtn,tv$st)
    if (!is.na(tv$linpars[3])) rtn <- cbind(rtn,tv$st^2)
    if (!is.na(tv$linpars[4])) rtn <- cbind(rtn,tv$st^3)
  }
  
  #Return:
  dimnames(rtn) <- NULL
  rtn  # T x (#of TV pars incl. delta0)
  
}  # End dg.dt

dg.dt2.H0 <- function(st){
  # purpose:	recursively calculates the partial derivative dg/dTheta2 in the linearized model (Max 5 linpars)
  #	Theta2 = parameters from the linearized TV component, d0,d1,d2,d3,d4
  # output:	Tx4 matrix, each row = dg/dd1, dg/dd2, dg/dd3, dg/dd4 
  # Note: dg/dd0 is always 1 so we skip it!

  rtn <- matrix(1,NROW(st),4)
  rtn[,1] <- st
  rtn[,2] <- st^2
  rtn[,3] <- st^3
  rtn[,4] <- st^4
  
  #Return:
  rtn # Tx4
}

#### ====== Derivative Functions & 'g' Calculators ======= ###


#### ================== Test Functions ================== ####

myTest.TV.noGARCH.TR2 <- function(e,tv,testorder=TESTorder$H0,debugmode=FALSE){
  # TR2-form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # H0: gt=delta0
  # H1: gt=delta0 + TV$pars
  # e    = raw data
  # shape = shape of TV under H0
  # speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # st = transition variable in TV and used in test
  
  #H03 Test => TV$linpars must be c(NA,1,1)
  #H02 Test => TV$linpars must be c(NA,1)
  #H01 Test => TV$linpars must be c(NA)
  
  TestStat <- NA
  
  # Test Method: Regress psi2_1 on 1/gt*(dgdt and dgdt2.H0) (here t=TV component in H0, t2=linearised (tested) TV component)
  
  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or... NCOL(dgdt) increases with the order of TV function. 
  dgdt <- dg.dt(tv)
  
  # 2. Calc derivatives of t2 (linearised component) under the null
  if (testorder==TESTorder$H0) {dgdt2 <- dg.dt2.H0(tv$st)}
  else if (testorder==TESTorder$H01) {dgdt2 <- as.matrix(tv$st)}
  else if (testorder==TESTorder$H02) {dgdt2 <- as.matrix((tv$st)^2)}
  else if (testorder==TESTorder$H03) {dgdt2 <- as.matrix((tv$st)^3)}
  else return(NA)
  
  g <- calculate_g(tv)
  X <- cbind(dgdt,dgdt2)/g
  
  # 3. Invert crossprod(X) to calculate SSR1
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.TR2() threw an error doing solve() - returning NA")
    rm(g,dgdt,dgdt2)
    if (debugmode) {
      # Return the X & X'X matrices for investigation!
      TestStat <- list()
      TestStat$X <- X
      TestStat$XX <- crossprod(X)
      return(TestStat)
    } else return(NA)
  }
  
  # 4. Calculate psi2_1 to calculate SSR0
  psi2_1 <- matrix(data=(e*e/g-1),nrow = tv$Tobs,ncol = 1)
  
  # 5. Calc the TestStat:
  SSR0 <- sum(psi2_1*psi2_1)    # Scalar
  SSR1 <- sum((psi2_1-X%*%Xinv%*%t(X)%*%psi2_1)^2)
  # Return:
  TestStat <- tv$Tobs*(SSR0-SSR1)/SSR0
  
}  #End: myTest.TV.noGARCH.TR2a()

myTest.TV.noGARCH.robust <- function(e,tv,testorder=TESTorder$H0,debugmode=FALSE){
  # Robust form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # Ho: gt=delta0
  # H1: gt=delta0 + TV-terms
  # pars = estimate of delta0
  # e    = raw data
  # tv$shape = shape of TV under H0
  # tv$speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # tv$st = transition variable in TV and used in test
  
  #H03 Test => TV$linpars must be c(NA,1,1)
  #H02 Test => TV$linpars must be c(NA,1)
  #H01 Test => TV$linpars must be c(NA)
  
  TestStat <- NA
  
  # Test Method:  Get squared standardised residuals minus one (psi2_1)
  # Regress r2 = 1/gt*dgdt2.Ho on r1 = 1/gt*dgdt, and compute residual vectors w(t)

  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or... NCOL(dgdt) increases with the order of TV function. 
  dgdt <- dg.dt(tv)
  
  # 2. Calc derivatives of t2 (linearised component) under the null
  if (testorder==TESTorder$H0) {dgdt2 <- dg.dt2.H0(tv$st)}
  else if (testorder==TESTorder$H01) {dgdt2 <- as.matrix(tv$st)}
  else if (testorder==TESTorder$H02) {dgdt2 <- as.matrix((tv$st)^2)}
  else if (testorder==TESTorder$H03) {dgdt2 <- as.matrix((tv$st)^3)}
  else return(err_output)
  
  g <- calculate_g(tv)
  X <- dgdt/g  
  
  # 3. Invert crossprod(X) 
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.robust() threw an error doing solve() - returning NA")
    rm(g,X,dgdt,dgdt2)
    if (debugmode) {
      # Return the X & X'X matrices for investigation!
      TestStat <- list()
      TestStat$X <- X
      TestStat$XX <- crossprod(X)
      return(TestStat)
    } else return(NA)
  }
  XXXX <- X%*%Xinv%*%t(X)
  Y <- as.matrix(dgdt2/g)
  W <- as.matrix(Y-XXXX%*%Y)

  #3. Regress 1 on (psi2-1)*w, and compute SSR
  psi2_1 <- as.vector(e*e/g-1)
  X <- psi2_1*W  #psi2_1 must be a vector for this!!
  
  #4. Compute test statistic:
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.Robust() threw an error doing solve() - returning NA")
    rm(psi2_1,g,W,Y,dgdt,dgdt2)
    if (debugmode) {
      # Return the X & X'X matrices for investigation!
      TestStat <- list()
      TestStat$X <- X
      TestStat$XX <- crossprod(X)
      return(TestStat)
    } else return(NA)
  }

  TestStat <- tv$Tobs-sum(diag(tv$Tobs)-(X%*%Xinv%*%t(X)))
  
  #Tidy up & release memory before returning:
  rm(tv,psi2_1,g,W,Y,X,Xinv,dgdt,dgdt2)
  
  #Return:
  TestStat
  
}  #End: myTest.TV.noGARCH.robust()

myTest.CCCvSTCC.LM <- function(e,H0,H1,testorder=1) {
  
  nGARCH <- H0$nGARCH
  nTV <- H0$nTV
  CCC <- H0$CCC 
  STCC <- H1
  
  if (is.null(CCC)) stop("LM Test: need CCC as input (H0)")
  if (is.null(STCC$st)) stop("LM Test: need transition variable for the test as input (H1)")
  
  N <- NCOL(e)  # Number of series
  Tobs <- NROW(e)
  P <- CCC$P
  st <- STCC$st
  testOrder <- testorder    # Can be 1 or 2, indicating single or double transition
  I <- diag(nrow = N,ncol = N) # identity matrix
  Pinv <- try(solve(P))
  if (!is.matrix(Pinv)) stop("LM Test: Can't invert the correlation matrix P")
  
  #Construct the K matrix: N^2 x N^2
  if(T){
    K <- NULL
    for (i in 1:N) {
      # block rows
      Krow <- NULL
      for (j in 1:N) {
        # block columns
        block <- matrix(0,N,N)
        block[j,i] <- 1
        Krow <- cbind(Krow,block)
      }
      K <- rbind(K,Krow)
    }
  }
  
  #Construct the U matrix:Dimensions = N^2 x N(N-1)/2
  if(T){
    U <- NULL
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        block <- matrix(0,N,N)
        block[i,j] <- block[j,i] <- 1
        Ucol <- as.vector(block)
        U <- cbind(U,Ucol)
      }
    }
  }
  
  # 2. Create the Information Matrix for the Correlation, and extract the SE corner of its inverse:
  if(T){
    
    # Initialise variables:
    IM_cor <- matrix(0,(testOrder+1)*N*(N-1)/2,(testOrder+1)*N*(N-1)/2)
    g <- matrix(1,Tobs,N) # for now, no TV estimated, so these are 1, T x N
    h <- matrix(1,Tobs,N) # place for estimated GARCH variances, ones at this stage, T x N
    if (!is.null(nGARCH)){
      # Extract the number of Garch pars: (Assume that Garch$type=0 has 1 par only, Garch$type=1 has 3, etc. )
      Num_garch_pars <- 0
      for (n in 1:N) Num_garch_pars <- Num_garch_pars + length(nGARCH[[n]]$pars)
      IM_garch <- matrix(0,Num_garch_pars,Num_garch_pars)
      IM_garch_cor <- matrix(0,Num_garch_pars,(testOrder+1)*N*(N-1)/2) 
    }
    
    # Calculate real values
    if (!is.null(nGARCH)){
      beta <- rep(0,N) # estimated beta coefficients from GARCH equations
      for (n in 1:N){
        h[,n] <- nGARCH[[n]]$condvars # univariate GARCH volatilities, T x N 
        if(nGARCH[[n]]$type==0) beta[n] <- 0 else beta[n] <- nGARCH[[n]]$pars[3]
      }
    }
    if (!is.null(nTV)){
      for (n in 1:N){
        g[,n] <- nTV[[n]]$condvars # univariate TV volatilities, T x N 
      }
    }
    w <- e/sqrt(g) # returns standardised by g(t), T x N
    z <- w/sqrt(h) # returns standardised by GARCH volatilities and TV g(t), T x N

    # partial derivatives of P w.r.t correlation parameters, T x 1 or T x 2
    # these are for alternative parameters, later add column of ones at the front
    if (testOrder==1) v_rho <- cbind(st) else if (testorder==2) v_rho <- cbind(st,st^2) else stop("LM Test: testOrder must be 1 or 2!")
    
    # score for rho
    zKRONz <- NULL
    for (t in 1:Tobs) zKRONz <- cbind(zKRONz, z[t,] %x% z[t,]) # (N^2 x T), each col = "z_t kron z_t"
    
    dlldrho_A <- -0.5*t(U)%*%( as.vector(Pinv)%x%t(rep(1,Tobs))-(Pinv%x%Pinv)%*%(zKRONz)  )%*%v_rho # N(N-1)/2 x testorder
    dlldrho_A <- vec(dlldrho_A) # testorder*N(N-1)/2 x 1, SUM OVER TIME
    v_rho <- cbind(rep(1,Tobs),v_rho) # T x 2 or T x 3, now add column of ones at the front
    
    # partial derivatives for GARCH and "x" matrices for GARCH parts    
    if (!is.null(nGARCH)){
      # partial derivatives of h1,h2,...,hN w.r.t garch_pars
      v_garch <- NULL
      h_scale <- NULL
      beta_scale <- NULL
      for (n in 1:N){
        if(nGARCH[[n]]$type==GARCHtype$constant) {
            # constant GARCH case
            v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1)))) 
            beta_scale <- c(beta_scale,beta[n])
            h_scale <- cbind(h_scale,h[,n])
          } else if(nGARCH[[n]]$type==GARCHtype$general) {
            # standard GARCH(1,1) case
            v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1),n]^2),c(0,h[1:(Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N    
            beta_scale <- beta_scale <- c(beta_scale,(beta[n] %x% rep(1,3)))
            h_scale <- cbind(h_scale,h[,n] %x% t(rep(1,3)))
          } else if (nGARCH[[n]]$type>GARCHtype$general) {
            # GJR or other forms not done yet
            print("ARRGH!")
          }
        #GLEN: the above needs an error code for garchtype>general
      }
      dhdt <- myFilter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
      # matrix containing the x_t's (garch part)
      x_garch <- -0.5*dhdt/(h_scale) # T x Num_garch_pars, each row = "x_it'", i=1,...,N
    }
    
    # partial derivatives for TV and "x" matrices for TV parts    
    if (!is.null(nTV)){
      v_tv <- NULL
      h_scale <- NULL
      g_scale <- NULL
      beta_scale <- NULL
      dgdt <- NULL
      dhdt <- NULL
      # partial derivatives of g1,g2,...,gN w.r.t tv_pars
      for (n in 1:N){
        dgdt <- dg.dt(nTV[[n]])  # T x #of TVpars in TV[n], includes d0's derivative
        if (nGARCH[[n]]$type>GARCHtype$none) dgdt<-dgdt[,2:NCOL(dgdt)] # drop delta0 if there is GARCH of any form
        g_scale <- cbind(g_scale,g[,n] %x% t(rep(1,NCOL(dgdt)))) # cbinds g(nt) as many times as g(n) has tv parameters
        if (!is.null(nGARCH)){
          if (nGARCH[[n]]$type==GARCHtype$general) {
            # standard GARCH(1,1) case
            v_tv <- cbind(v_tv,-nGARCH[[n]]$pars[2]*c(0,1/g[1:(Tobs-1)])*c(0,w[1:(Tobs-1)]^2)*dgdt)
            beta_scale <- c(beta_scale,(beta[n] %x% rep(1,NCOL(dgdt))))
            h_scale <- cbind(h_scale,h[,n] %x% t(rep(1,NCOL(dgdt))))
          } else if(nGARCH[[n]]$type>GARCHtype$general) {
            # GJR or other type not done yet!
            print("ARRGH2!")
          } else {
            # if no GARCH or constant only GARCH, fill v_tv and beta_scale with zeros
            v_tv <- cbind(v_tv,0 %x% dgdt)
            beta_scale <- c(beta_scale,0 %x% rep(1,NCOL(dgdt)))
            h_scale <- cbind(h_scale,h[,n] %x% t(rep(1,NCOL(dgdt))))
          }
        }
      }
      # THIS is needed only if there is GARCH (constant garch spaces filled with zeros)
      if (!is.null(nGARCH)){
        dhdt <- myFilter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
      }
  
      # matrix containing the x_t's (tv part)
      if (is.null(nGARCH)){
        x_tv <- -0.5*dgdt/(g_scale)  # T x Num_tv_pars, each row = x_it, i=1...N 
      } else {
        x_tv <- -0.5*dhdt/(h_scale)-0.5*dgdt/(g_scale)  # T x Num_tv_pars, each row = x_it, i=1...N 
      }
    }
    
    
    # Information matrix:
    # Blocks involving GARCH only or GARCH and Correlation:
    if (!is.null(nGARCH)){
      # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
      One_33 <- matrix(1,3,3)
      One_31 <- matrix(1,3,1)
      One_13 <- matrix(1,1,3)
      One_11 <- matrix(1,1,1)
      I.P.Pinv <- I + P*Pinv
      I.P.Pinv_scale <- NULL
      scaleFactor <- NULL
      for (i in 1:N) {
        # i = row index
        I.P.Pinv_scale_row <- NULL
        for (j in 1:N) {
          # j = col index
          if(nGarch[[i]]$type==1 && nGarch[[j]]$type==1) scaleFactor <- One_33 else
            if(nGarch[[i]]$type==1 && nGarch[[j]]$type==0) scaleFactor <- One_31 else
              if(nGarch[[i]]$type==0 && nGarch[[j]]$type==1) scaleFactor <- One_13 else
                if(nGarch[[i]]$type==0 && nGarch[[j]]$type==0) scaleFactor <- One_11 else stop("Something went wrong!")
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
        I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
      }
      
      IM_garch <- ((t(x_garch)%*%x_garch) * I.P.Pinv_scale) / Tobs
    
      # IM_garch_cor, 3N x (testorder+1)*N(N-1)/2
      Mhelp1 <- NULL
      for (n in 1:N){
        Mhelp1 <- rbind(Mhelp1, (Pinv[n,]%x%I[n,]+I[n,]%x%Pinv[n,])) # N x N^2
      }
      Mhelp1 <- -0.5*(Mhelp1%*%U)
      Mhelp1_scale <- NULL  # (Num_Garch_Pars x N(N-1)/2)
      for (n in 1:N) {
        if (nGarch[[n]]$type==0) Mhelp1_scale <- rbind(Mhelp1_scale,Mhelp1[n,]) else
          if (nGarch[[n]]$type==1) Mhelp1_scale <- rbind(Mhelp1_scale, (Mhelp1[n,] %x% One_31))
      }

      Mhelp2 <- t(t(v_rho)%*%x_garch)/Tobs # Num_garch_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
      IM_garch_cor <- NULL
      for (i in 1:NCOL(v_rho)){
        IM_garch_cor <- cbind(IM_garch_cor,((Mhelp2[,i]%x%t(rep(1,(N*(N-1))/2))) * Mhelp1_scale))  # (Num_Garch_Pars x (testOrder+1)*N(N-1)/2),  SUM OVER TIME
      }
    }
    
    # Blocks involving TV only or TV and Correlation:
    if (!is.null(nTV)){
      # to be written...
      IM_tv <- NULL
      I.P.Pinv <- I + P*Pinv
      I.P.Pinv_scale <- NULL
      scaleFactor <- NULL
      for (i in 1:N) {
        # i = row index
        ii <- length(nTV[[i]]$pars[!is.na(nTV[[i]]$pars)])
        if (nGarch[[i]]$type==GARCHtype$none) ii <- ii+1   #add delta0 to the count
        I.P.Pinv_scale_row <- NULL
        for (j in 1:N) {
          # j = col index
          jj <- length(nTV[[j]]$pars[!is.na(nTV[[j]]$pars)])
          if (nGarch[[j]]$type==GARCHtype$none) jj <- jj+1   #add delta0 to the count
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
        I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
      }
      IM_tv <- ((t(x_tv)%*%x_tv) * I.P.Pinv_scale) / Tobs  # sum(num_tv_pars) x sum(num_tv_pars)
      
      # IM_tv_cor, (#of tv pars in n,n=1...N)N x (testorder+1)*N(N-1)/2
      Mhelp1 <- NULL
      for (n in 1:N){
        Mhelp1 <- rbind(Mhelp1, (Pinv[n,]%x%I[n,]+I[n,]%x%Pinv[n,])) # N x N^2
      }
      Mhelp1 <- -0.5*(Mhelp1%*%U)  # N x N(N-1)/2
      Mhelp1_scale <- NULL  # (Total Num_tv_Pars x N(N-1)/2)
      for (n in 1:N) {
        ii <- length(nTV[[n]]$pars[!is.na(nTV[[n]]$pars)])
        if (nGarch[[n]]$type==GARCHtype$none) ii <- ii+1   #add delta0 to the count
        scaleFactor <- matrix(1,ii,1)
        Mhelp1_scale <- rbind(Mhelp1_scale, (Mhelp1[n,] %x% scaleFactor))
      }
      Mhelp2 <- t(t(v_rho)%*%x_tv)/Tobs # Num_tv_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
      IM_tv_cor <- NULL
      for (i in 1:NCOL(v_rho)){
        IM_tv_cor <- cbind(IM_tv_cor,((Mhelp2[,i]%x%t(rep(1,(N*(N-1))/2))) * Mhelp1_scale))  # (Num_tv_Pars x (testOrder+1)*N(N-1)/2),  SUM OVER TIME
      }
    }
    # Blocks involving TV and GARCH:
    if (!is.null(nGARCH) && !is.null(nTV)){
      I.P.Pinv <- I + P*Pinv  # N x N
      I.P.Pinv_scale <- NULL
      for (i in 1:N){
        # i = row index
        ii <- length(nTV[[i]]$pars[!is.na(nTV[[i]]$pars)])
        if (nGarch[[i]]$type==GARCHtype$none) ii <- ii+1   #add delta0 to the count
        I.P.Pinv_scale_row <- NULL
        for (j in 1:N){
          # j = col index
          jj <- length(nGARCH[[j]]$pars[!is.na(nGARCH[[j]]$pars)])
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
        I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
      }
      
      IM_tv_garch <- NULL
      IM_tv_garch <- (t(x_tv)%*%x_garch)*(I.P.Pinv%x%scaleFactor)/Tobs  #  (total # of tvpars)x(total # of garchpars)
    }
    # Block involving Correlation only:
    # IM_cor (testorder+1)*N(N-1)/2, SUM OVER TIME
    Mhelp3 <- t(U)%*%(Pinv%x%Pinv+(Pinv%x%I)%*%K%*%(Pinv%x%I))%*%U # N(N-1)/2 x N(N-1)/2
    Mhelp4 <- t(v_rho)%*%v_rho/Tobs # 2x2 or 3x3, SUM OVER TIME
    IM_cor <- 0.25*(Mhelp4%x%Mhelp3) # (testorder+1)*N(N-1)/2 x (testorder+1)*N(N-1)/2
    
    # IM
    #IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
    # block corresponding to correlations of the inverse of the IM matrix:
    if (!is.null(nGARCH) && !is.null(nTV)){
      IM_tv_inv <- try(solve(IM_tv))
      if (!is.matrix(IM_tv_inv)) stop("LM Test: Can't invert the information matrix IM_tv")
      MHelp <- IM_garch - t(IM_tv_garch)%*%IM_tv_inv%*%IM_tv_garch
      MHelp_inv <- try(solve(MHelp))
      if (!is.matrix(MHelp_inv)) stop("LM Test: Can't invert the information matrix MHelp_inv")
      MHelp2 <- IM_cor-cbind(t(IM_tv_cor),t(IM_garch_cor))%*%rbind(cbind(IM_tv_inv+IM_tv_inv%*%IM_tv_garch%*%MHelp_inv%*%t(IM_tv_garch)%*%IM_tv_inv,-IM_tv_inv%*%IM_tv_garch%*%MHelp_inv),cbind(-MHelp_inv%*%t(IM_tv_garch)%*%IM_tv_inv,MHelp_inv))%*%rbind(IM_tv_cor,IM_garch_cor)
      IM_inv <- solve(MHelp2)
    }
    if (is.null(nGARCH) && !is.null(nTV)){
      IM_tv_inv <- try(solve(IM_tv))
      if (!is.matrix(IM_tv_inv)) stop("LM Test: Can't invert the information matrix IM_tv")
      IM_inv <- solve(IM_cor-t(IM_tv_cor)%*%IM_tv_inv%*%IM_tv_cor)
    }
    if (!is.null(nGARCH) && is.null(nTV)){
      IM_garch_inv <- try(solve(IM_garch))
      if (!is.matrix(IM_garch_inv)) stop("LM Test: Can't invert the information matrix IM_garch")
      IM_inv <- solve(IM_cor-t(IM_garch_cor)%*%IM_garch_inv%*%IM_garch_cor)
    }
    if (is.null(nGARCH) && is.null(nTV)){
      IM_inv <- solve(IM_cor)
    }
    if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")
    # block corresponding to the corr.parameters that are set to zero under null
    SE_dim <- testOrder*(N*(N-1))/2
    block_start <- NCOL(IM_inv)-SE_dim+1
    block_end <- NCOL(IM_inv)
    IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]
    # LM statistic
    LM <- (1/Tobs)*t(dlldrho_A)%*%IM_inv_SE%*%dlldrho_A
    # Return:
    LM
  }
  

}  # End: myTest.CCCvSTCC.LM.new2()

#### ================== Test Functions ================== ###



#### ==================  Log-Liklihood Functions  ================== ####
myLogLik.tv_univar <- function(optimpars,e,tv,return_ll=TRUE) {
  # model: TV (no GARCH), univariate
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars   -- parameter vector including: (d0,delta(i),speed(i),loc1(i),<loc2(i)>) for all G(i)
  #        e           -- returns (not standardised)
  #        tv          -- TV object 
  # output:       return_ll   -- return value indicator:
  #               TRUE  : sum of loglikelihood(t) values over t=1...T
  #               FALSE : conditional variances, T rows
  # Note:  tv$var_target: TRUE => delta0 is free => include delta0 in optimpars
  
  # Create a new TV list within the scope of this function.
  # Extract the parameters from optimpars back into the tv() object - then check constraints:
  TV <- getEstimatedTVpars(tv,optimpars)

  #### ======== constraint checks ======== ####
  
  # Check 1. Check that delta0 is positive
  if (TV$Estimated$delta0 < 0) return(err_output) 
  
  # Checks 2 - 6 (if we have any more TV parameters, Note: no constraints on $linpars )
  if (TV$shape[1] > TRshape$none) {
    
    st <- TV$st
    shape <- TV$shape
    speedoption <- TV$speedoption
    numCol <- NCOL(TV$Estimated$pars)
    vecSpeed <- .subset(TV$Estimated$pars,"speed",1:numCol) 
    vecLoc1 <- .subset(TV$Estimated$pars,"loc1",1:numCol)
    vecLoc2 <- .subset(TV$Estimated$pars,"loc2",1:numCol)
    
    # Check 2: Check the boundary values for speed params:
    #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
    maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
    if (max(vecSpeed) > maxSpeed) return(err_output)
    if (min(vecSpeed) < 0) return(err_output)
    
    # Check 3: Check the loc1 locations fall within min-max values of st
    # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
    if (min(vecLoc1) < min(st)) return(err_output)
    if (max(vecLoc1) > max(st)) return(err_output)
    
    # Check 4: Check that loc1.1 < loc1.2 .. locN.1 < locN.2 for all G(i)
    # Method: Subtract loc1_pos vector from loc2_pos vector and ensure it is positive:
    tmp <- vecLoc2 - vecLoc1
    # Note: tmp will contain NA wherever a loc2 element was NA - we can strip these out:
    if (sum(tmp < 0,na.rm = TRUE) > 0) return(err_output)
    
    # Check 5: Check the loc2 locations fall within min-max values of st
    # Confirm we have at least one valid numeric loc 2, before checking min & max:
    if (any(!is.na(vecLoc2))) {
      if (min(vecLoc2,na.rm = TRUE) < min(st)) return(err_output)
      if (max(vecLoc2,na.rm = TRUE) > max(st)) return(err_output)  
    }
    
    # Check 6: Check that loc1.1 < loc2.1 where 2 locations exist... for all G(i)
    # We do need to have at least 2 locations for this error check
    if (NROW(vecLoc1) > 1) { 
      v1 <- head(vecLoc1,-1)
      v2 <- tail(vecLoc1,-1)
      if (sum(v2-v1 < 0) > 0) return(err_output)
    }
    
  }  # End: If (!is.na(TV$shape[1])) )
  
  #### ======== calculate loglikelihood ======== ####
  
  # recursively (over t=1...Tobs) compute g(t)=delta0 [+ delta1*G1(t) [+ delta2*G2(t) [+ ... ]]]
  g <- calculate_g(TV)
  
  #Check that g is positive! - No negative or NA elements allowed:
  if (min(g) < 0) return(err_output)
  
  #Return:
  if (return_ll) return( sum(-0.5*log(2*pi) - 0.5*log(g) - 0.5*(e*e)/g) ) else return (g) 
  
  
}  #End: myLogLik.tv_univar

myLogLik.garch_univar <- function(optimpars,w,garch,return_ll=TRUE){
  # model: GARCH
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars = parameter vector: (<omega>,<alpha>,beta,<delta>)
  #                    Note: The actual optimpars parameters passed-in depend on the Garch$type and Garch$var_target
  #                    garch$pars will always contain 4 elements as follows:
  #                       omega <- optimpars[1]
  #                       alpha <- optimpars[2]
  #                       beta <- optimpars[3]
  #                       delta <- optimpars[4]
  #        w = returns/sqrt(g) = e/sqrt(g)
  #        type = 1:GARCH or 2:GJR or 3:GJR with alpha=0
  
  
  # Check if variance targetting is being used. If not, set default to optimise ALL parameters.
  if (is.null(garch$var_target)) garch$var_target <- FALSE
  
  # Extract parameters from optimpars & do error-checking:
  garchPars <- constructGARCHpars(garch,optimpars)
  
  #### ======== constraint checks ======== ####
  if (garchPars["omega"] <= 0) return(err_output)
  if (garchPars["alpha"] < 0) return(err_output)
  if (garchPars["beta"] < 0) return(err_output)
  if (garchPars["gamma"] < 0) return(err_output)
  
  
  #### ======== calculate loglikelihood ======== ####
  Tobs <- length(w)
  h <- rep(0,Tobs)
  h[1] <- sum(w*w)/Tobs
  # min(w[t-1],0) => Only return negative values or 0
  for(t in 2:Tobs) h[t] <- garchPars["omega"] + garchPars["alpha"]*(w[t-1])^2 + garchPars["beta"]*h[t-1] + garchPars["gamma"]*(min(w[t-1],0))^2
  
  #Return:
  if (return_ll)  return( sum(-0.5*log(2*pi) -0.5*log(h) -0.5*(w*w)/h) ) else return (h)
  
}  #End: myLogLik.garch_univar()

myLogLik.stcc <- function(optimpars,z,stcc,return_ll=TRUE){
  # model: STCC
  # input: pars        -- c(speed,loc1,[loc2])
  #        z           -- volatility standardised returns (matrix TxN)
  #        stcc        -- list containing all the other parameters

  STCC <- stcc
  speedoption <- STCC$speedoption
  shape <- STCC$shape  
  st <- STCC$st
  
  Tobs <- NROW(z)
  N <- NCOL(z)

  #### ======== constraint checks ======== ####
  
  # Check 1: Confirm we have a valid shape & extract speed & location:
  if (shape==TRshape$none) return(err_output)
  if(shape==TRshape$double) numTRpars <- 3 else numTRpars <- 2
  TRpars <- tail(optimpars,numTRpars)
  speed <- TRpars[1]
  loc1 <- TRpars[2]
  if(numTRpars==3) loc2 <- TRpars[3] else loc2 <- NA
  
  # Check 2: Check the boundary values for speed params:
  #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
  if (speed > maxSpeed) return(err_output)
  if (speed < 0) return(err_output)
  
  # Check 3: Check the locations fall within min-max values of st
  # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
  if (loc1 < min(st)) return(err_output)
  if (loc1 > max(st)) return(err_output)
  if (numTRpars==3) {
    if (loc2 < min(st)) return(err_output)
    if (loc2 > max(st)) return(err_output)  
  }

  
  #### ======== calculate loglikelihood ======== ####
  # - - - CCC - - -
  if (STCC$type==STCCtype$CCC){
    vP <- c(0,0,0) #pars
    mP <- myUnVecl(vP)
    eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
    if (min(eig$values) <= 0) return(err_output)
    
    # - - - P(t) and loglik-value
    llt <- rep(0,Tobs)
    mPinv <- solve(mP)
    const <- -0.5*log(det(mP))
    for(t in seq(1,Tobs)) llt[t] <- const - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
  }
  
  # - - - STCC - - -
  if (STCC$type==STCCtype$STCC) {
    numCovPars <- NROW(myVecl(STCC$P1))
    vP1 <- optimpars[1:numCovPars]
    vP2 <- optimpars[(numCovPars+1):(2*numCovPars)]
    
    mP1 <- myUnVecl(vP1)
    eig1 <- eigen(mP1,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig1$values) <= 0) return(err_output)
    mP2 <- myUnVecl(vP2)
    eig2 <- eigen(mP2,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig2$values) <= 0) return(err_output)
    
    # - - - Calculate G(t) - - -
    st_c <- switch(shape,st-loc1,(st-loc1)*(st-loc2),(st-loc1)*(st-loc1))
    G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
    Gt <- matrix(1/G_inv,nrow = Tobs,ncol = 1)
    
    # - - - P(t) and loglik-value
    #Pt <- matrix(0,nrow=Tobs,ncol=(N*(N-1)/2))
    llt <- NULL
    #calcPt <- function(X,P1,P2) ((1-X)*P1 + X*P2)
    #Pt <- t(apply(Gt,MARGIN = 1,FUN = calcPt, P1=vP1, P2=vP2))
    Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
    for(t in seq(1,Tobs))
    {
      #Pt[t,] <- (1-Gt[t])*vP1 + Gt[t]*vP2  #This line has been replaced by the apply() function above
      mPt <- myUnVecl(Pt[t,])
      mPtinv <- chol2inv(chol(mPt))
      llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(mPtinv)%*%z[t,])
    } # End: for loop

  } # End: if(STCC$type==STCCtype$STCC)
  
  # Return: 
  if (return_ll) return(sum(llt)) else {
      if (STCC$type==STCCtype$CCC) Pt <- matrix(vP,nrow=Tobs,ncol=length(vP),byrow=TRUE)
      return(Pt)
    }

}  #End: myLogLik.stcc()

myLogLik.multivar.TVGARCHSTCC <- function(optimpars,e,ntv,ngarch,stcc,focus,var_target=TRUE,return_ll=TRUE){
  ###  Inputs:  ###
  # optimpars
  # e
  # ntv
  # ngarch
  # stcc
  # focus
  # var_target      
  # return_ll

  Tobs <- NROW(e)   # time dimension
  N <- NCOL(e)      # number of variables
  llt_const <- -0.5*log(2*pi)*N
  llt <- vector("numeric",Tobs)
  # We need matrices of g(t)_n & h(t)_n to standardise the data & calculate the LL
  g_n <- h_n <- matrix(data=0,nrow = Tobs,ncol=N)
  parsLength <- vector("integer",N)
  z <- w <- e

  #### ======== Focus: TV ======== ####  
  if (focus=="TV") {

    for (n in 1:N) {
      # 1. We map the optimpars back into TV objects
      # Note: if var_target = TRUE, then optimpars has all the delta0's at the start, otherwise none
      ntv[[n]]$var_target <- var_target
      
      #parsLength[n] <- length(ntv[[n]]$Estimated$pars[!is.na(ntv[[n]]$Estimated$pars)])
      parsLength[n] <- length(ntv[[n]]$Estimated$parsVector) - 1
      
      if (var_target) ntv[[n]]$Estimated$delta0 <- optimpars[n] else ntv[[n]]$Estimated$delta0 <- ntv[[n]]$Estimated$lastdelta0
      if (var_target) endPos <- N + sum(parsLength) else endPos <- sum(parsLength)
      startPos <- endPos - parsLength[n] + 1
      ntv[[n]]$Estimated$pars <- constructTVpars(optimpars[startPos:endPos],ntv[[n]]$shape)

      TV <- ntv[[n]]
      
      ### ======== TV: constraint checks ======== ###
      # Check 1. Check that delta0 is positive
      if (TV$Estimated$delta0 < 0) return(err_output) 
      
      # Checks 2 - 6 (if we have any more TV parameters, Note: no constraints on $linpars )
      if (.subset(TV$shape,1) > TRshape$none) {
        
        st <- TV$st
        shape <- TV$shape
        speedoption <- TV$speedoption
        numCol <- NCOL(TV$Estimated$pars)
        vecSpeed <- .subset(TV$Estimated$pars,"speed",1:numCol)
        vecLoc1 <- .subset(TV$Estimated$pars,"loc1",1:numCol)
        vecLoc2 <- .subset(TV$Estimated$pars,"loc2",1:numCol)
        
        # Check 2: Check the boundary values for speed params:
        #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
        maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
        if (max(vecSpeed) > maxSpeed) return(err_output)
        if (min(vecSpeed) < 0) return(err_output)
        
        # Check 3: Check the loc1 locations fall within min-max values of st
        # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
        if (min(vecLoc1) < min(st)) return(err_output)
        if (max(vecLoc1) > max(st)) return(err_output)
        
        # Check 4: Check that loc1.1 < loc1.2 .. locN.1 < locN.2 for all G(i)
        # Method: Subtract loc1_pos vector from loc2_pos vector and ensure it is positive:
        tmp <- vecLoc2 - vecLoc1
        # Note: tmp will contain NA wherever a loc2 element was NA - we can strip these out:
        if (sum(tmp < 0,na.rm = TRUE) > 0) return(err_output)
        
        # Check 5: Check the loc2 locations fall within min-max values of st
        # Confirm we have at least one valid numeric loc 2, before checking min & max:
        if (any(!is.na(vecLoc2))) {
          if (min(vecLoc2,na.rm = TRUE) < min(st)) return(err_output)
          if (max(vecLoc2,na.rm = TRUE) > max(st)) return(err_output)  
        }
        
        # Check 6: Check that loc1.1 < loc2.1 where 2 locations exist... for all G(i)
        # We do need to have at least 2 locations for this error check
        if (NROW(vecLoc1) > 1) { 
          v1 <- head(vecLoc1,-1)
          v2 <- tail(vecLoc1,-1)
          if (sum(v2-v1 < 0) > 0) return(err_output)
        }
        
      }  # End: If (!is.na(TV$shape[1])) )
      
      
      # 2. We recalculate the 'g' for all n = 1..N  (g_n) using the optimpars
      g_n[,n] <- calculate_g(ntv[[n]])  
      # Update the nTV object, so we can return it!
      ntv[[n]]$Estimated$condvars <- g_n[,n]
    
      # 3. Now we get the 'h' for all n = 1..N  (h_n) from the passed-in nGarch object
      w[,n] <- e[,n]/sqrt(g_n[,n])
      
      # 3. Recalculate 'h'
      garchPars <- constructGARCHpars(ngarch[[n]],ngarch[[n]]$Estimated$pars)
      h_n[1,n] <- sum(w[,n]*w[,n])/Tobs
      for(t in 2:Tobs) h_n[t,n] <- garchPars["omega"] + garchPars["alpha"]*(w[t-1,n])^2 + garchPars["beta"]*h_n[t-1,n] + garchPars["gamma"]*(min(w[t-1,n],0))^2
      # Update the nGarch object, so we can return it!
      ngarch[[n]]$Estimated$condvars <- h_n[,n]

    } # End: for (n in 1:N)

  } # End: focus = TV

  #### ======== Focus: GARCH ======== ####  
  if (focus=="GARCH") {
    
    for (n in 1:N) {
      # 1. We map the optimpars back into a single column matrix, 'garchPars'
      
      if (ngarch[[n]]$type == GARCHtype$GJR) parsLength[n] <- 4 else parsLength[n] <- 3
      endPos <- sum(parsLength)
      startPos <- endPos - parsLength[n] + 1
      # Note: The constructGARCHpars() function will use the $var_target field on the garch[[n]] objects
      ngarch[[n]]$var_target <- var_target
      garchPars <- constructGARCHpars(ngarch[[n]],optimpars[startPos:endPos])
      
      ### ======== GARCH: constraint checks ======== ###
      if (garchPars["omega"] <= 0) return(err_output)
      if (garchPars["alpha"] < 0) return(err_output)
      if (garchPars["beta"] < 0) return(err_output)
      #if (garchPars["gamma"] < 0) return(err_output)
      ## =============   End of Error Checks   ============= ##
      
      
      # 2. Use the passed in 'g' data to generate w:
      g_n[,n] <- ntv[[n]]$Estimated$condvars
      w[,n] <- e[,n]/sqrt(g_n[,n])
      
      # 3. Recalculate 'h'
      h_n[1,n] <- sum(w[,n]*w[,n])/Tobs
      for(t in 2:Tobs) h_n[t,n] <- garchPars["omega"] + garchPars["alpha"]*(w[t-1,n])^2 + garchPars["beta"]*h_n[t-1,n] + garchPars["gamma"]*(min(w[t-1,n],0))^2
      # Update the nGarch object, so we can return it!
      ngarch[[n]]$Estimated$condvars <- h_n[,n]
      ngarch[[n]]$Estimated$pars <- garchPars

    } # End: for (n in 1:N)
    
  } # End: focus = GARCH
  
  #### ======== Focus: STCC ======== ####  
  if (focus=="STCC") {
    
    if(any(is.na(stcc$P1))) return( err_output)
    
    # 1. We map the optimpars back into STCC Pars
    speedoption <- stcc$speedoption
    shape <- stcc$shape  
    st <- stcc$st
    
    ### ======== STCC: constraint checks ======== ###
    
    if(shape==TRshape$double) numTRpars <- 3 else numTRpars <- 2
    TRpars <- tail(optimpars,numTRpars)
    speed <- TRpars[1]
    loc1 <- TRpars[2]
    if(numTRpars==3) loc2 <- TRpars[3] else loc2 <- NA
    
    # Check 2: Check the boundary values for speed params:
    #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
    maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
    if (speed > maxSpeed) return(err_output)
    if (speed < 0) return(err_output)
    
    # Check 3: Check the locations fall within min-max values of st
    # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
    if (loc1 < min(st)) return(err_output)
    if (loc1 > max(st)) return(err_output)
    if (numTRpars==3) {
      if (loc2 < min(st)) return(err_output)
      if (loc2 > max(st)) return(err_output)  
    }
    
    ###  ALL constraint checks PASSED  ###
    
    # - - - CCC - - -
    if (stcc$type==STCCtype$CCC){
      vP <- c(0,0,0) #pars
      mP <- myUnVecl(vP)
      eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
      if (min(eig$values) <= 0) return(err_output)
      
      # - - - P(t) and loglik-value
      llt <- rep(0,Tobs)
      mPinv <- solve(mP)
      const <- -0.5*log(det(mP))
      for(t in seq(1,Tobs)) llt[t] <- const - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
    }
    
    # - - - STCC - - -
    if (stcc$type==STCCtype$STCC) {
      numCovPars <- NROW(myVecl(stcc$P1))
      vP1 <- optimpars[1:numCovPars]
      vP2 <- optimpars[(numCovPars+1):(2*numCovPars)]
      
      mP1 <- myUnVecl(vP1)
      eig1 <- eigen(mP1,symmetric=TRUE,only.values=TRUE)
      # Check for SPD - positive-definite check:
      if (min(eig1$values) <= 0) return(err_output)
      mP2 <- myUnVecl(vP2)
      eig2 <- eigen(mP2,symmetric=TRUE,only.values=TRUE)
      # Check for SPD - positive-definite check:
      if (min(eig2$values) <= 0) return(err_output)
      
      # - - - Calculate G(t) - - -
      st_c <- switch(shape,st-loc1,(st-loc1)*(st-loc2),(st-loc1)*(st-loc1))
      G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
      Gt <- matrix(1/G_inv,nrow = Tobs,ncol = 1)
      
      #calcPt <- function(X,P1,P2) ((1-X)*P1 + X*P2)
      #Pt <- t(apply(Gt,MARGIN = 1,FUN = calcPt, P1=vP1, P2=vP2))
      stcc$Estimated$condcorrs <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))

      # Finally, get the h_n & g_n from the passed in objects:
      for (n in 1:N) {
        h_n[,n] <- ngarch[[n]]$Estimated$condvars 
        g_n[,n] <- ntv[[n]]$Estimated$condvars
      }
      
    }
    
  } # End: focus = STCC

 
  # Finally, calculate ll(t) using g_n, h_n & STCC$rho
  ### ======== calculate loglikelihood ======== ###
  
  z <- e/(sqrt(h_n)*sqrt(g_n))
  P <- stcc$Estimated$condcorrs
  llt_const <- -0.5*log(2*pi)*N
  llt <- vector("numeric",Tobs)
  
  for(t in seq_along(1:Tobs)) {
    mPt <- myUnVecl(P[t,])
    mPtinv <- solve(mPt)
    llt[t] <- llt_const -0.5*sum(log(h_n[t,])) -0.5*sum(log(g_n[t,])) -0.5*log(det(mPt)) -0.5*(t(z[t,])%*%(mPtinv)%*%z[t,])
  }

  # Return:
  if (return_ll) return(sum(llt)) else {
    RTN <- list()
    RTN$value <- sum(llt)
    RTN$ntv <- ntv
    RTN$ngarch <- ngarch
    try(stcc$Estimated$parsVector <- c(myVecl(stcc$Estimated$P1),myVecl(stcc$Estimated$P2),stcc$Estimated$TRpars))
    RTN$stcc <- stcc
    #
    RTN
  }
  
}  #End: myLogLik.multivar.TVGARCHSTCC()

#### ==================  Log-Liklihood Functions  ================== ###


#### ==============  Estimation Functions  ============== ####

#===============================================================================#
#      Wrapper Functions below - These call the lower level functions above.
#   - Try to use these functions instead of calling the ones above directly -
#===============================================================================#

EstimateTV <- function(e,tv,calcHess=FALSE,verbose=FALSE) {
  
  TV <- tv
  
  ###
  ### ---   Set up 'optimpars' and call optim()  --- ###
  ###
  
  # First we check for the simple case of just delta0 provided, no TV$pars or TV$linpars:
  if (TV$shape[1]==TRshape$none && NROW(TV$linpars) == 1) {
    #TV Ordser 0 function - delta0 is the only param, so skip the optim()
    TV$Estimated$delta0 <- var(e) #Quick estimate for delta0 only
    TV$Estimated$lastdelta0 <- TV$Estimated$delta0
    TV$Estimated$condvars <- rep(TV$Estimated$delta0,NROW(e))
    TV$Estimated$value <- myLogLik.tv_univar(TV$Estimated$delta0,e,TV,return_ll=TRUE)
    TV$Estimated$error <- FALSE
    return(TV)
    
  } else {
    #else call optim() to get the estimate
    if (verbose) TV$optimcontrol$trace <- 10
    
    optimpars <- NULL
    # if Var_Target is ON (TRUE), then delta0 is a free param and will be passed to optim,
    # otherwise, we will not pass it - the log.lik() will use the last estimated value
    if (TV$var_target) optimpars <- c(optimpars,TV$delta0)

    if(TV$shape[1] > TRshape$none) optimpars <- c(optimpars,getTVparsVector(TV$pars))
    if(NROW(TV$linpars) > 1) optimpars <- c(optimpars,TV$linpars[!is.na(TV$linpars)])  

    # Now call optim:
    tmp <- NULL
    try(tmp <- optim(optimpars,myLogLik.tv_univar,e,TV,gr=NULL,method="BFGS",control=TV$optimcontrol,hessian=calcHess))      
  }

  ###    
  ### ---  Interpret response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    TV$Estimated$value <- err_output
    TV$Estimated$error <- TRUE
    return(TV)
  }
  
  if (tmp$convergence==0) { 
    #Optim converged successfully => we expect tmp$par to have good estimates!
    TV$Estimated$value <- tmp$value
    TV$Estimated$error <- FALSE
    # Now get the conditional variances
    TV$Estimated$condvars <- myLogLik.tv_univar(tmp$par,e,TV,return_ll=FALSE)
    #Update the TV object paramters using optimised pars:
    TV <- getEstimatedTVpars(TV,tmp$par)
    TV$Estimated$lastdelta0 <- tmp$par[1]
    TV$Estimated$parsVector <- getTVparsVector(TV$Estimated$pars)
    TV$Estimated$parsVector <- c(TV$Estimated$delta0,TV$Estimated$parsVector)

    if (calcHess) TV$Estimated$hessian <- tmp$hessian
    } 
  else { 
      #Failed to converge 
      TV$Estimated$value <- err_output
      TV$Estimated$error <- TRUE
      TV$Estimated$optimoutput <- tmp
      warning("EstimateTV() - failed to converge. Check the optim controls & starting params ")
  }

  # Return all the optimiser details
  if(verbose) TV$Estimated$optimoutput <- tmp
  rm(tmp)
  
  #Return:
  TV
}  #End: EstimateTV()

calcStderr_TV <- function(e,tv) {
  library(numDeriv)
  numDeriv::hessian()
  
  TV <- tv
  TV$Estimated$stderr <- NULL
  
  if(is.null(TV$Estimated$hessian)) {
    warning("This method of generating Standard Errors is unreliable. \nPlease re-estimate using the 'calcHess=TRUE' parameter")
    optimpars <- NULL
    if (TV$var_target) optimpars <- c(optimpars,TV$Estimated$delta0)
    if(TV$shape[1] > TRshape$none) optimpars <- c(optimpars,getTVparsVector(TV$Estimated$pars))
    TV$Estimated$hessian <- optimHess(optimpars,myLogLik.tv_univar,gr=NULL,e,TV,control=TV$optimcontrol)
  }

  # Using LU decomp:
  TV$Estimated$stderr <- sqrt(-diag(solve(TV$Estimated$hessian)))
  # Using Choleski
  # try(TV$Estimated$stderr <- sqrt(-diag(chol2inv(chol(hess)))))  
  if (is.null(TV$Estimated$stderr)) {
    msgWarning <- "Failed to calculate Standard Error"
    warning(msgWarning)
  } else TV$Estimated$stderr <- round(TV$Estimated$stderr,6)
  
  #Return
  TV
}
calcParamStats_TV <- function(TV,sig_level_percent=5) {
  
  if(sig_level_percent < 1) {
    warning("Please enter an integer % value, e.g. 5 for 5%")
    return()
    }
  
  TV$Estimated$tStats <- vector(mode="numeric",length=length(TV$Estimated$parsVector))
  TV$Estimated$PValues <- vector(mode="numeric",length=length(TV$Estimated$parsVector))
  TV$Estimated$fStat <- NULL
  
  pars <- TV$Estimated$parsVector
  errs <- TV$Estimated$stderr
  
  sig_level_percent = 5
  
  sig_value = (100 - sig_level_percent)/100 * 2
  degFreedom <- TV$Tobs - length(pars) - 1
  
  # Calculate t-stat = test if param is significantly different from zero
  TV$Estimated$tStats <- round(pars/errs,6)
  
  # Calculate P-Values based on the above t-Stats, for the significance level provided
  FUN1 <- function(x) {round(sig_value*pt(-abs(x),df=degFreedom),6)}
  TV$Estimated$PValues <- vapply(TV$Estimated$tStats, FUN1, vector("numeric",length = 1))
  
  # Add an asterix to each P-Value (in the name), to help users understand the parameter significance
  FUN2 <- function(x) {ifelse((abs(x) < sig_level_percent/100), "*", "") }
  names(TV$Estimated$PValues) <- paste0(names(TV$Estimated$PValues),vapply(TV$Estimated$PValues, FUN2, vector("character",length = 1)))
  
  
  # ## Loop logic for the vapply functions above
  # for(i in 1:(length(pars))) {
  # 
  #   TV$Estimated$tStat[i] <- pars[i]/errs[i]
  #   TV$Estimated$PValues[i] <- sig_value*pt(-abs(TV$Estimated$tStat[i]),df=degFreedom)
  #   # Indicate point-estimate significance 
  #   if((abs(TV$Estimated$PValues[i])) < sig_level_percent/100 ) names(TV$Estimated$pars[i]) = "*" else names(TV$Estimated$pars[i]) = " "
  #   
  # }
  
  return(TV)
  
}

calcModelStats_TV <- function(TV,sig_level_percent=100) {
  
  if(sig_level_percent < 1) {
    warning("Please enter an integer % value, e.g. 5 for 5%")
    return()
  }
  
  TV$Estimated$R2 <- NULL
  TV$Estimated$Residuals <- vector(mode="numeric",length=TV$Tobs)
  TV$Estimated$fStat <- NULL
  
  #
  #
  
  return(TV)
  
}

EstimateGARCH <- function(w,garch,calcHess=FALSE,verbose=FALSE) {
  
  ###
  ### ---   Set up 'optimpars' and call optim()  --- ###
  ###
  
  GARCH <- garch

  # First we check for the simple case of Garch$type = constant
  if (GARCH$type == GARCHtype$constant) {
    #No Garch - just constant unconditional variance, so skip the optim()
    GARCH$Estimated$pars[1] <- var(w) 
    GARCH$Estimated$condvars <- rep(GARCH$Estimated$pars[1],NROW(w))
    GARCH$Estimated$value <- myLogLik.garch_univar(GARCH$pars,w,GARCH,return_ll=TRUE)
    GARCH$Estimated$error <- FALSE
    return(GARCH)
  } else {
    #else call optim() to get the estimate
    if (verbose) GARCH$optimcontrol$trace <- 1
    # if Var_Target is ON (TRUE), then omega (GARCH$pars[1]) is NOT a free param and will be calculated,
    # otherwise, we will pass it to optim(), which will find the best value
    if (GARCH$var_target) optimpars <- tail(GARCH$pars,-1) else optimpars <- GARCH$pars
    tmp <- NULL
    try (tmp <- optim(optimpars,myLogLik.garch_univar,w,GARCH,gr=NULL,method="BFGS",control=GARCH$optimcontrol,hessian=calcHess))
  }

  ###    
  ### ---  Interpret response from optim --- ###
  ###

  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    GARCH$Estimated$value <- err_output
    GARCH$Estimated$error <- TRUE
    return(GARCH)
  }
  
  if (tmp$convergence==0) 
    {
      #Optim converged successfully => we expect tmp$par to have good estimates!
      GARCH$Estimated$pars <- constructGARCHpars(GARCH,tmp$par)   
      GARCH$Estimated$parsVector <- as.vector(tmp$par)
      names(GARCH$Estimated$parsVector) <- c("omega","alpha","beta","gamma")
      GARCH$Estimated$condvars <- myLogLik.garch_univar(tmp$par,w,GARCH,return_ll=FALSE)
      GARCH$Estimated$value <- tmp$value
      GARCH$Estimated$error <- FALSE
      if (calcHess) GARCH$Estimated$hessian <-tmp$hessian
    } else 
      {
        #Failed to converge 
        GARCH$Estimated$value <- err_output
        GARCH$Estimated$error <- TRUE
        GARCH$Estimated$optimoutput <- tmp
        warning("EstimateGARCH() - failed to converge. Check the optim controls & starting params ")
        return(GARCH)
      }
  
  # Return all the optimiser details
  if(verbose) GARCH$Estimated$optimoutput <- tmp
  rm(tmp)
  
  #Return:
  GARCH
}  #End: EstimateGARCH()

calcStderr_GARCH <- function(w,garch) {
  GARCH <- garch
  GARCH$Estimated$stderr <- NULL
  
  if(is.null(GARCH$Estimated$hessian)) {
    warning("This method of generating Standard Errors is unreliable. \nPlease re-estimate using the 'calcHess=TRUE' parameter")
    optimpars <- NULL
    # if Var_Target is ON (TRUE), then omega (GARCH$pars[1]) is NOT a free param and will be calculated,
    # otherwise, we will pass it to optim(), which will find the best value
    if (GARCH$var_target) optimpars <- tail(GARCH$Estimated$pars,-1) else optimpars <- GARCH$Estimated$pars
    GARCH$Estimated$hessian <- optimHess(optimpars,myLogLik.garch_univar,gr=NULL,w,GARCH,control=GARCH$optimcontrol)
  }

  # Using LU decomp:
  #Use try(..) here as solve can fail!
  GARCH$Estimated$stderr <- sqrt(-diag(solve(GARCH$Estimated$hessian)))
  # Using Choleski
  # try(TV$Estimated$stderr <- sqrt(-diag(chol2inv(chol(hess)))))  
  if (is.null(GARCH$Estimated$stderr)) {
    msgWarning <- "Failed to calculate Standard Error"
    warning(msgWarning)
  } else GARCH$Estimated$stderr <- round(GARCH$Estimated$stderr,6)
  
  #Return
  GARCH
  
}

calcParamStats_GARCH <- function(GARCH,sig_level_percent=5) {

  if(sig_level_percent < 1) {
    warning("Please enter an integer % value, e.g. 5 for 5%")
    return()
  }
  
  GARCH$Estimated$tStats <- vector(mode="numeric",length=length(GARCH$Estimated$pars))
  GARCH$Estimated$PValues <- vector(mode="numeric",length=length(GARCH$Estimated$pars))
  GARCH$Estimated$fStat <- NULL
  
  pars <- as.vector(GARCH$Estimated$pars)
  errs <- GARCH$Estimated$stderr
  
  sig_level_percent = 5
  
  sig_value = (100 - sig_level_percent)/100 * 2

  # Calculate t-stat = test if param is significantly different from zero
  GARCH$Estimated$tStats <- round(pars/errs,6)
  
  # Calculate P-Values based on the above t-Stats, for the significance level provided
  FUN1 <- function(x) {round(sig_value*pnorm(-abs(x)),6)}
  GARCH$Estimated$PValues <- vapply(GARCH$Estimated$tStats, FUN1, vector("numeric",length = 1))
  
  # Add an asterix to each P-Value (in the name), to help users understand the parameter significance
  FUN2 <- function(x) {ifelse((abs(x) < sig_level_percent/100), "*", "") }
  names(GARCH$Estimated$PValues) <- paste0(names(GARCH$Estimated$PValues),vapply(GARCH$Estimated$PValues, FUN2, vector("character",length = 1)))
  
  #Return
  GARCH
  
}

EstimateSTCC <- function(z,stcc,calcHess=FALSE,verbose=FALSE) {
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###

  STCC <- stcc
  numCovPars <- NROW(myVecl(STCC$P1))
  optimpars <- c(myVecl(STCC$P1),myVecl(STCC$P2),STCC$TRpars)
  if (verbose) STCC$optimcontrol$trace <- 10
  tmp <- NULL
  try(tmp <- optim(optimpars,myLogLik.stcc,z,STCC,gr=NULL,method="BFGS",control=STCC$optimcontrol,hessian=calcHess))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STCC$Estimated$value <- err_output 
    STCC$Estimated$error <- TRUE
    return(STCC)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    if(STCC$shape==TRshape$double) numTRpars <- 3 else numTRpars <- 2
    if (!is.null(STCC$P1)) STCC$Estimated$P1 <- myUnVecl(tmp$par[1:numCovPars])
    if (!is.null(STCC$P2)) STCC$Estimated$P2 <- myUnVecl(tmp$par[(numCovPars+1):(2*numCovPars)])
    if (!is.null(STCC$TRpars)) STCC$Estimated$TRpars <- tail(tmp$par,numTRpars)
    # TODO: fix this parsVector code - remove try()
    try(STCC$Estimated$parsVector <- c(myVecl(STCC$Estimated$P1),myVecl(STCC$Estimated$P2),STCC$Estimated$TRpars))
    
    STCC$Estimated$value <- tmp$value
    STCC$Estimated$condcorrs <- myLogLik.stcc(tmp$par,z,STCC,return_ll=FALSE)
    STCC$Estimated$error <- FALSE
    if (calcHess) STCC$Estimated$hessian <- tmp$hessian

  } else { 
    #Failed to converge
    STCC$Estimated$error <- TRUE
    STCC$Estimated$value <- err_output
    STCC$Estimated$optimoutput <- tmp
  }

  if (verbose) STCC$Estimated$optimoutput <- tmp
  #Return:
  STCC

}  #End: EstimateSTCC()

calcStderr_STCC <- function(e,stcc) {
  STCC <- stcc
  STCC$Estimated$stderr <- NULL
  
  if(is.null(STCC$Estimated$hessian)) {
    warning("This method of generating Standard Errors is unreliable. \nPlease re-estimate using the 'calcHess=TRUE' parameter")
    optimpars <- NULL
    optimpars <- c(myVecl(STCC$Estimated$P1),myVecl(STCC$Estimated$P2),STCC$Estimated$TRpars)
    STCC$Estimated$hessian <- optimHess(optimpars,myLogLik.tv_univar,gr=NULL,e,STCC,control=STCC$optimcontrol)
  }
  
  # Using LU decomp:
  STCC$Estimated$stderr <- sqrt(-diag(solve(STCC$Estimated$hessian)))
  # Using Choleski
  # try(TV$Estimated$stderr <- sqrt(-diag(chol2inv(chol(hess)))))  
  if (is.null(STCC$Estimated$stderr)) {
    msgWarning <- "Failed to calculate Standard Error"
    warning(msgWarning)
  } else STCC$Estimated$stderr <- round(STCC$Estimated$stderr,6)
  
  #Return
  STCC
}

calcParamStats_STCC <- function(STCC,sig_level_percent=5) {
  
  if(sig_level_percent < 1) {
    warning("Please enter an integer % value, e.g. 5 for 5%")
    return()
  }
  
  pars <- STCC$Estimated$parsVector
  errs <- STCC$Estimated$stderr
  
  STCC$Estimated$tStats <- vector(mode="numeric",length=length(pars))
  STCC$Estimated$PValues <- vector(mode="numeric",length=length(pars))
  STCC$Estimated$fStat <- NULL

  sig_level_percent = 5
  sig_value = (100 - sig_level_percent)/100 * 2
  
  # Calculate t-stat = test if param is significantly different from zero
  STCC$Estimated$tStats <- round(pars/errs,6)
  
  # Calculate P-Values based on the above t-Stats, for the significance level provided
  FUN1 <- function(x) {round(sig_value*pnorm(-abs(x)),6)}
  STCC$Estimated$PValues <- vapply(STCC$Estimated$tStats, FUN1, vector("numeric",length = 1))
  
  # Add an asterix to each P-Value (in the name), to help users understand the parameter significance
  FUN2 <- function(x) {ifelse((abs(x) < sig_level_percent/100), "*", "") }
  names(STCC$Estimated$PValues) <- paste0(names(STCC$Estimated$PValues),vapply(STCC$Estimated$PValues, FUN2, vector("character",length = 1)))

  #Return
  STCC
  
}

Estimate_MTVGARCH <- function(e,ntv,ngarch,stcc,focus,var_target=TRUE){
  ###  Inputs:  ###
  # e
  # ntv
  # ngarch
  # stcc
  # focus
  # var_target      
  
  Tobs <- NROW(e)   # time dimension
  N <- NCOL(e)      # number of variables
  RTN <- list()     # list object to hold the return objects; ll_value,m_pars,ntv,ngarch,stcc
  RTN$value <- NA
  RTN$ntv <- ntv
  RTN$ngarch <- ngarch
  RTN$stcc <- stcc
  RTN$error <- FALSE
  RTN$optimoutput <- list()
  RTN$m_pars <- NA
  
  optimpars <- NULL
  
  # We need matrices of g(t)_n & h(t)_n to standardise the data & calculate the LL
  g_n <- h_n <- matrix(data = 0,nrow = Tobs,ncol = N)
  
  if (focus=="TV") {

    # # Put all the delta0's into optimpars first
    # if (var_target) for (n in 1:N) optimpars <- c(optimpars,ntv[[n]]$Estimated$delta0)
    # # Now add the other parameters
    # for (n in 1:N) optimpars <- c(optimpars,getTVparsVector(ntv[[n]]$Estimated$pars))  
    # # Set the optim controls
    # if (!is.null(optimcontrols)) TVoptimcontrol <- optimcontrols$TV else TVoptimcontrol <- list(fnscale = -1,maxit=1000,reltol = 1e-7)

    # Put all the delta0's into optimpars first
    parScale <- NULL
    if (var_target) for (n in 1:N) {
      optimpars <- c(optimpars,ntv[[n]]$delta0)
      parScale <- c(parScale,ntv[[n]]$optimcontrol$parscale[1]) 
    }
    # Now add the other parameters & build optimcontrol
    for (n in 1:N) {
      optimpars <- c(optimpars,getTVparsVector(ntv[[n]]$pars))  
      parScale <- c(parScale,tail(ntv[[n]]$optimcontrol$parscale,-1)) 
    }
    nDeps <- rep(1e-4,length(parScale))
    TVoptimcontrol <- list(fnscale = -1,ndeps=nDeps,parscale=parScale,maxit=2000,reltol=1e-7)
    
    ## Call optim
    tmp = optim(optimpars,myLogLik.multivar.TVGARCHSTCC,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE,
                    gr=NULL,method="BFGS",control=TVoptimcontrol)
    
    ## Interpret response from optim:
    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
      RTN$value <- err_output
      RTN$error <- TRUE
      return(RTN)
    }
    
    if (tmp$convergence==0) {
      RTN <- myLogLik.multivar.TVGARCHSTCC(tmp$par,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)
      RTN$optimoutput <- tmp
    } else stop("Error: optim failed to converge!")
    
    
    } # End: focus = TV
  
  if (focus=="GARCH") {
    
    # # Add the parameters to optimpars
    # for (n in 1:N) optimpars <- c(optimpars,ngarch[[n]]$Estimated$pars)
    # # Set the optim controls
    # if (!is.null(optimcontrols)) GARCHoptimcontrol <- optimcontrols$GARCH else GARCHoptimcontrol <- list(fnscale = -1,maxit=1000,reltol = 1e-7)
    
    # Add the parameters to optimpars
    parScale <- NULL
    nDeps <- NULL
    for (n in 1:N) {
      parVec <- as.vector(ngarch[[n]]$pars)
      optimpars <- c(optimpars,parVec)
      parScale <- c(parScale,ngarch[[n]]$optimcontrol$parscale)
      nDeps <- c(nDeps,ngarch[[n]]$optimcontrol$ndeps)
    }
    GARCHoptimcontrol <- list(fnscale = -1,ndeps=nDeps,parscale=parScale,maxit=2000,reltol=1e-8)
    
    ## Call optim
    tmp = optim(optimpars,myLogLik.multivar.TVGARCHSTCC,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE,
                    gr=NULL,method="BFGS",control=GARCHoptimcontrol)
    
    ## Interpret response from optim:
    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
      RTN$value <- err_output
      RTN$error <- TRUE
      return(RTN)
    }
    
    if (tmp$convergence==0) {
      RTN <- myLogLik.multivar.TVGARCHSTCC(tmp$par,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)
      RTN$optimoutput <- tmp
    } else stop("Error: optim failed to converge!")
    
  } # End: focus = GARCH
  
  if (focus=="STCC") {
    
    # Add the parameters to optimpars
    #optimpars <- c(optimpars,myVecl(stcc$Estimated$P1),myVecl(stcc$Estimated$P2),stcc$Estimated$TRpars)
    optimpars <- c(optimpars,myVecl(stcc$P1),myVecl(stcc$P2),stcc$TRpars)
    # # Set the optim controls
    # if (!is.null(optimcontrols)) STCCoptimcontrol <- optimcontrols$STCC else STCCoptimcontrol <- list(fnscale = -1,maxit=1000,reltol = 1e-7)
    
    parScale <- rep(1,length(optimpars))
    nDeps <- c(rep(1e-4,(length(optimpars) - 2)), 1e-5, 1e-4)
    STCCoptimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
    
    ## Call optim
    tmp = optim(optimpars,myLogLik.multivar.TVGARCHSTCC,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE,
                    gr=NULL,method="BFGS",control=STCCoptimcontrol)
    
    ## Interpret response from optim:
    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
      RTN$value <- err_output
      RTN$error <- TRUE
      return(RTN)
    }
    
    if (tmp$convergence==0) {
      RTN <- myLogLik.multivar.TVGARCHSTCC(tmp$par,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)
      RTN$optimoutput <- tmp
      
    } else stop("Error: optim failed to converge!")
    
  } # End: focus = STCC
  
  
  # Build Return object:
  m_tvpars <- m_garchpars <- NULL
  for (n in 1:N) {
    m_tvpars <- c(m_tvpars,ntv[[n]]$Estimated$delta0,as.vector(ntv[[n]]$Estimated$pars))
    m_garchpars <- c(m_garchpars,ngarch[[n]]$Estimated$pars)
  }
  m_tvpars <- m_tvpars[!is.na(m_tvpars)]
  RTN$m_pars <- c(m_tvpars,m_garchpars,myVecl(stcc$Estimated$P1),myVecl(stcc$Estimated$P2),stcc$Estimated$TRpars)
  
  # Return:
  RTN 
  
  
}  #End: Estimate_MTVGARCH()



## TODO: Fix the two following functions!
EstimateARCH <- function(w,arch,calcHess=FALSE) {
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  ARCH <- arch
  #if (is.null(GARCH$var_target)) GARCH$var_target <- FALSE  #Default value
  if (is.null(ARCH$optimcontrol)) stop("ARCH$optimcontrol is mandatory.")
  
  ## First, set the optimpars, based on Garch$type & Garch$var_target:
  optimpars <- ARCH$pars
  ## Second, set the optim Controls, based on Garch$var_target:
  myControl <- ARCH$optimcontrol
  tmp <- NULL
  try (tmp <- optim(optimpars,myLogLik.arch_univar,w,ARCH,gr=NULL,method="BFGS",control=myControl,hessian=calcHess))
  
  ###    
  ### ---  Interpret the response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    ARCH$value <- -111333 #Set the value = the ErrorOutput value used by the myLogLik.garch_univar() function
    ARCH$error <- TRUE
    return(ARCH)
  }
  
  if (tmp$convergence==0) {
    #Optim converged successfully => we expect tmp$par to have good estimates!
    ARCH$pars <- tmp$par
    ARCH$condvars <- myLogLik.arch_univar(tmp$par,w,ARCH,return_ll=FALSE)
    ARCH$value <- tmp$value
    ARCH$error <- FALSE
    if (calcHess) {
      ARCH$hess <-tmp$hessian
      ARCH$stderr <- NULL
      try(ARCH$stderr <- sqrt(-diag(solve(tmp$hessian))))
      if (is.null(ARCH$stderr)) warning("Failed to calculate ARCH$stderr - probably a numerical solve() problem")
    }
    
  } else {
    #Failed to converge - return specific error code in the $value
    ARCH$value <- -111000-tmp$convergence
    ARCH$error <- TRUE
    return(ARCH)
  }
  
  rm(tmp)
  
  #Return:
  ARCH
}  #End: EstimateARCH()

EstimateTVGARCH <- function(e,series,maxLoops=20,verbose=FALSE) {
  ###
  ### ---  Call optim for TV & GARCH in a loop until Converged --- ###
  ###
  
  TV <- series$tv
  GARCH <- series$garch
  
  #Set the variance targetting:
  TV$var_target <- series$var_target
  GARCH$var_target <- series$var_target
  
  # Initialise Repeat..loop control values:
  loopcount <- 0
  ll_const <- -0.5*log(2*pi)
  e2 <- e*e
  TV_Converged <- FALSE
  GARCH_Converged <- FALSE
  
  repeat {
    
    # parameter values before optimization
    if (!series$var_target) G_prevPars <- GARCH$pars else G_prevPars <- GARCH$pars[-1]
    G_prevPars <- G_prevPars[!is.na(G_prevPars)]
    G_prevData <- e/sqrt(TV$condvars)
    
    GARCH <- EstimateGARCH(G_prevData,GARCH)
    
    if (!series$var_target) currentPars <- GARCH$pars else currentPars <- GARCH$pars[-1]
    currentPars <- currentPars[!is.na(currentPars)]
    # Check if the change in param value is small enough to exit:
    if (series$var_target) {
      if (all((G_prevPars-currentPars) < series$garchparLimit[-1])) GARCH_Converged <- TRUE } 
    else {
      # var_target = FALSE  
      if (all(abs(G_prevPars-currentPars) < series$garchparLimit)) GARCH_Converged <- TRUE 
    }
    if (verbose) cat("\nGARCH value:",GARCH$value, "Converged=",GARCH_Converged )
    
    if (verbose) cat("\n garchPars-GARCH$pars:",(G_prevPars-currentPars))
    if (verbose) cat("\n garchparLimit:",series$garchparLimit)
    
    
    # parameter values before optimization
    if (!series$var_target) TV_prevPars <- TV$pars else TV_prevPars <- c(TV$delta0,TV$pars)
    TV_prevPars <- TV_prevPars[!is.na(TV_prevPars)]
    
    TV <- EstimateTV(e/sqrt(GARCH$condvars),TV)
    
    if (!series$var_target) currentPars <- TV$pars else currentPars <- c(TV$delta0,TV$pars)
    currentPars <- currentPars[!is.na(currentPars)]
    # Check if the change in param value is small enough to exit:
    if (!series$var_target) {
      if (all(abs(TV_prevPars-currentPars) < series$tvparLimit[-1])) TV_Converged <- TRUE }
    else {
      # Var_target = TRUE
      if (all(abs(TV_prevPars-currentPars) < series$tvparLimit)) TV_Converged <- TRUE
    }
    if (verbose) cat("\nTV value:",TV$value, "Converged=", TV_Converged)
    
    if (verbose) cat("\n tv$pars-TV$pars:",(TV_prevPars-currentPars))
    if (verbose) cat("\n tvparLimit:",series$tvparLimit)
    
    g <- TV$condvars
    h <- GARCH$condvars
    series$value <- sum(ll_const - 0.5*log(g) - 0.5*log(h) - 0.5*e2/(g*h))
    if (verbose) cat("\nSeries value:",series$value," in Loop:",loopcount)
    if (verbose) cat("\n")
    
    #Exit time:
    if ((loopcount==maxLoops) || (TV_Converged && GARCH_Converged)) {
      # Finally - Calculate the hessian matrix & std Errors
      
      GARCH$pars[!is.na(GARCH$pars)] <- G_prevPars
      GARCH <- EstimateGARCH(G_prevData,GARCH,calcHess = TRUE)
      
      if (!series$var_target) TV$pars[!is.na(TV$pars)] <- TV_prevPars else {
        TV$delta0 <- TV_prevPars[1]
        TV$pars[!is.na(TV$pars)] <- tail(TV_prevPars,-1)
      } 
      TV <- EstimateTV(e/sqrt(GARCH$condvars),TV,calcHess = TRUE)
      
      if (loopcount==maxLoops) {
        cat("\n\nMaximum Loops reached.  Convergence state is:")
        cat("\nTV_Converged:",TV_Converged,"\nGARCH_Converged:",GARCH_Converged)
      } else cat("\n\nTV & GARCH Convergence achieved in:",loopcount,"loops")
      
      break 
    }
    
    #Progress Bar:
    if (!verbose) cat(".")
    
    loopcount <- loopcount + 1
    
  } # End repeat..loop
  
  # Update the series with the optimised values
  series$tv <- TV
  series$garch <- GARCH
  
  #Return:
  series
  
} # End: EstimateTVGARCH()

#### ==================  Estimation Functions  ================== ###


#### ==================  Simulate Prob Distribution  ================== ####

GenTestStatDist <- function(e,tv,refdata,reftest,usetest=TESTtype$ALL,testorder=TESTorder$H0,saveas=NULL,numloops=1100) {
  
  # Make an explicit local copy of the TV object
  TVloc <- tv
  
  # 1. It is useful to know how long this process takes to run, so...
  tmr <- proc.time()
  timestamp(prefix = "Starting Simulation on ",suffix = "\nPlease be patient as this may take a while...\n")
  
  # 2. Create Output folder (if not there) & set Save filename
  if (!dir.exists(file.path(getwd(),"Output"))) dir.create(file.path(getwd(),"Output"))
  saveAs <- file.path("Output",saveas)
  
  # 3. Load the generated data with Garch and add the 'g' from our TV object
  Tobs <- NROW(e)
  RefData_WithGarch <- refdata[1:Tobs,1:numloops]
  
  ## Add 'g' into the Reference Data.  
  RefData_WithGarch <- RefData_WithGarch*sqrt(TVloc$Estimated$condvars)
  
  # 4. Select Tests to run
  if (usetest==TESTtype$TR2) {
    runTR2 <- TRUE
    runRobust <- FALSE
  } else if (usetest==TESTtype$Robust)  {
    runTR2 <- FALSE
    runRobust <- TRUE
  } else {
    runTR2 <- runRobust <- TRUE
  }
  
  # 5. Setup the matrix to store the simulation results - depends on the Order of TV function
  numTestResults <- as.integer(runTR2)*3 + as.integer(runRobust)*3   # ([TestStat,Pval,RefStat] x (1 or 2) Tests) 
  LogLiklihood_value <- 1
  if (NROW(getTVparsVector(TVloc$pars)) >  1) numpars <- 1 + NROW(getTVparsVector(TVloc$pars)) else numpars <- 1 
  if (NROW(TVloc$linpars) > 1) numpars <- numpars + NROW(TVloc$linpars[!is.na(TVloc$linpars)])
  numResultCols <- 1 + numTestResults + LogLiklihood_value + numpars  # The first column is used for the b-loop index 'b'
  testStatDist <- matrix(NA,nrow=numloops,ncol=numResultCols)
  
  # 6. Collect the required functions & variables to pass into each parallel environment
  export_functions <- c("calculate_g","dg.dt","dg.dt2.H0","myLogLik.tv_univar",
                        "getEstimatedTVpars","getTVparsVector","constructTVpars",
                        "TRshape","TRspeedopt","TESTorder","err_output", 
                        "EstimateTV","myTest.TV.noGARCH.TR2","myTest.TV.noGARCH.robust")
  
  # 7. Perform the simulation - in parallel
  testStatDist <- foreach(b = 1:numloops, .inorder=FALSE, .combine=rbind, .verbose = FALSE, .export = export_functions) %do% {
  
    runSimrow <- c(b,rep(NA,numTestResults))  
    sim_e <- as.vector(RefData_WithGarch[,b])
    
    
    TV <- EstimateTV(sim_e,TVloc)  # Note: The tv must always be the same params that are passed in!  Only the sim_e changes!!
    
    if (!TV$Estimated$error) {
      # Now run the requested Tests:
      if (runTR2) {
        simTEST1 <- myTest.TV.noGARCH.TR2(sim_e,TV,testorder)         
        runSimrow[2:4] <- c(simTEST1,as.integer(simTEST1 > reftest$LMTR2),reftest$LMTR2) 
      }
      if (runRobust) {
        simTEST2 <- myTest.TV.noGARCH.robust(sim_e,TV,testorder) 
        #TODO: There is a bug here with the indexing - only works when both tests are run - FIX
        runSimrow[5:7] <- c(simTEST2,as.integer(simTEST2 > reftest$LMRobust),reftest$LMRobust)
      }
    } else {
      # Optim failed, but we need these fields to fill in the runsimrow vector:
      TV$Estimated$delta0 <- TVloc$delta0
      TV$Estimated$pars <- TVloc$pars
      TV$Estimated$linpars <- TVloc$linpars[!is.na(TVloc$linpars)]
      TV$Estimated$value <- NaN
    } # End: if (!TV$Estimated$error) 

    result <- c(runSimrow,TV$Estimated$value,TV$Estimated$delta0) 
    if (TV$shape[1] > TRshape$none) result <- c(result,getTVparsVector(TV$Estimated$pars))
    if (NROW(TV$linpars) > 1) result <- c(result,TV$Estimated$linpars)
    #Return:
    result
    
  } # End: foreach(b = 1:numloops,...
  
  # 8. Add column names to results
  loopDataNames <- c("b")
  if (runTR2) loopDataNames <- c(loopDataNames,"Stat_TR2","Pval_TR2","Ref$LMTR2")
  if (runRobust) loopDataNames <- c(loopDataNames,"Stat_Robust","Pval_Robust","Ref$LMRobust")
  loopDataNames <- c(loopDataNames,"Estimated_LL")
  TVparNames <- c("TV$delta0")
  if(TV$shape[1] > TRshape$none) TVparNames <- c(TVparNames,paste0("TV$par",as.character(seq(1,length(getTVparsVector(TV$Estimated$pars))))))
  if(NROW(TV$linpars) > 1) TVparNames <- c(TVparNames,paste0("TV$linpar",as.character(seq(1,length(TV$Estimated$linpars)))))
  loopDataNames <- c(loopDataNames,TVparNames)
  try(colnames(testStatDist) <- loopDataNames)
  
  # 9. Save the distribution
  if (!is.null(saveas)) saveRDS(testStatDist,saveas)
  
  # 10. Extract Test P_Values from Results & express as %
  Test <- list()
  Test$p_TR2 <- 100*mean(testStatDist[,"Pval_TR2"],na.rm = TRUE)
  Test$p_ROB <- 100*mean(testStatDist[,"Pval_Robust"],na.rm = TRUE)
  Test$Ref_TR2 <- testStatDist[1,"Ref$LMTR2"]
  Test$Stat_TR2 <- testStatDist[,"Stat_TR2"]
  Test$Ref_Robust <- testStatDist[1,"Ref$LMRobust"]
  Test$Stat_Robust <- testStatDist[,"Stat_Robust"]
  Test$FailCount_TR2 <- length(Test$Stat_TR2[is.na(Test$Stat_TR2)])
  Test$FailCount_Robust <- length(Test$Stat_Robust[is.na(Test$Stat_Robust)])
  # Troubleshooting:
  #Test$dist <- testStatDist
  
  # 11. Print the time taken to the console:
  cat("Simulation Time:",(proc.time()-tmr)[3],"seconds\n")
  
  # 12. Attempt to release memory:
  rm(RefData_WithGarch,TVloc)
  
  # Return:
  Test
}

#### ==================  Simulate Prob Distribution  ================== ###


#### ==================  Data Generation  ================== ####

GenerateRefData <- function(corr,ngarch,ntv,e=NULL,e_discard=NULL,seedstart=0,numseries=10,tobs=1000,saveas="GenRefData.RDS")
{
  
## ----- Comments to explain how this function works ----- ##
##  This function is designed to generate a matrix of N data series, to be used as reference data with known qualities
##  Garch, TV, and Correlation information can be passed to the function as a parameter. None or all can be included.
##  If e is provided we will use it as the data - matrix of noise expected
##  If e_discard is provided we will use it as the garch discard data - matrix of noise expected
##  If e is not provided, then e & e_discard will be generated using the other parameters passed in.
  
# v   : this will store iid data
# z   : this will store correlated iid data 
# w   : this will store correlated data + GARCH 
# e   : this will store correlated data + GARCH + TV 

  
  Tobs <- tobs
  N <- numseries
  refData <- matrix(NA,Tobs,N)  # this will store the Tobs x N final generated data
  
  if (!is.null(e)) {
    v <- e
    v_discard <- e_discard
    discardObs <- NROW(v_discard)
  } else {
  
    # We need to have some data to discard, when generating Garch:
    discardObs <- 1500
    
    # Make some noise to start with!!
    set.seed(seedstart)
    dataVec <- rnorm((Tobs+discardObs)*N)
    v_vec <- dataVec[1:(Tobs*N)]
    # Use the tail of this vector for our discard:
    discard_vec <- tail(dataVec,(discardObs*N))
    
    v <- matrix(v_vec,nrow=Tobs,ncol=N)
    v_discard <- matrix(discard_vec,nrow=discardObs,ncol=N)
  }
  
  # Generate Correlated iid Data
  z <- v
  if (!is.null(corr)) {
    
      # - - - CCC - - -
      if (corr$type=="CCC"){
        CCC <- corr$CCC
        P <- CCC$P
        eig <- eigen(P)       # eigenvalues and eigenvectors
        Psqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
        z <- t(Psqrt%*%t(v))
      }
      
      # - - - STCC - - -
      else if (corr$type=="STCC"){
        STCC <- corr$STCC
        rhovec1 <- as.vector(STCC$P1)
        rhovec2 <- as.vector(STCC$P2[(N*(N-1)/2+1):(N*(N-1))])
        speed <- STCC$pars[1]
        loc1 <- STCC$pars[2]
        loc2 <- STCC$pars[3]
        
        Gt <- myG(speed,loc1,loc2,STCC$st,STCC$shape,STCC$speedoption)
        P <- matrix(0,Tobs,N)
        for (t in 1:Tobs){
          P[t,] <- (1-Gt[t])*rhovec1 + Gt[t]*rhovec2
          mP <- myUnVecl(P[t,])
          eig <- eigen(mP)           # eigenvalues and eigenvectors
          mPtsqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
          z[t,] <- mPtsqrt%*%v[t,]
        }
      }

  } #End: Generate Correlated Data
  
  
  # Now add Garch and/or TV to the data (loop over N-series)
  for (n in 1:N) {
    
    w <- z
    # Generate GARCH in the data:
    if (!is.null(ngarch)) {
      garch <- ngarch[[n]]
      # Set Garch params based on type:
      if (garch$type==0){
        par_o <- garch$pars[1]
        par_a <- 0
        par_b <- 0
        par_d <- 0
      }else if (garch$type==1){
        par_o <- garch$pars[1]
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- 0
      }else if (garch$type==2){
        par_o <- garch$pars[1]
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- garch$pars[4]
      }else if (garch$type==3){
        par_o <- calculateOmega(garch$pars[2],garch$pars[3],garch$pars[4])
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- garch$pars[4]
      } else stop("Garch$type is missing or invalid!")
      
      # Initialise the first data points: e(t-1), h(t-1) using the v_discard[<discard>,N]
      ht_1 <- ht <- 1
      et_1 <- et <- v_discard[1,n]
      for (t in 2:discardObs) {
        ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 + par_d*(min(et_1,0)^2)
        et_1 <- et <- sqrt(ht)*v_discard[t,n]
      }
      
      # Now generate the actual Garch data:  
      ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      w[1,n] <- sqrt(ht)*z[1,n]
      for (t in 2:Tobs) {
        ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 + par_d*(min(w[t-1,n],0)^2)
        w[t,n] <- sqrt(ht)*z[t,n]
      }
    } #End: Generate Garch Data
    
    e <- w
    # Generate TV in the data:
    if (!is.null(ntv)) {
      tv <- ntv[[n]]
      # Note: calculation of 'g' depends on tv$st to get the length of data series
      tv$Tobs <- NROW(e)
      tv$st <- seq(1,tv$Tobs)/tv$Tobs
      if (NROW(tv$pars) == 1) {
        # - - - NO TV (only constant delta0) - - -
        #g <- tv$delta0
        g <- rep(tv$delta0,Tobs)
        e[,n] <- w[,n]*sqrt(g)
      } else {
        # - - - Full TV - - -
        g <- calculate_g(tv)
        e[,n] <- w[,n]*sqrt(g)
      } 
    } #End: Generate TV Data 
    
    # Update the return matrix with this data series:
    refData[,n] <- e[,n]
    
    
  } #End: for (n in 1:N) 
  
  
  if(!is.null(saveas)) saveRDS(refData,saveas)
  #Return:
  refData
  
}  # End: GenerateRefData()


## TODO: Check & fix all functions below
GenerateRefData_WithGarch <- function(garch,seedstart=1,numseries=1000,tsamples=5000,
                                      saveas="RefData_WithGarch.RDS") {
  
  ###
  #   Note: R uses column-major as the default order for matrices
  #   Therefore, it makes programming/life easier if we build our reference data this way!
  ###
  
  par_o <- garch$pars[1]
  if (garch$type==3) par_a <- 0 else par_a <- garch$pars[2]
  par_b <- garch$pars[3]
  if (garch$type==2) par_d <- garch$pars[4] else par_d <- 0
  
  refData <- matrix(NA,tsamples,numseries)
  
  for (b in seq(1,numseries))
  {
    set.seed(seedstart+b)
    z <- rnorm(tsamples)
    discard <- 2000    #Discard the first 1000 iterations of h(t) to get a good starting value
    z_discard <- rnorm(discard)
    # recursive ht and et
    sim_e_garch <- rep(0,tsamples)
    # Set the first data points, then loop through the remainder
    ht_1 <- ht <- 1
    et_1 <- et <- sqrt(ht)*z_discard[1]
    for (t in 2:discard) {
      ht_1 <- ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      et_1 <- et <- sqrt(ht)*z_discard[t]
    }
    #
    ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
    sim_e_garch[1] <- sqrt(ht)*z[1]
    for (t in 2:tsamples) {
      ht_1 <- ht <- par_o + par_a*(sim_e_garch[t-1])^2 + par_b*ht_1 + par_d*(min(sim_e_garch[t-1],0)^2)
      sim_e_garch[t] <- sqrt(ht)*z[t]
    }
    refData[,b] <- sim_e_garch
  } # End of for (b in seq(1,numloops))
  
  saveRDS(refData,saveas)
  #Return:
  #rm(refData)
  refData
  
} # End: GenerateRefData_WithGarch()


GenerateRefData_WithTVGarch <- function(garch,tv,seedstart=0,numloops=5000,tsamples=5000,saveas="GeneratedRefData_WithTVGarch.RDS")
{
  
  ###
  #   Note: R uses column-major as the default order for matrices
  #   Therefore, it makes programming/life easier if we build our reference data this way!
  ###
  
  if (garch$type==1){
    par_o <- garch$pars[1]
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- 0
  }else if (garch$type==2){
    par_o <- garch$pars[1]
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- garch$pars[4]
  }else if (garch$type==3){
    par_o <- calculateOmega(garch$pars[2],garch$pars[3],garch$pars[4])
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- garch$pars[4]
  }
  
  TV <- tv
  refData <- matrix(NA,tsamples,numloops)
  
  #z    # this will store iid data
  #e    # this will store data with GARCH
  #eps  # this will store correlated data with GARCH and TV
  
  for (b in seq(1,numloops))
  {
    set.seed(seedstart+b)
    z <- rnorm(tsamples)
    if (is.null(GARCH)) e<-z
    else{
      discard <- round(tsamples*0.25)    #Discard the first 25% of the data set
      z_discard <- rnorm(discard)
      # recursive ht and et
      e <- rep(0,tsamples)
      # Set the first data points, then loop through the remainder
      ht_1 <- ht <- 1
      et_1 <- et <- sqrt(ht)*z_discard[1]
      for (t in 2:discard) {
        ht_1 <- ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
        et_1 <- et <- sqrt(ht)*z_discard[t]
      }
      #
      ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      e[1] <- sqrt(ht)*z[1]
      for (t in 2:tsamples) {
        ht_1 <- ht <- par_o + par_a*(e[t-1])^2 + par_b*ht_1 + par_d*(min(e[t-1],0)^2)
        e[t] <- sqrt(ht)*z[t]
      }
    } # GARCH data created => e complete
    
    if (is.null(TV)) eps <- e
    else{
      # - - - NO TV at all (no delta0)- - -
      if (is.null(TV$delta0)){
        gt <- 1
        eps <- e
      } else if (NROW(TV$pars) == 1) {
        # - - - NO TV (only constant delta0)- - -
        gt <- TV$delta0
        eps <- e*sqrt(gt)
      } else {
        # - - - Full TV - - -
        gt <- calculate_g(TV)
        eps <- e*sqrt(gt)
      }
    } # TV data created -> eps complete
    
    refData[,b] <- eps
  } # End of for (b in seq(1,numloops))
  
  saveRDS(refData,saveas)
  #Return:
  refData
  
}  # End: GenerateRefData_WithTVGarch



#### ==================  Data Generation  ================== ###


