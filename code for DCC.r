

## -- generateRefData -- ####
setGeneric(name="generateRefData",
           valueClass = "matrix",
           signature = c("nr.series","nr.obs","tvObj","garchObj","corrObj"),
           def =  function(nr.series,nr.obs,tvObj,garchObj,corrObj)
           {
             ## TODO: Improve function to handle TV & GARCH individually
             refData <- matrix(NA,nrow = nr.obs, ncol = nr.series)

             # Step1: Generate iid Data & create Correlation
             # u = un-correlated data
             u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)

             # e = correlated data (defaults to 'u' if there is no correlation object)
             e <- u
			 
			 
			corrType <- class(corrObj)
			 
			 
			 # - - - DCC - - -
			 # pass in Qbar matrix and "a" alpha & "b" beta for DCC
			 # need to create more than just nr.obs, add discard amount then drop the discard part before returning
			 discard <- 2000
             discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
             u <- rbind(discardData,u)
			 e <- u
			 endRow <- discard + nr.obs
			 Qt_1 <- Qbar
			 for (t in 2:endRow)){	# starting from t=2! (set Qt[1]=Qbar)
				Qt <- (1-a-b)*Qbar + a*e[t-1,,drop=FALSE]%*%t(e[t-1,,drop=FALSE]) + b*Qt_1 # N x N
				#scale Qt by inverse of its sqrt diagonals (from front and back) to make it correlation matrix
				Pt <- diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)%*%Qt%*%diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)  # N x N
				# create DCC correlated data
				Pt.sqrt <- sqrt_mat1(Pt)
				e[t,] <- t( mPt.sqrt %*% t(u[t,,drop=FALSE]) )
				# for the next loop round, store the "previous" Qt
				Qt_1 <- Qt 
			 }
			#return e
			 startRow <- discard + 1
             refData <- e[(startRow:endRow), ]
             



             # - - - CCC - - -
             if (corrType[1] == "ccc_class"){
               if(is.null(corrObj$Estimated)) {
                 P <- corrObj$P
               } else P <- corrObj$Estimated$P
               P.sqrt <-  sqrt_mat1(P)
               e <- t(P.sqrt %*% t(u))
             }

             if (corrType[1] == "stcc1_class"){
               if(is.null(corrObj$Estimated)) {
                 corrObj$Estimated$P1 <- corrObj$P1
                 corrObj$Estimated$P2 <- corrObj$P2
                 corrObj$Estimated$pars <- corrObj$pars
               }
               Pt <- .calc.Pt(corrObj)
               for (t in 1:corrObj@Tobs){
                 mPt <- unVecL(Pt[t,,drop=FALSE])
                 mPt.sqrt <- sqrt_mat1(mPt)
                 e[t,] <- t( mPt.sqrt %*% t(u[t,,drop=FALSE]) )
               }
             }

             #End: Generate Correlated Data

             # Step2: Inject GARCH into Data

             garchObj$pars["omega",1] <- ( 1 - garchObj$pars["alpha",1] - garchObj$pars["beta",1] )
             discard <- 2000
             discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
             refData <- rbind(discardData,e)
             endRow <- discard + nr.obs

             for (b in 1:nr.series){
               w <- z <- refData[,b]
               ht_1 <- 1
               w[1] <- z[1]
               for (t in 2:endRow) {
                 ht <- garchObj$pars["omega",1] + garchObj$pars["alpha",1]*(w[t-1])^2 + garchObj$pars["beta",1]*ht_1
                 if(garchObj$type == garchtype$gjr) { ht <- ht + garchObj$pars["gamma",1]*(min(w[t-1],0)^2) }
                 ht_1 <- ht
                 w[t] <- sqrt(ht)*z[t]
               }
               refData[,b] <- as.numeric(w)
             }

             # Discard the first 2000
             startRow <- discard + 1
             refData <- refData[(startRow:endRow), ]

             # Step3: Inject TV into Data
             gt <- get_g(tvObj)
             refData <- refData*sqrt(gt)

             #Return:
             refData

           }
)


# Simulation code

# N=2,5,10 series, Tobs=1000 length, R=500 replications

# rho = 0.0, 0.5, 0.9 used for Qbar matrix (constant matrix for DCC intercept)

# a = 0.100 b = 0.80 used for DCC dynamics parameters
# a = 0.050 b = 0.90 used for DCC dynamics parameters
# a = 0.025 b = 0.95 used for DCC dynamics parameters

# generate refData with Qbar = toeplitz(rho^seq.int(0, N-1)) and a and b.

# estimate CCC model

# test CCC vs TVCC using st=(1:Tobs)/Tobs


