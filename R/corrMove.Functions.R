##################################################
#DATA IMPORT
##################################################

# Currently this just fixes the timestamps.
as.corrData <- function(object, timeformat="",timezone="UTC")
{

  # fastPOSIXct doesn't have a timeformat argument and as.POSIXct doesn't accept this argument if empty/NA/NULL ???
  if(class(object$timestamp)=="character" && timeformat=="") { TIMES <- fasttime::fastPOSIXct(object$timestamp,tz=timezone) }
  else if(timeformat=="") { TIMES <- as.POSIXct(object$timestamp,tz=timezone) }
  else { TIMES <- as.POSIXct(object$timestamp,tz=timezone,format=timeformat) }
  object$timestamp <- TIMES
  cdDF <- as.data.frame(object)
  #class(cdDF) <- "corrData"
  return(cdDF)
}



##################################################
#MODEL ESTIMATION
##################################################

############################################################################
#Calculate the difference terms in the estimators for s and r
############################################################################
d <- function(dz, vn, ti) {dz - ti %o% vn}



############################################################################
#Calculate the sufficient statistic s from the data
############################################################################
s <-function(N, S, dz, vn, ti) {(1/(N*S)) * sum( d(dz, vn, ti)^2 / ti )}



############################################################################
#Calculate the sufficient statistic r from the data
############################################################################
r <- function(N, S, dz, vn, ti) {(1/(S)) * sum( (1/ti) *  ((1/N) * rowSums( d(dz, vn, ti) ))^2 )}



##################################################
#MODEL SELECTION
##################################################

############################################################################
#IC calculator for 00 model. IC=1 gives AICc, while IC=2 gives BIC
############################################################################
icMod00 <- function(N, S, C, sg, IC){

  n <- 2*N*S

  #Set IC to 1 for AICc, or 2 for BIC
  ic <- switch(IC,

               #AICc
               n*(log(sg) + n/(n-2)) - 2*C,

               #BIC
               log(N)*1 - 2*(C-(n/2)*(log(sg) + 1))
  )

  return(ic)
}



############################################################################
#IC calculator for U0 model. IC=1 gives AICc, while IC=2 gives BIC
############################################################################
icModU0 <- function(N, S, C, sg, sgU, IC){

  n <- 2*N*S

  ic <- switch(IC,

               #AICc
               n*(log(sgU) + ((n+2)*(n-2))/(n*(n-4))) - 2*C,

               #BIC
               log(N)*3 - 2*(C - (n/2)*(log(sgU) + sg/sgU))
  )

  return(ic)
}



############################################################################
#IC calculator for 0D model. IC=1 gives AICc, while IC=2 gives BIC
############################################################################
icMod0D <- function(N, S, C, loD, lpD, IC){

  n0 <- 2*S
  nP <- 2*(N-1)*S

  ic <- switch(IC,

               #AICc
               n0*(log(loD) + n0/(n0-2)) + nP*(log(lpD) + nP/(nP-2)) - 2*C,

               #BIC
               log(N)*2 - 2*(C - (n0/2)*(log(loD) + 1) - (nP/2)*(log(lpD) + 1))

  )

  return(ic)
}



############################################################################
#IC calculator for UD model. IC=1 gives AICc, while IC=2 gives BIC
############################################################################
icModUD <- function(N, S, C, loD, lpD, loDU, lpDU, IC){

  n0 <- 2*S
  nP <- 2*(N-1)*S

  ic <- switch(IC,

               #AICc
               n0*(log(loDU)+((n0+2)*(n0-2))/(n0*(n0-4)))+nP*(log(lpDU)+nP/(nP-2))-2*C,

               #BIC
               log(N)*4-2*(C-(n0/2)*(log(loDU) + loD/loDU)-(nP/2)*(log(lpDU)+lpD/lpDU))

  )

  return(ic)

}



##################################################
#INDEX ESTIMATION
##################################################

############################################################################
#Calculate the bias of an MCI estimator
############################################################################
biasEta <- function(eta, P, biasP, M, biasM, covPM, varM){
  eta*(biasP/P - biasM/M - covPM/(P*M) + varM/M^2)
}



############################################################################
#MCI bias correction and interval estimation based on the Fisher transform
############################################################################
zTrans <- function(eta, etaBias, varEta){

  z <- atanh(eta)
  dz <- -(1-eta^2)*etaBias

  varzdz <- varEta * (1/(1-eta^2) + 2 * eta * etaBias)^2

  zCiL <- z + dz - 1.96*sqrt(varzdz)
  zCiU <- z + dz + 1.96*sqrt(varzdz)

  etaDb <- tanh(z + dz)
  etaDbCiL <- tanh(zCiL)
  etaDbCiU <- tanh(zCiU)

  return(c(etaDb, etaDbCiL, etaDbCiU))

}



############################################################################
#Estimate the bias-corrected set of MCIs with CIs from the 00 model
############################################################################
eta00 <- function(msRes, cntTm=NA){

  #Everything is zero in this model
  etaDif <- c(0, 0, 0)
  etaDft <- c(0, 0, 0)
  etaTot <- c(0, 0, 0)

  #Package results for return
  #return(data.frame(cbind(etaDif, etaDft, etaTot)))
  #return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))
  if(is.na(cntTm)){
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))
  } else {
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes, cntTm))
  }

}



############################################################################
#Estimate the bias-corrected set of MCIs with CIs from the U0 model
############################################################################
etaU0 <- function(N, S, mt, Tt, vU2, sgU, msRes, cntTm=NA){

  #Index num and denom components
  Pdif <- 0
  Pdft <- mt*vU2
  M <- 2*sgU + mt*vU2

  #Index component covariance Jacobian
  cmp.jac <- matrix(c(0, mt, 2, mt), nrow=2, ncol=2, byrow=T)

  #Model degrees of freedom
  d <- 2*(N*S-1)

  #Model component variances
  varSgU <- 2 * sgU^2 / d
  varVU2 <- 4 * vU2 * sgU / (N * Tt) + 4 * (sgU / (N * Tt))^2

  #Model sufficient statistic variances
  sfs.var <- matrix(c(varSgU, 0, 0, varVU2), nrow=2, ncol=2, byrow=T)

  #Index estimator component covariances by applying the Jacobian to the suff   stat vars
  cmp.cov <- cmp.jac %*% sfs.var %*% t(cmp.jac)

  #Biases of index estimator components
  biasPdif <- 0
  biasPdft <- 2*(mt/(N*Tt))*sgU
  biasM <- 2*(mt/(N*Tt))*sgU

  #Biased index estimators
  etaDifBiased <- 0
  etaDftBiased <- Pdft/M
  etaTotBiased <- etaDftBiased

  #Jacobian for transforming index estimator components into index covariances
  ind.jac <- (1/M)*matrix(c(1, -etaDftBiased, 1, -etaTotBiased), nrow=2, ncol
                          =2, byrow=T)

  #Index covariances
  ind.cov <- ind.jac %*% cmp.cov %*% t(ind.jac)

  #etaDif estimator bias
  biasEtaDif <- 0

  #etaDft estimator bias
  biasEtaDft <- biasEta(etaDftBiased, Pdft, biasPdft, M, biasM, cmp.cov[1, 2],
                        cmp.cov[2,2])

  #etaTot estimator bias
  biasEtaTot <- biasEtaDft

  #Debiased index estimates plus 95% CI end points
  etaDif <- c(0, 0, 0)
  etaDft <- zTrans(etaDftBiased, biasEtaDft, ind.cov[1,1])
  etaTot <- zTrans(etaTotBiased, biasEtaTot, ind.cov[2,2])

  #Package results for return
  if(is.na(cntTm)){
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))
  } else {
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes, cntTm))
  }

}



############################################################################
#Estimate the bias-corrected set of MCIs with CIs from the 0D model
############################################################################
eta0D <- function(N, S, loD, lpD, msRes, cntTm=NA){

  #Index num and denom components
  Pdif <- (2/N)*(loD - lpD)
  Pdft <- 0
  M <- (2/N)*(loD + (N - 1)*lpD)

  #Index component covariance Jacobian
  cmp.jac <- matrix(c(2/N, -2/N, 2/N, 2*(N-1)/N), nrow=2, ncol=2, byrow=T)

  #Model degrees of freedom
  no <- 2*S
  np <- 2*(N-1)*S

  #Model component variances
  varLoD <- 2 * loD^2 / no
  varLpD <- 2 * lpD^2 / np

  #Model sufficient statistic variances
  sfs.var <- matrix(c(varLoD, 0, 0, varLpD), nrow=2, ncol=2, byrow=T)

  #Index estimator component covariances by applying the Jacobian to the        sufficient statistic vars
  cmp.cov <- cmp.jac %*% sfs.var %*% t(cmp.jac)

  #Biases of index estimator components
  biasPdif <- 0
  biasPdft <- 0
  biasM <- 0

  #Biased index estimators
  etaDifBiased <- Pdif/M
  etaDftBiased <- 0
  etaTotBiased <- etaDifBiased

  #Jacobian for transforming index estimator components into index covariances
  ind.jac <- (1/M)*matrix(c(1, -etaDifBiased, 1, -etaTotBiased), nrow=2, ncol
                          =2, byrow=T)

  #Index covariances
  ind.cov <- ind.jac %*% cmp.cov %*% t(ind.jac)

  #etaDif estimator bias
  biasEtaDif <- biasEta(etaDifBiased, Pdif, biasPdif, M, biasM, cmp.cov[1, 2],
                        cmp.cov[2,2])

  #etaDft estimator bias
  biasEtaDft <- 0

  #etaTot estimator bias for DU model
  biasEtaTot <- biasEtaDif

  #Debiased index ests plus 95% CI end points
  etaDif <- zTrans(etaDifBiased, biasEtaDif, ind.cov[1,1])
  etaDft <- c(0, 0, 0)
  etaTot <- zTrans(etaTotBiased, biasEtaTot, ind.cov[2,2])

  #Package results for return
  #return(data.frame(cbind(etaDif, etaDft, etaTot)))
  #return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))

  if(is.na(cntTm)){
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))
  } else {
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes, cntTm))
  }

}



############################################################################
#Estimate the bias-corrected set of MCIs with CIs from the UD model
############################################################################
etaUD <- function(N, S, mt, Tt, vU2, loDU, lpDU, msRes, cntTm=NA){

  #Index num and denom components
  Pdif <- (2/N)*(loDU - lpDU)
  Pdft <- mt*vU2
  M <- (2/N)*(loDU + (N - 1)*lpDU) + mt*vU2

  #Index component covariance Jacobian
  cmp.jac <- matrix(c(2/N, -2/N, 0, 0, 0, mt, 2/N, 2*(N-1)/N, mt), nrow=3,
                    ncol=3, byrow=T)

  #Model degrees of freedom
  no <- 2*S
  np <- 2*(N-1)*S

  #Model component variances
  varLoDU <- 2 * loDU^2 / (no - 2)
  varLpDU <- 2 * lpDU^2 / np
  varVU2 <- 4 * vU2 * loDU / (N * Tt) + 4 * (loDU / (N * Tt))^2

  #Sufficient statistic variances
  sfs.var <- matrix(c(varLoDU, 0, 0, 0, varLpDU, 0, 0, 0, varVU2), nrow=3,
                    ncol=3, byrow=T)

  #Index estimator component covariances by applying the Jacobian to the suff   stat vars
  cmp.cov <- cmp.jac %*% sfs.var %*% t(cmp.jac)

  #Biases of index estimator components
  biasPdif <- 0
  biasPdft <- 2*(mt/(N*Tt))*loDU
  biasM <- 2*(mt/(N*Tt))*loDU

  #Biased index estimators for DU model
  etaDifBiased <- Pdif/M
  etaDftBiased <- Pdft/M
  etaTotBiased <- etaDifBiased + etaDftBiased

  #Jacobian for transforming index estimator components into index covariances
  ind.jac <- (1/M)*matrix(c(1, 0, -etaDifBiased, 0, 1, -etaDftBiased, 1, 1,
                            -etaTotBiased), nrow=3, ncol=3, byrow=T)

  #Index covariances
  ind.cov <- ind.jac %*% cmp.cov %*% t(ind.jac)

  #etaDif estimator bias
  biasEtaDif <- biasEta(etaDifBiased, Pdif, biasPdif, M, biasM, cmp.cov[1, 3],
                        cmp.cov[3,3])

  #etaDft estimator bias
  biasEtaDft <- biasEta(etaDftBiased, Pdft, biasPdft, M, biasM, cmp.cov[2, 3],
                        cmp.cov[3,3])

  #etaTot estimator bias for DU model
  biasEtaTot <- biasEta(etaDftBiased, Pdif + Pdft, biasPdif + biasPdft, M,
                        biasM, cmp.cov[1, 3] + cmp.cov[2, 3], cmp.cov[3,3])

  #Debiased index estimates and 95% CI end points
  etaDif <- zTrans(etaDifBiased, biasEtaDif, ind.cov[1,1])
  etaDft <- zTrans(etaDftBiased, biasEtaDft, ind.cov[2,2])
  etaTot <- zTrans(etaTotBiased, biasEtaTot, ind.cov[3,3])

  #Package results for return
  if(is.na(cntTm)){
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes))
  } else {
    return(list(data.frame(cbind(etaDif, etaDft, etaTot)), msRes, cntTm))
  }

}


############################################################################
#A function to perform a moving window analysis with model selection
#within each window
############################################################################
# mciWindow <- function(dataAll, dateFormat, W, IC=1){
#
#   #Cast times as POSIXct
#   if(class(dataAll[,1])[1]!="POSIXct"){
#     dataAll[,1] <- as.POSIXct(strptime(dataAll[, 1], dateFormat))
#   }
#
#   #Extract vector of all times in the data
#   tmsAll <- dataAll[,1]
#
#   #Get the start date
#   t0 <- tmsAll[1]
#
#   #Total number of observations in the dataset
#   nObs <- nrow(dataAll)
#
#   #Vector to store model selection results
#   msRes <- NULL
#
#   #Counter to keep track of window number
#   cntWins <- 0
#
#   #List for storing the mci results
#   mcis <- list()
#
#   #Set the window width in terms of number of observations in the data,
#   #nObs: 0 < W <= nObs. To do a full dataset analysis, set W = nObs.
#   #For a moving window analysis, choose W such that W < nObs.
#   numWins <- nObs - (W - 1)
#
#   #cnt <- as.POSIXct(rep(NA, numWins))
#
#   #Stop execution if an inappropriate window size is used
#   if(W > nObs) {
#     stop("The window size W cannot be greater than the number of observations")
#   } else if(W < 5) {
#     stop("The window size W must be at least 5")
#   } else {
#
#     #Loop over all possible windows in the dataset
#     for(i in W:nObs){
#
#       #Get the positions of the lower and upper limits of the current
#       #window
#       lLim <- i - (W - 1)
#       uLim <- i
#
#       #Extract the current window from the data. This will include rows that
#       #are missing data (coded by NAs)
#       dataWin <- dataAll[lLim:uLim, ]
#       data <- dataWin[complete.cases(dataWin), ]
#
#       #Increment the window counter
#       cntWins <- cntWins + 1
#
#       #Get the dimensions of the dataset
#       nCol <- ncol(data)
#       N <- (nCol-1)/2
#       S <- nrow(data) - 1
#
#       #Times are expected to be in POSIXct format
#       #tms <- as.POSIXct(strptime(data[ , 1], dateFormat)) #, tz="EDT"
#       tms <- data[ , 1]
#       ti <- diff(tms)
#       units(ti) <- "secs"
#       ti <- as.numeric(ti)
#
#       #Get the regular timestep from the data
#       mt <- median(ti)
#
#       #Calculate total T
#       Tt <- sum(ti)
#
#       #Calculate the position of the center of the current window
#       winCnt <- i - (W - 1)/2
#
#       td1 <- diff(c(t0, tmsAll[floor(winCnt)]))
#       td2 <- diff(c(t0, tmsAll[ceiling(winCnt)]))
#       td <- (td1 + td2)/2
#       cntTm <- t0 + td
#       #cnt[cntWins] <- cntTm
#
#       #Get the X and Y positions for all individuals
#       x <- data[, 2:(N+1)]
#       y <- data[, (N + 2):(2*N + 1)]
#
#       #Calculate the x and y position shifts for all individuals
#       x1 <- x[1:S,]
#       x2 <- x[2:(S + 1),]
#       dx <- x2 - x1
#
#       y1 <- y[1:S,]
#       y2 <- y[2:(S + 1),]
#       dy <- y2 - y1
#
#       #Mean vector is just 0's for the no drift models
#       v0 <- rep(0, N)
#
#       #Estimate the X and Y uniform drift values
#       vxU <- (1/(N*Tt)) * c(sum( dx ))
#       vyU <- (1/(N*Tt)) * c(sum( dy ))
#
#       #Calculate the squared drift
#       vU2 <- vxU^2 + vyU^2
#
#       #Create vectors of the X and Y uniform drift values of length N
#       vyU <- rep(vyU, N)
#       vxU <- rep(vxU, N)
#
#       #Estimate the sufficient statistics s and r for the no drift models
#       s0 <- (1/2)*(s(N, S, dx, v0, ti) + s(N, S, dy, v0, ti))
#       r0 <- (1/2)*(r(N, S, dx, v0, ti) + r(N, S, dy, v0, ti))
#
#       #Estimate the sufficient statistics s and r for the uniform drift              models
#       sU <- (1/2)*(s(N, S, dx, vxU, ti) + s(N, S, dy, vyU, ti))
#       rU <- (1/2)*(r(N, S, dx, vxU, ti) + r(N, S, dy, vyU, ti))
#
#
#       ##################################################
#       #Model-specific parameter estimates
#
#       #Lambda ests for correlated diffusion, uniform drift
#       loDU <- N*S*rU / (S - 1)
#       lpDU <- N*(sU - rU) / (N - 1)
#
#       #Variance est for uncorrelated diffusion, uniform drift
#       sgU <- (N*S / (N*S-1)) * sU
#
#       #Lambda ests for correlated diffusion, no drift
#       loD <- N*r0
#       lpD <- N*(s0 - r0) / (N - 1)
#
#       #Variance est for uncorrelated diffusion, no drift
#       sg <- s0
#       ##################################################
#
#
#       ##################################################
#       #Model selection
#
#       #Model independent constant
#       C <- -N*sum(log(2*pi*ti))
#
#       #No drift, uncorrelated diffusion
#       ic00 <- icMod00(N, S, C, sg, IC)
#
#       #Uniform drift, uncorrelated diffusion
#       icU0 <- icModU0(N, S, C, sg, sgU, IC)
#
#       #No drift, correlated diffusion
#       ic0D <- icMod0D(N, S, C, loD, lpD, IC)
#
#       #Uniform drift, correlated diffusion
#       icUD <- icModUD(N, S, C, loD, lpD, loDU, lpDU, IC)
#
#       #Identify the selected model
#       icVals <- c(ic00, icU0, ic0D, icUD)
#       msRes <- which.min(icVals)
#       #msRes <- 4
#       #Contingent on model selection, choose appropriate index estimators
#       ests <- switch(msRes,
#
#                      #No drift, uncorrelated diffusion
#                      eta00(msRes, cntTm),
#
#                      #Uniform drift, uncorrelated diffusion
#                      etaU0(N, S, mt, Tt, vU2, sgU, msRes, cntTm),
#
#                      #No drift, correlated diffusion
#                      eta0D(N, S, loD, lpD, msRes, cntTm),
#
#                      #Uniform drift, correlated diffusion
#                      etaUD(N, S, mt, Tt, vU2, loDU, lpDU, msRes, cntTm)
#
#       )
#
#       #Store results for each window in a list indexed by window number
#       mcis[[cntWins]] <- ests
#
#
#     } #end for
#
#   } #end W check if
#
#   return(mcis)
#
# } #end mciWindow



############################################################################
#A function that returns an information criterion for associate with
#a given piece of data. Performs model selection on the given piece
#of data, and gives the result from the lowest IC model
#Can choose between AICc (IC=1) or BIC (IC=2)
############################################################################
getIC <- function(data, IC=1){

  #Eliminate any row that has an NA anywhere in it
  data <- data[complete.cases(data), ]

  #Check that W is still respected after removing NAs
  #NOTE: THIS FUNCTION SHOULD BE REVISED TO EXPLICITLY ACCEPT W
  #AS AN ARGUMENT
  if(dim(data)[1] >= 2){

    #Get the dimensions of the no-NA dataset
    nCol <- ncol(data)
    N <- (nCol-1)/2
    S <- nrow(data) - 1

    #Times are expected to be in POSIXct format
    #if(class(data[,1])[1]!="POSIXct"){
     # data[,1] <- as.POSIXct(strptime(data[, 1], dateFormat))
    #}

    tms <- data[ ,1]

    ti <- diff(tms)
    units(ti) <- "secs"
    ti <- as.numeric(ti)

    #Get the regular timestep from the data
    mt <- median(ti)

    #Calculate total T
    Tt <- sum(ti)

    #Calculate the position of the center of the current window
    #winCnt <- i - (W - 1)/2

    #td1 <- diff(c(t0, tmsAll[floor(winCnt)]))
    #td2 <- diff(c(t0, tmsAll[ceiling(winCnt)]))
    #td <- (td1 + td2)/2
    #cntTm <- t0 + td
    #cnt[cntWins] <- cntTm

    #Get the X and Y positions for all individuals
    x <- data[, 2:(N+1)]
    y <- data[, (N + 2):(2*N + 1)]

    #Calculate the x and y position shifts for all individuals
    x1 <- x[1:S,]
    x2 <- x[2:(S + 1),]
    dx <- x2 - x1

    y1 <- y[1:S,]
    y2 <- y[2:(S + 1),]
    dy <- y2 - y1

    #Mean vector is just 0's for the no drift models
    v0 <- rep(0, N)

    #Estimate the X and Y uniform drift values
    vxU <- (1/(N*Tt)) * c(sum( dx ))
    vyU <- (1/(N*Tt)) * c(sum( dy ))

    #Calculate the squared drift
    vU2 <- vxU^2 + vyU^2

    #Create vectors of the X and Y uniform drift values of length N
    vyU <- rep(vyU, N)
    vxU <- rep(vxU, N)

    #Estimate the sufficient statistics s and r for the no drift models
    s0 <- (1/2)*(s(N, S, dx, v0, ti) + s(N, S, dy, v0, ti))
    r0 <- (1/2)*(r(N, S, dx, v0, ti) + r(N, S, dy, v0, ti))

    #Estimate the sufficient statistics s and r for the uniform drift              models
    sU <- (1/2)*(s(N, S, dx, vxU, ti) + s(N, S, dy, vyU, ti))
    rU <- (1/2)*(r(N, S, dx, vxU, ti) + r(N, S, dy, vyU, ti))

    ##################################################
    #Model-specific parameter estimates

    #Lambda ests for correlated diffusion, uniform drift
    loDU <- N*S*rU / (S - 1)
    lpDU <- N*(sU - rU) / (N - 1)

    #Variance est for uncorrelated diffusion, uniform drift
    sgU <- (N*S / (N*S-1)) * sU

    #Lambda ests for correlated diffusion, no drift
    loD <- N*r0
    lpD <- N*(s0 - r0) / (N - 1)

    #Variance est for uncorrelated diffusion, no drift
    sg <- s0
    ##################################################


    ##################################################
    #Model selection

    #Model independent constant
    C <- -N*sum(log(2*pi*ti))

    #No drift, uncorrelated diffusion
    ic00 <- icMod00(N, S, C, sg, IC)

    #Uniform drift, uncorrelated diffusion
    icU0 <- icModU0(N, S, C, sg, sgU, IC)

    #No drift, correlated diffusion
    ic0D <- icMod0D(N, S, C, loD, lpD, IC)

    #Uniform drift, correlated diffusion
    icUD <- icModUD(N, S, C, loD, lpD, loDU, lpDU, IC)


    #Identify the selected model
    icVals <- c(ic00, icU0, ic0D, icUD)
    #aicVals <- c(aic00, aicU0)
    msRes <- which.min(icVals)
    #msRes <- 4

    return(icVals[msRes])

  } else {

    #If W is violated after removing NAs, return NA
    return(NA)
  }
}



############################################################################
#A function to find the single best partition within a chunk of data
#based on IC differences. Returns NA if no such partition exists.
############################################################################

icPart <- function(data, W=2, IC=1){

  #Get the length of the dataset passed to the function
  nPart <- nrow(data)

  #Check if it is possible to introduce a partition without violating
  #the window constraint. If not, return NA
  if(nPart >= (2*W-1)){

    pCnt <- 0
    partIC <- NULL

    #Loop over all possible partition points
    for(j in W:(nPart-W+1)){

      pCnt <- pCnt + 1

      #Check if there is any missing data at the focal partition point.
      #If so, skip it and return NA, if not, evaluate it.
      if(sum(!is.na(data[j, ]))==length(data[j, ])){

        #Split the data into two pieces at the focal partition point j
        part1 <- data[1:j, ]
        part2 <- data[j:nPart, ]

        #Get the IC vals for each piece and add them together

        ic1 <- getIC(part1, IC)
        ic2 <- getIC(part2, IC)

        partIC[pCnt] <- ic1 + ic2

      } else {

        partIC[pCnt] <- NA

      }

    } #end j for

    #Find the lowest IC value associated with the candidate partition points.
    sPartLoc <- which.min(partIC)

    minIC <- partIC[sPartLoc]

    #Check if the lowest IC partition point is better than leaving the        #data chunk intact.
    if((getIC(data, IC) - minIC) > 0){

      #If the best partition point beats the intact data chunk, return its
      #abosulte location and associated IC val.
      loc <- (W-1) + sPartLoc
      return(data.frame(loc=loc, ic=minIC))

    } else {

      #If it is not possible to beat the IC of the intact data chunk,
      #return NAs.
      return(data.frame(loc=NA, ic=NA))

    } #end minAIC if


  } else {

    #If it is not possible to introduce a partition point without
    #violating W, return NAs.
    return(data.frame(loc=NA, ic=NA))

  } #end nPart if

}



############################################################################
#A function to find partitions in the dataset via IC differences
#and a window-based stopping rule
############################################################################
findPrts <- function(data, W, IC=1){

  #Set the absolute minimum window width to the level that, below which
  #models can no longer be fit due to insufficient data
  MW <- 2

  #Get the length of the dataset
  dNd <- dim(data)[1]

  #Start with full dataset IC, then modify to reflect changes via partitioning
  dIcNew <- getIC(data, IC)
  dIcOld <- dIcNew + 1

  #Give the start and end of the dataset as the initial "partition points"
  pPS <- c(1, dNd)

  while(dIcNew < dIcOld){

    #Update the looping criterion
    dIcOld <- dIcNew

    cIC <- NULL
    cPS <- NULL
    cMW <- NULL

    #Loop over partitions
    for(i in 1:(length(pPS)-1)){

      #Get the start and end points of the current partitions
      pSt <- pPS[i]
      pNd <- pPS[i+1]

      #Try to partition the current chunk of data
      prt <- icPart(data[pSt:pNd, ], MW, IC)

      #If the attempted partitioning worked (ie, did not produce NAs)
      #Update the IC diffs, partition location lists, and min width of
      #the resulting new partition.
      #If the partitioning failed, return all NAs.
      if(!is.na(prt$ic)){
        icFC <- getIC(data[pSt:pNd, ], IC)
        cIC <- c(cIC, icFC - prt$ic)
        cps <- pSt + prt$loc - 1
        cPS <- c(cPS, cps)
        cMW <- c(cMW, min(c(pNd-cps+1, prt$loc)))
      } else {
        cIC <- c(cIC, NA)
        cPS <- c(cPS, NA)
        cMW <- c(cMW, NA)
      }

    }

    #Find the biggest AIC difference
    #PROBABLY NEED TO DEAL WITH CASE OF MULTIPLE IDENTICAL MAX DELTA IC VALS     #HERE
    mIC <- which.max(cIC)

    #If this position if not an NA, and the resulting partitions do not
    #violate the min window size, then accept the new partition and
    #update the partition locations list, and also the IC value for
    #the whole dataset.
    if(length(mIC) == 1){
      if(!is.na(cMW[mIC]) & (cMW[mIC] >= W)){
        pPS <- sort(c(pPS, cPS[mIC]))
        dIcNew <- dIcOld - cIC[mIC]
      } else {
        dIcNew <- dIcOld
      }
    } else {
      dIcNew <- dIcOld
    }

  }

  #Return the list of partition points
  return(pPS)

}



############################################################################
#A function to get index estimates for a partitioned dataset
#and format them for plotting
############################################################################
corrMove <- function(data, prts, IC=1){

  #Empty list to store results
  mcis <- list()
  prt <- NULL

  #Loop over all partitions
  for(i in 1:(length(prts)-1)){

    #Translate the partition locs into partition start and end points
    pSt <- prts[i]
    pNd <- prts[i+1]

    #Extract the current partition from the full dataset
    part <- data[pSt:pNd, ]

    #Eliminate any row that has an NA anywhere in it
    part <- part[complete.cases(part), ]

    #Get the dimensions of the dataset
    nCol <- ncol(part)
    N <- (nCol-1)/2
    S <- nrow(part) - 1

    #Times are expected to be in POSIXct format
    #Cast times as POSIXct if necessary
    #if(class(part[,1])[1]!="POSIXct"){
     # part[,1] <- as.POSIXct(strptime(part[, 1], dateFormat))
    #}

    tms <- part[ , 1]
    ti <- diff(tms)
    units(ti) <- "secs"
    ti <- as.numeric(ti)

    #Get the regular timestep from the data
    mt <- median(ti)

    #Calculate total T
    Tt <- sum(ti)

    #Get the X and Y positions for all individuals
    x <- part[, 2:(N+1)]
    y <- part[, (N + 2):(2*N + 1)]

    #Calculate the x and y position shifts for all individuals
    x1 <- x[1:S,]
    x2 <- x[2:(S + 1),]
    dx <- x2 - x1

    y1 <- y[1:S,]
    y2 <- y[2:(S + 1),]
    dy <- y2 - y1

    #Mean vector is just 0's for the no drift models
    v0 <- rep(0, N)

    #Estimate the X and Y uniform drift values
    vxU <- (1/(N*Tt)) * c(sum( dx ))
    vyU <- (1/(N*Tt)) * c(sum( dy ))

    #Calculate the squared drift
    vU2 <- vxU^2 + vyU^2

    #Create vectors of the X and Y uniform drift values of length N
    vyU <- rep(vyU, N)
    vxU <- rep(vxU, N)

    #Estimate the sufficient statistics s and r for the no drift models
    s0 <- (1/2)*(s(N, S, dx, v0, ti) + s(N, S, dy, v0, ti))
    r0 <- (1/2)*(r(N, S, dx, v0, ti) + r(N, S, dy, v0, ti))

    #Estimate the sufficient statistics s and r for the uniform drift models
    sU <- (1/2)*(s(N, S, dx, vxU, ti) + s(N, S, dy, vyU, ti))
    rU <- (1/2)*(r(N, S, dx, vxU, ti) + r(N, S, dy, vyU, ti))


    ##################################################
    #Model-specific parameter estimates

    #Lambda ests for correlated diffusion, uniform drift
    loDU <- N*S*rU / (S - 1)
    lpDU <- N*(sU - rU) / (N - 1)

    #Variance est for uncorrelated diffusion, uniform drift
    sgU <- (N*S / (N*S-1)) * sU

    #Lambda ests for correlated diffusion, no drift
    loD <- N*r0
    lpD <- N*(s0 - r0) / (N - 1)

    #Variance est for uncorrelated diffusion, no drift
    sg <- s0
    ##################################################


    ##################################################
    #Model selection

    #Model independent constant
    C <- -N*sum(log(2*pi*ti))

    #No drift, uncorrelated diffusion
    ic00 <- icMod00(N, S, C, sg, IC)

    #Uniform drift, uncorrelated diffusion
    icU0 <- icModU0(N, S, C, sg, sgU, IC)

    #No drift, correlated diffusion
    ic0D <- icMod0D(N, S, C, loD, lpD, IC)

    #Uniform drift, correlated diffusion
    icUD <- icModUD(N, S, C, loD, lpD, loDU, lpDU, IC)

    #Identify the selected model
    icVals <- c(ic00, icU0, ic0D, icUD)
    msRes <- which.min(icVals)

    #Contingent on model selection, choose appropriate index estimators
    ests <- switch(msRes,

                   #No drift, uncorrelated diffusion
                   eta00(msRes),

                   #Uniform drift, uncorrelated diffusion
                   etaU0(N, S, mt, Tt, vU2, sgU, msRes),

                   #No drift, correlated diffusion
                   eta0D(N, S, loD, lpD, msRes),

                   #Uniform drift, correlated diffusion
                   etaUD(N, S, mt, Tt, vU2, loDU, lpDU, msRes)

    )

    #Adjust the end point of the partitions appropriately so
    #that index values do not overlap at the partition points.
    pNd <- if(i == length(prts)-1){
      pNd
    } else  {
      pNd-1
    }

    #Assign the index values of the focal partition to each point
    #within the partition.
    for(j in pSt:pNd){
      #MCI and model selection results
      mcis[[j]] <- ests

      #Keep track of the partitions
      prt[j] <- i
    }


  }

  #Return the index estimates and CIs. Old output format. Retain for now for
  #Compatibility with old plotting code.
  #return(mcis)

  #New, human readable output code starts here. This could be simplified by
  #refactoring the code that builds up the early and intermediate results.

  #Vectors to temporarily store extracted results
  etaDifPt <-NULL
  etaDifLci <-NULL
  etaDifUci <-NULL

  etaDftPt <-NULL
  etaDftLci <-NULL
  etaDftUci <-NULL

  etaTotPt <-NULL
  etaTotLci <-NULL
  etaTotUci <-NULL

  selModNum <- NULL
  selModCode <- NULL

  for(k in 1:length(mcis)){

    #Unpack the mcis list and extract the components for each index. This is etaDif
    eDif <- mcis[[k]][[1]][[1]]

    #Store the different indices in different vectors
    etaDifPt[k] <-eDif[1]
    etaDifLci[k] <-eDif[2]
    etaDifUci[k] <-eDif[3]

    #As above, but for etaDft
    eDft <- mcis[[k]][[1]][[2]]

    etaDftPt[k] <-eDft[1]
    etaDftLci[k] <-eDft[2]
    etaDftUci[k] <-eDft[3]

    #As above but for etaTot
    eTot <- mcis[[k]][[1]][[3]]

    etaTotPt[k] <-eTot[1]
    etaTotLci[k] <-eTot[2]
    etaTotUci[k] <-eTot[3]

    #Extract the selected model and label as described in the paper.
    modNum <- mcis[[k]][[2]]
    selModNum[k] <- modNum
    modCode <-switch(modNum, "UU", "CU", "UC", "CC")
    selModCode[k] <- modCode

  }


  #Create a dataframe with all relevant info and a simple
  #structure: 1 row per timestamp.
  mcisOut <- data.frame(timestamp=data$timestamp, etaDif.MLE=etaDifPt,
      etaDif.CI.Low=etaDifLci, etaDif.CI.Upp=etaDifUci, etaDft.MLE=etaDftPt,
      etaDft.CI.Low=etaDftLci, etaDft.CI.Upp=etaDftUci, etaTot.MLE=etaTotPt,
      etaTot.CI.Low=etaTotLci, etaTot.CI.Upp=etaTotUci, sel.model.num=selModNum, sel.mod.code=selModCode, partition
      =prt)

  class(mcisOut) <- "corrMove"

  return(mcisOut)

}
############################################################################
############################################################################

########################################################################
#Plot method for corrMove objects
########################################################################
plot.corrMove <- function(x, ...){

  #Set a buffer for the y-axis range above the max and below the min
  buf <-0.0125

  with(x, {
    #Find the min y value over all panels
    minVal <- min(c(min(etaDif.CI.Low), min(etaDft.CI.Low), min(etaTot.CI.Low)))

    #Find the max y value over all panels
    maxVal <- max(c(max(etaDif.CI.Upp), max(etaDft.CI.Upp), max(etaTot.CI.Upp)))

    #Set the y-axis limits such that everything to be plotted will fit
    yMin <- minVal - buf
    yMax <- maxVal + buf

    #Define 4 colors for the 4 different models. This is the
    #ColorBrewer 'Dark2' palette for 4 categories
    cols <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a')

    #Set up a matrix for the different panels to be plotted.
    #The 4th panel is a dummy plot that just holds the legend
    m <- matrix(c(1, 2, 3, 4), nrow=4, ncol=1, byrow=TRUE)
    layout(mat=m, heights=c(0.3, 0.3,0.3,0.1))

    #Set the margins such that white space between panels is reduced.
    par(mar=c(2.5, 5, 1, 1))

    #Plot the point estimate and CIs over time for etaDif.
     plot(timestamp, etaDif.MLE, pch=20, col=cols[sel.model.num],
          ylim=c(yMin, yMax), ylab="Diffusive corr.", xlab="", cex.lab=1.5,
          cex.axis=1.5)
    points(timestamp, etaDif.CI.Low, pch="-", col=cols[sel.model.num])
    points(timestamp, etaDif.CI.Upp, pch="-", col=cols[sel.model.num])

    #Plot the point estimate and CIs over time for etaDft.
    plot(timestamp, etaDft.MLE, pch=20, col=cols[sel.model.num],
         ylim=c(yMin, yMax), ylab="Drift corr.", xlab="",  cex.lab=1.5,
         cex.axis=1.5)
    points(timestamp, etaDft.CI.Low, pch="-", col=cols[sel.model.num])
    points(timestamp, etaDft.CI.Upp, pch="-", col=cols[sel.model.num])

    #Plot the point estimate and CIs over time for etaTot.
    plot(timestamp, etaTot.MLE, pch=20, col=cols[sel.model.num],
         ylim=c(yMin, yMax), ylab="Total corr.", xlab="", cex.lab=1.5,
         cex.axis=1.5)
    points(timestamp, etaTot.CI.Low, pch="-", col=cols[sel.model.num])
    points(timestamp, etaTot.CI.Upp, pch="-", col=cols[sel.model.num])

    #Get the model numbers and codes of all models represented in the plot
    modNums <- unique(sel.model.num)
    modCodes <- as.character(unique(sel.mod.code))

    #Reset the margins for the last (dummy) panel.
    par(mar=c(0, 5, 0, 1))

    #Create a dummy plot that just holds the legend
    plot(1, type="n", axes=FALSE, xlab="", ylab="")

    #Add legend. Legend changes based on which models are actually represent
    #in the plot.
    legend("center", legend=modCodes, fill=cols[modNums], ncol=length(modNums),
           cex=1.75)
  }
  )

}
