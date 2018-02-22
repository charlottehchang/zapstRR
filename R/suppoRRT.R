### SUPPORT FUNCTIONS FOR RRT ANALYSES ###

#' zapstRR: A package for performing statistical inference on univariate and multivariate Randomized Response Technique (RRT) data
#'
#' zapstRR provides three categories of functions: RRT analysis, power analysis, and data simulation.
#'
#' @section Package Functions:
#' The RRT analysis functions are as follows:
#' \itemize{
#'      \item \code{\link{RRunivariate}}: Prevalence and logistic regression for single RRT items
#'      \item \code{\link{RRsumscore}}: Prevalence and ordinal regression for sum scores of RRT items
#'      \item \code{\link{RRirt}}: Prevalence and multivariate regression for two or more RRT items using item response theory
#'      \item \code{\link{RRbootstrap}}: Nonparametric bootstrap estimates of confidence intervals for the three models above.
#'      \item \code{\link{RRgof}}: Parametric bootstrap \eqn{\chi^2} goodness-of-fit estimates for the three RRT models.
#' }
#' In addition, \code{\link{powerRRT}} performs power analysis for study designs and \code{\link{RRsimulate}} can simulate example RRT data. Two datasets are provided: \code{\link{huntRRdata}}, a case study on illegal hunting in Southwest China, and \code{\link{RRdata}}, which investigates sensitive behaviors among university students.
#'
#' @aliases zapstRR-package zapstRR
#'
#' @details
#' \tabular{ll}{
#' Package: \tab zapstRR\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.0.900\cr
#' Date: \tab 2017-10-10\cr
#' Depends: \tab R (>= 3.3.0)\cr
#' Imports: \tab stats, grDevices, graphics\cr
#' Suggests: \tab knitr\cr
# Encoding: \tab UTF-8\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' URL: \tab \url{https://github.com/charlottehchang/zapstRR}\cr
#' Vignette: \tab \url{http://rpubs.com/chwchang/318110}\cr
#' }
#'
#' @docType package
#' @name zapstRR-package
#' @title zapstRR: ZoologicAl Package (Zap!) for Randomized Response Technique Analysis
#' @author Charlotte H. Chang \email{chchang@@nimbios.org} and Maarten J.L.F. Cruyff \email{m.cruyff@@uu.nl}
#' @keywords package
#'
#' @import stats
#'
NULL

#=============================================================================
# FUNCTION FOR SIMULATING RRT DATA (UNIVARIATE and MULTIVARIATE)
#=============================================================================

#' @title Simulate RRT data for three models: univariate, sum score, and item response theory
#'
#' @description \code{RRsimulate} simulates univariate (single item) and multivariate (multiple items) RRT data. Probabilities (\code{pi.1, p00, p11}) must be within the interval (0, 1).
#'
#' @param model A character string, one of \code{"univariate","sumscore",} or \code{"irt"}.
#' @param pi.1 A numeric vector with the true prevalence of sensitive trait(s) (each value must be within \code{c(0,1)}).
#' @param p00 Probability of answering \code{no} conditional on a true \code{no} (must be within c(0,1))
#' @param p11 Probability of answering \code{yes} conditional on a true \code{yes}
#' @param n Sample size (number of respondents)
#' @param x.betas A numeric vector of predictor variable regression coefficients (defaults to \code{c(0)})
#' @param seed See \code{\link[base]{Random}}.
#'
#' @examples
#' ## UNIVARIATE SIMULATION - with no informative predictor variable
#' UNIdata <- RRsimulate(model="univariate", pi.1 = 0.35, p00=5/6, p11=5/6, n=200)
#'        ## Simulate with an informative predictor
#' UNIdataLR <- RRsimulate(model="univariate",pi.1 = 0.35, p00=5/6, p11=5/6, x.betas=c(2), n=200)
#' head(UNIdataLR)
#'
#' ## SUM SCORE SIMULATION
#' pi.1s <- c(0.2, 0.35, 0.4)
#' p00 <- p11 <- 5/6
#' SSdata <- RRsimulate(model="sumscore",pi.1 = pi.1s, p00=p00, p11=p11, x.betas=c(-1, 3), n=200)
#' head(SSdata)
#'
#' ## ITEM RESPONSE THEORY SIMULATION
#' pi.1s <- c(0.15, 0.2, 0.3)
#' p00 <- p11 <- 5/6
#' IRTdata <- RRsimulate(model="irt", pi.1 = pi.1s, p00=p00, p11=p11, x.betas=c(-1, 3), n=200)
#' head(IRTdata)
#'
#' @return \code{RRsimulate} returns a data frame with y (RRT responses) and x (predictor variable) observations.
#'
#' @seealso \code{\link{RRunivariate}} for single-item RRT estimation, \code{\link{RRsumscore}} for sum score modeling, and \code{\link{RRirt}} for item randomized response theory modeling. The latter two models (sum score and IRT) are suitable for multivariate RRT data (two or more question prompts).
#' @export
RRsimulate <- function(model=c("univariate","sumscore","irt"), pi.1=c(0.2), p00=5/6, p11=5/6, n=100, x.betas=c(0), seed=1) {

  set.seed(seed)
  options(warn=-1)

  ## Check input
  if (!is.numeric(pi.1) | !is.numeric(p00) | !is.numeric(p11)) {
    stop("All pi.1, p00, and p11 values must be numeric and within the interval (0, 1).")
  }

  if (!is.numeric(x.betas)) {
    if (is.data.frame(x.betas)) {
      if (!is.numeric(as.matrix(x.betas))) {
        stop("All values for x.betas must be numeric.")
      }
    } else {
      stop("Specify a numeric vector for x.betas")
    }
  }

  pi.check <- findInterval(c(pi.1,p00,p11), c(0,1))
  pi.check <- ifelse(pi.check != 1, 0, 1)
  if ( sum(pi.check) != length(c(pi.1,p00,p11))) {
    stop("All pi.1, p00, and p11 values must be within the interval (0,1).")
  }

  ## Cond. classification matrix and predictor variable
  P <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)
  x.var <- replicate(length(x.betas), rnorm(n)) # generate predictor variables

  if (model=="Univariate" | model=="univariate") { # univariate RRT
    if (length(pi.1) > 1) {
      message("For the univariate model, pi.1 should only have one value, choosing first value.")
      pi.1 <- pi.1[1]
    }
    # Generate logistic probabilities for possessing sensitive trait
    alpha <- log(pi.1/(1-pi.1)) # intercept, should recapitulate pi.1 rate absent x predictors
    uni.prev <- exp(alpha + x.var%*%x.betas)/(1 + exp(alpha+ x.var%*%x.betas))
    uni.p <- cbind(1-uni.prev, uni.prev) %*% P
    # RRT process for single item
    y <- apply(uni.p, 1, function(x) {rmultinom(1, 1, x)})
    y <- y[2,] # observed yes'es
    uni.df <- data.frame(y=y, x=x.var)
    return(uni.df)
  } else if (model=="Sumscore" | model=="sumscore") { # sum score RRT
    if (length(pi.1) < 2) {
      message("pi.1 must have two or more values for the sumscore model, adding a random value")
      pi.1[2] <- runif(1)
    }
    pi.ordered <- sort(pi.1, index.return=TRUE) # Get pi.1 into increasing values
    pi.1 <- pi.ordered$x

    # Generate ordinal probabilities for 1..K sensitive traits
    b <- c(log(pi.1/(1-pi.1)),x.betas)
    alpha <- matrix(b[1:length(pi.1)], n, length(pi.1), byrow=T)
    cum.prev <- lg(alpha-c(x.var%*%x.betas))
    ss.p <- cbind(cum.prev, 1) - cbind(0, cum.prev)

    # Generate transition matrix for sum score RRT misclassification
    M  <- length(pi.1)
    P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)
      # Kronecker product
    Q <- P%x%P
    if(M>2){
      for(i in 3:M){
        Q <- Q%x%P
      }
    }
      # Combining rows and columns - Indices
    lr <- list()
    for(i in 1:M){
      lr[[i]] <- 0:1
    }
    r  <- expand.grid(lr)
    rt <- rowSums(r)
      # Combining columns
    Q0 <- Q[, rt==0]
    QM <- Q[, rt==M]
    Qx <- NULL
    for(i in 1:(M-1)){
      Qx <- cbind(Qx, rowSums(Q[, rt==i])/sum(rt==i))
    }
    Q <- cbind(Q0, Qx, QM)
      # Combining rows
    Q0 <- Q[rt==0, ]
    QM <- Q[rt==M, ]
    Qx <- NULL
    for(i in 1:(M-1)){
      Qx <- rbind(Qx, colSums(Q[rt==i, ]))
    }
    Q   <- rbind(Q0, Qx, QM) # final Sum Score transition matrix

    # RRT process for sum score variables
    ss.rr <- ss.p %*% t(Q)
    ss.y <- t(apply(ss.rr, 1, function(x) {rmultinom(1, 1, x)}))
    ss.y <- ss.y %*% c(0:M)
    yobs <- NULL
    for (i in 1:nrow(ss.y)) {
      yobs <- rbind(yobs, as.numeric(c(replicate(ss.y[i], 1),replicate(M-ss.y[i], 0))))
    }
    ss.df <- data.frame(yobs, x.var)
    names(ss.df) <- c(paste0("Y",pi.ordered$ix),paste0("x",1:length(x.betas)))
    return(ss.df)
  } else if (model=="IRT" | model=="irt") { # IRT RRT
    P <- matrix(c(1-p00, p00, p11, 1-p11), 2, 2)

    if (length(pi.1) < 2) {
      message("pi.1 must have two or more values for the irt model, adding a random value.")
      pi.1[2] <- runif(1)
    }
    pi.ordered <- sort(pi.1, decreasing=TRUE, index.return=TRUE)
    pi.1 <- pi.ordered$x

    theta <- x.var %*% x.betas
    gamma <- log(pi.1/(1-pi.1))
    irt.prev <- unlist(lapply(gamma, function(x) {lgp(theta+x)}))
    irt.prev <- matrix(irt.prev, nrow=n, byrow=F)
    irt.inds <- 1:ncol(irt.prev)
    irt.y <- NULL
    for (i in 1:length(irt.inds[c(TRUE,FALSE)])) {
      ind <- irt.inds[c(TRUE,FALSE)][i]
      irt.rr <- irt.prev[,c(ind:(ind+1))]%*%P
      irt.y <- cbind(irt.y, t(apply(irt.rr, 1, function(x) {rmultinom(1,1, x)}))[,2])
    }
    irt.df <- data.frame(irt.y, x.var)
    names(irt.df) <- c(paste0("y",pi.ordered$ix), paste0("x",1:length(x.betas)))
    return(irt.df)
  } else {
    stop("model must be one of the following: 'univariate', 'sumscore', or 'irt'.")
  }
}

#=============================================================================
# FUNCTION FOR OBTAINING BOOTSTRAP ESTIMATES
#=============================================================================
#' @title Bootstrap estimates for RRT models
#'
#' @description \code{RRbootstrap} generates bootstrap estimates for univariate (\code{\link{RRunivariate}}) and multivariate (\code{\link{RRsumscore}} and \code{\link{RRirt}}) RRT models, both with and without covariates. \code{RRbootstrap} is intended to allow users to deal with issues arising from non-invertible Hessians where standard errors cannot otherwise be computed. We caution, however, that when standard errors cannot be computed, that often gestures to issues with the underlying data (too few observations or traits that are too rare) and/or poorly-specified models.
#'
#' @param model A character string, one of \code{"univariate","sumscore",} or \code{"irt"}.
#' @param y Matrix or data frame of observed responses.
#' @param x Predictor variable values (matrix or data frame).
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no).
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes).
#' @param noncompliance Should the model specify question noncompliance (one-sided lying to avoid saying "yes")? \code{FALSE} (default) or \code{TRUE}. Noncompliance is only included in the \code{irt} model.
#' @param reps The number of samples with replacement that \code{RRbootstrap} takes to generate the estimates.
#' @param seed See \code{\link[base]{Random}}.
#'
#' @examples
#' ## BOOTSTRAP ESTIMATES FOR ITEM RESPONSE THEORY
#' # Step 1: Load or simulate multivariate RRT data
#' pi.1s <- c(0.15, 0.2, 0.3)
#' p00 <- p11 <- 5/6
#' IRTdata <- RRsimulate(model="irt", pi.1 = pi.1s, p00=p00, p11=p11, x.betas=c(-1, 3), n=200)
#'
#' # Step 2: Run RRbootstrap (NB: increase reps for actual analyses)
#' IRT.boot <- RRbootstrap("irt", y=IRTdata[,1:3], p00=p00, p11=p11, reps=50)
#' # Not run:
#' # IRT.boot <- RRbootstrap("irt", y=IRTdata[,1:3], x=IRTdata[,4:5], p00=p00, p11=p11, reps=500)
#' @return \code{RRbootstrap} returns a list \code{$bootresult} that provides the mean, 50\%, 2.5\%, and 97.5\% quantile bootstrap estimates for parameters generated by each RRT model.
#'
#' @seealso \code{\link{RRunivariate}} for single-item RRT estimation, \code{\link{RRsumscore}} for sum score modeling, and \code{\link{RRirt}} for item randomized response theory modeling. The latter two models (sum score and IRT) are suitable for multivariate RRT data (two or more question prompts).
#' @export
RRbootstrap <- function(model=c("univariate","sumscore","irt"), y, x=NULL, p00, p11, noncompliance = F, reps=500, seed=1){

  quants   <- function(B)quantile(B, probs=c(.5, .025, .975))

  set.seed(seed)

  ## Data prep
  if (is.null(dim(y))) {
    y <- data.matrix(y)
  }

  if(is.null(x)){
    yx <- y
  }else{
    if (is.null(dim(x))) {
      x <- as.data.frame(x)
    }

    if(is.null(colnames(x))){
      if (ncol(x) > 2) {
        colnames(x) <- paste0("X",1:ncol(x))
      } else {
        colnames(x) <- c("X1")
      }
      suppressMessages( blank.name("x") )
    }

    ## Dealing with numeric & non-numeric predictor data
    x        <- unclass(x)
    x        <- data.frame(x, stringsAsFactors = TRUE)
    x        <- model.matrix(~.+0, x)

    yx <- cbind(y, x)
  }

  if(model=="univariate"){

    if(is.null(x)){
      yx <- data.matrix(y)
    }else{

      if (is.null(dim(x))) {
        x <- as.data.frame(x)
      }

      if(is.null(colnames(x))){
        if (ncol(x) > 2) {
          colnames(x) <- paste0("X",1:ncol(x))
        } else {
          colnames(x) <- c("X1")
        }
        suppressMessages( blank.name("x") )
      }

      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~ ., x) # Include intercept

      yx <- cbind(y, x)
    }

    parnames <- rownames(RRunivariate(y=y, x=x, p00=p00, p11=p11, bootin=TRUE)$estimates)
    out <- matrix(0, reps, length(parnames))
    id <- 1:nrow(yx)

    ## Bootstrap: Univariate
    for(k in 1:reps){

      yx.boot <- data.matrix(yx[sample(id, replace=T), ])

      x.boot <- NULL

      if(!is.null(x)) {
        x.boot <- yx.boot[, -(1:ncol(y))]
        x.boot <- data.matrix(x.boot)
        if (!is.null(colnames(x))) {
          colnames(x.boot) <- colnames(x)
        } else {
          colnames(x.boot) <- paste0("X",1:ncol(x.boot))
        }
      }

      out[k, ] <- RRunivariate(y=yx.boot[, 1], x=x.boot, p00=p00, p11=p11, bootin=TRUE)$estimates$est

    }

  }else if(model=="sumscore"){

    if(is.null(x)){
      yx <- y
    }else{

      if (is.null(dim(x))) {
        x <- as.data.frame(x)
      }

      if(is.null(colnames(x))){
        if (ncol(x) > 2) {
          colnames(x) <- paste0("X",1:ncol(x))
        } else {
          colnames(x) <- c("X1")
        }
        suppressMessages(blank.name("x"))
      }
      ## Dealing with numeric & non-numeric predictor data
      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~.+0, x)

      yx <- cbind(y, x)
    }

    parnames <- rownames(RRsumscore(y=y, x=x, p00=p00, p11=p11, bootin=TRUE)$estimates)
    out <- matrix(0, reps, length(parnames))
    id <- nrow(yx)

    for(k in 1:reps){

      yx.boot <- data.matrix(yx[sample(id, replace=T), ])
      x.boot  <- NULL
      if(!is.null(x)) {
        x.boot <- yx.boot[, -(1:ncol(y))]
        x.boot <- data.matrix(x.boot)
        if (!is.null(colnames(x))) {
          colnames(x.boot) <- colnames(x)
        } else {
          colnames(x.boot) <- paste0("X",1:ncol(x.boot))
        }
      }

      out[k, ] <- RRsumscore(y=yx.boot[, 1:ncol(y)], x=x.boot, p00=p00, p11=p11, bootin=TRUE)$estimates$est

    }

  }else if(model=="irt" & !noncompliance){

    if(is.null(x)){
      yx <- y
    }else{
      if (is.null(dim(x))) {
        x <- as.data.frame(x)
      }

      if(is.null(colnames(x))){
        if (ncol(x) > 2) {
          colnames(x) <- paste0("X",1:ncol(x))
        } else {
          colnames(x) <- c("X1")
        }
        suppressMessages( blank.name("x") )
      }

      ## Dealing with numeric & non-numeric predictor data
      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~.+0, x)

      yx <- cbind(y, x)
    }

    parnames <- rownames(RRirt(y=y, x=x, p00=p00, p11=p11,bootin=TRUE)$estimates)
    out <- matrix(0, reps, length(parnames))
    id <- 1:nrow(yx)

    ## Bootstrap: IRT
    for(k in 1:reps){

      yx.boot <- data.matrix(yx[sample(id, replace=T), ])

      x.boot <- NULL

      if(!is.null(x)) {
        x.boot <- yx.boot[, -(1:ncol(y))]
        x.boot <- data.matrix(x.boot)
        if (!is.null(colnames(x))) {
          colnames(x.boot) <- colnames(x)
        } else {
          colnames(x.boot) <- paste0("X",1:ncol(x.boot))
        }
      }

      out[k, ] <- RRirt(y=yx.boot[, 1:ncol(y)], x=x.boot, p00=p00, p11=p11, bootin=TRUE)$estimates$est

    }

  }else if(model=="irt" & noncompliance){

    parnames <- rownames(RRirt(y=y, x=x, p00=p00, p11=p11, noncompliance=T, bootin=TRUE)$estimates)
    out <- matrix(0, reps, length(parnames))
    id <- 1:nrow(yx)

    ## IRT noncompliance bootstrap
    for(k in 1:reps){

      yx.boot <- data.matrix(yx[sample(id, replace=T), ])

      x.boot <- NULL

      if(!is.null(x)) {
        x.boot <- yx.boot[, -(1:ncol(y))]
        x.boot <- data.matrix(x.boot)
        if (!is.null(colnames(x))) {
          colnames(x.boot) <- colnames(x)
        } else {
          colnames(x.boot) <- paste0("X",1:ncol(x.boot))
        }
      }

      out[k, ] <- RRirt(y=yx.boot[, 1:ncol(y)], x=x.boot, p00=p00, p11=p11,
                         noncompliance=T, bootin=TRUE)$estimates$est

    }
  }

  return(list(bootresult = as.data.frame(cbind(mean=round(colMeans(out), 3),
                                               quan=round(t(apply(out, 2, quants)), 3)),
                                         row.names=parnames)))

}

#=============================================================================
# GOODNESS OF FIT CALCULATION
#=============================================================================

#' @title Goodness of Fit Calculation for RRT models
#'
#' @description \code{RRgof} calculates a parametric bootstrap \eqn{\chi^2} statistic for the goodness-of-fit (MacKenzie and Bailey 2004) of the three RRT models presented in this package: univariate (\code{\link{RRunivariate}}) and multivariate (\code{\link{RRsumscore}} and \code{\link{RRirt}}) RRT models, both with and without covariates. \code{RRgof} should be used to evaluate the model fit generated by one of the RRT models (\code{\link{RRunivariate}}, \code{\link{RRsumscore}}, or \code{\link{RRirt}}); it does not return parameter estimates for the models.
#'
#' @param model A character string, one of \code{"univariate","sumscore",} or \code{"irt"}.
#' @param y Matrix or data frame of observed responses.
#' @param x Predictor variable values (matrix or data frame).
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no).
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes).
#' @param noncompliance Should the model specify question noncompliance (one-sided lying to avoid saying "yes")? \code{FALSE} (default) or \code{TRUE}. Noncompliance is only included in the \code{irt} model (see \code{\link{RRirt}} for more information).
#' @param reps The number of parametric bootstraps performed to assess goodness of fit (default=500)
#' @param seed See \code{\link[base]{Random}}.
#'
#' @examples
#' ## GOODNESS OF FIT ESTIMATES
#' ### UNIVARIATE MODEL EXAMPLES
#' data(huntRRdata)
#' p00 <- p11 <- 5/6
#' reps <- 100
#' gof.univariate <- RRgof(model="univariate", y=huntRRdata['Bulbul'], p00=p00, p11=p11, reps=reps)
#' gof.univariate$gof.test.result
#' gof.univariate.GLM <- RRgof(model="univariate", y=huntRRdata['Bulbul'], x=huntRRdata['Age'], p00=p00, p11=p11, reps=reps)
#' gof.univariate.GLM$gof.test.result
#' ### IRT MODEL EXAMPLES (with noncompliance)
#' data(huntRRdata)
#' p00 <- p11 <- 5/6
#' reps <- 10 # Increase the number of reps for your own calculations (10 chosen for convenience)
#' gof.IRT.GLM <- RRgof(model="irt", y=huntRRdata[,c(1:3)], x=huntRRdata['Age'], noncompliance=TRUE, p00=p00, p11=p11, reps=reps)
#' @return \code{RRbootstrap} returns a list. \code{$gof.test.result} provides a data frame with the model string, \eqn{\chi^2_0} value (\code{chiSquare0}), degrees of freedom (\code{df}), p-value observed for original data against the model fit (\code{chiSquarePval}), and the parametric bootstrap p-value generated from resampling observations based on the original model (\code{boot.p.value}). \code{$X2boot} provides the bootstrapped \eqn{\chi^2_{boot}} values.
#'
#' @seealso \code{\link{RRunivariate}} for single-item RRT estimation, \code{\link{RRsumscore}} for sum score modeling, and \code{\link{RRirt}} for item randomized response theory modeling. The latter two models (sum score and IRT) are suitable for multivariate RRT data (two or more question prompts).
#' 
#' @references MacKenzie, Darryl I and Bailey, Larissa L. (2004).
#' \emph{Assessing the fit of site-occupancy models.}
#'  Journal of Agricultural, Biological, and Environmental Statistics. 9(3):300-318.
#' @export
RRgof <- function(model=c("univariate","sumscore","irt"), y, x=NULL, p00, p11, noncompliance = F, reps=500, seed=1){
  
  y  <- data.matrix(y)
  n  <- nrow(y)
  P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)
  set.seed(seed=seed, kind="L'Ecuyer-CMRG")
  
  if(model=="univariate"){
    
    if(is.null(x)){
      
      if(is.null(names(y))){
        parnames <-  "P(Y)"
      }else{
        parnames <- paste0("P(", colnames(y),")")
      }
      
      # Calculate reference X2
      result <- RRunivariate(y=y, p00=p00, p11=p11)$estimates 
      prev <- cbind(as.numeric(result[1]), 1-as.numeric(result[1]))
      pRR <- prev %*% P
      Mobs <- rbind(table(y), rev( round(pRR*nrow(y), 0) ) )
      chisq.0 <- suppressWarnings( chisq.test(Mobs) )
      
      # Generate bootstrap X2 observations
      y.obs <- NULL
      for (i in 1:reps) {
        y.intermediate <- rmultinom(nrow(y), 1, as.numeric(pRR))[1,]
        result.int <- RRunivariate(y=y.intermediate, p00=p00, p11=p11)$estimates
        prev.int <- cbind( as.numeric(result.int[1]), 1-as.numeric(result.int[1]))
        pRR.int <- prev.int %*% P
        y.obs <- cbind(y.obs, rmultinom(nrow(y), 1, as.numeric(pRR.int))[1,])
      }
      
      # Calculate X2 bootstrap
      chisq.vals <- c()
      for (i in 1:reps) {
        Mboot <- rbind(table(y), table(y.obs[,i]))
        rownames(Mboot) <- c("Expected","Observed")
        chi.test <- suppressMessages( chisq.test(Mboot) )
        chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
      }
    } else {
      if (is.null(dim(x))) {
        x <- as.data.frame(x)
      }
      
      if(is.null(colnames(x))){
        if (ncol(x) > 2) {
          colnames(x) <- paste0("X",1:ncol(x))
        } else {
          colnames(x) <- c("X1")
        }
        suppressMessages(blank.name("x"))
      }
      
      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~., x) # Include intercept
      
      # Calculate reference X2
      result <- RRunivariate(y=y, x=x, p00=p00, p11=p11, bootin = TRUE)$estimates$est
      prev <- cbind(lg(x%*%result), 1-lg(x%*%result))
      pRR <- prev %*% P
      yRR <- apply(pRR, 1, function(x) {rmultinom(1, 1, x)[1,]})
      Mobs <- rbind(table(y), table(yRR))
      chisq.0 <- suppressMessages( chisq.test(Mobs) )
      
      # Generate bootstrap X2 observations
      y.obs <- NULL
      for (i in 1:reps) {
        y.intermediate <- apply(pRR, 1, function(x) {rmultinom(1, 1, x)[1,]})
        result.int <- RRunivariate(y=y.intermediate, x=x, p00=p00, p11=p11, bootin = TRUE)$estimates$est
        prev.int <- cbind( lg(x%*%result.int), 1-lg(x%*%result.int))
        pRR.int <- prev.int %*% P
        y.obs <- cbind(y.obs, apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)[1,]}))
      }
      
      # Calculate X2 bootstrap
      chisq.vals <- c()
      for (i in 1:reps) {
        Mboot <- rbind(table(y), table(y.obs[,i]))
        rownames(Mboot) <- c("Expected","Observed")
        chi.test <- suppressMessages( chisq.test(Mboot) )
        chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
      }
    }
  } else if (model=="sumscore") {

    M  <- ncol(y)
    Q <- P%x%P
    if(M>2){
      for(i in 3:M){
        Q <- Q%x%P
      }
    }
    
    sy <- rowSums(y)
    ny <- matrix(0, M+1, 1)
    for(i in 0:M){
      ny[i+1] <- sum(sy==i)
    }
    py <- ny/n
    
    
    lr <- list()
    for(i in 1:M){
      lr[[i]] <- 0:1
    }
    r  <- expand.grid(lr)
    rt <- rowSums(r)
    
    Q0 <- Q[, rt==0]
    QM <- Q[, rt==M]
    Qx <- NULL
    for(i in 1:(M-1)){
      Qx <- cbind(Qx, rowSums(Q[, rt==i])/sum(rt==i))
    }
    
    Q <- cbind(Q0, Qx, QM)
    
    Q0 <- Q[rt==0, ]
    QM <- Q[rt==M, ]
    Qx <- NULL
    for(i in 1:(M-1)){
      Qx <- rbind(Qx, colSums(Q[rt==i, ]))
    }
    Q   <- rbind(Q0, Qx, QM)
    
    if(is.null(x)){
      # Calculate reference X2
      result <- RRsumscore(y=y, p00=p00, p11=p11, bootin = TRUE)$estimates$est
      pRR <- Q %*% result
      Mobs <- rbind( table(sy),round(as.numeric(pRR)*189,0) )
      rownames(Mobs) <- c("Expected","Observed")
      chisq.0 <- suppressMessages( chisq.test(Mobs) )
      
      # Generate bootstrap X2 observations
      chisq.vals <- c()
      for (i in 1:reps) {
        ss.int <- t(rmultinom(nrow(y), 1, as.numeric(pRR)))[,-1]
        ss.int <- ss.int %*% c(1:M)
        y.intermediate <- t( apply(ss.int, 1, function(x) {c(rep(1, x), rep(0, M-x))}) )
        result.int <- RRsumscore(y=y.intermediate, p00=p00, p11=p11)$estimates$est
        pRR.int <- Q %*% result.int
        
        Mboot <- rbind(table(sy), round(as.numeric(pRR.int)*189, 0))
        rownames(Mboot) <- c("Expected","Observed")
        chi.test <- suppressMessages( chisq.test(Mboot) )
        chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
      }
    } else if (!is.null(x)) {
        if(is.null(colnames(x))){
          if (ncol(x) > 2) {
            colnames(x) <- paste0("X",1:ncol(x))
          } else {
            colnames(x) <- c("X1")
          }
          suppressMessages(blank.name("x"))
        }
        ## Dealing with numeric & non-numeric predictor data
        x        <- unclass(x)
        x        <- data.frame(x, stringsAsFactors = TRUE)
        x        <- model.matrix(~.+0, x)
        
        # Calculate reference X2
        result <- RRsumscore(y=y, x=x, p00=p00, p11=p11, bootin = TRUE)$estimates$est
        alpha <- result[1:M]
        beta <- tail(result, -M)
        xb <- x%*%beta
        cp <- cbind( vapply( alpha, function(a) {lg(a-c(xb))}, rep(1,nrow(y))))
        pSS <- cbind(cp, 1) - cbind(0, cp)
        pRR <- pSS %*% t(Q)
        ySS <- t( apply(pRR, 1, function(x) {rmultinom(1, 1, x)}) )[,-1]
        ySS <- ySS %*% c(1:M)
        Mobs <- rbind( table(sy), table(ySS) )
        rownames(Mobs) <- c("Expected","Observed")
        chisq.0 <- suppressMessages( chisq.test(Mobs) )
        
        # Generate bootstrap X2 observations
        chisq.vals <- c()
        for (i in 1:reps) {
          ss.int <- t( apply(pRR, 1, function(x) {rmultinom(1, 1, x)}) )[,-1]
          ss.int <- ss.int %*% c(1:M)
          y.intermediate <- t( apply(ss.int, 1, function(x) {c(rep(1, x), rep(0, M-x))}) )
          result.int <- RRsumscore(y=y.intermediate, x=x, p00=p00, p11=p11, bootin=TRUE)$estimates$est
            # Calculating probabilties
          alpha.int <- result.int[1:M]
          beta.int <- tail(result.int, -M)
          xb.int <- x%*%beta.int
          cp.int <- cbind( vapply( alpha.int, function(a) {lg(a-c(xb.int))}, rep(1,nrow(y))))
          pSS.int <- cbind(cp.int, 1) - cbind(0, cp.int)
          pRR.int <- pSS.int %*% t(Q)
            # Generating observations
          ySS.int <- t( apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)}) )[,-1]
          ySS.int <- ySS.int %*% c(1:M)
            # Calculating X2-boot
          Mboot <- rbind(table(sy), table(ySS.int))
          rownames(Mboot) <- c("Expected","Observed")
          chi.test <- suppressMessages( chisq.test(Mboot) )
          chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
        }
      } 
  } else if (model=="irt") {
      M <- ncol(y) 
      if (is.null(x) & !noncompliance) {
        # Calculate reference X2
        result <- RRirt(y=y, p00=p00, p11=p11)$estimates$est
        prev <- cbind(result, 1-result)
        pRR <- prev %*% P
        pRR <- pRR[,c(2,1)]
        Mobs <- rbind( c( apply(y, 2, table)) , c(apply(pRR, 1, function(x) {round( x*nrow(y), 0)})) )
        chisq.0 <- suppressMessages( chisq.test(Mobs) )
        
        # Generate bootstrap X2 observations
        chisq.vals <- c()
        for (i in 1:reps) {
          y.intermediate <- apply(pRR, 1, function(x) {rmultinom(nrow(y), 1, x)[2,]})
          result.int <- RRirt(y=y.intermediate, p00=p00, p11=p11)$estimates$est
          prev.int <- cbind( result.int, 1-result.int)
          pRR.int <- prev.int %*% P
          pRR.int <- pRR.int[,c(2,1)]
          
            # Calculating bootstrap X2
          Mboot <- rbind( c(apply(y, 2, table)), c(apply(pRR.int, 1, function(x) {round( x*nrow(y), 0)})))
          rownames(Mboot) <- c("Expected","Observed")
          chi.test <- suppressMessages( chisq.test(Mboot) )
          chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
        }
      } else if (is.null(x) & noncompliance) {
        # Calculate reference X2
        result <- RRirt(y=y, p00=p00, p11=p11, noncompliance = TRUE)$estimates$est
        result <- head(result, M)
        prev <- cbind(result, 1-result)
        pRR <- prev %*% P
        pRR <- pRR[,c(2,1)]
        Mobs <- rbind( c( apply(y, 2, table)) , c(apply(pRR, 1, function(x) {round( x*nrow(y), 0)})) )
        chisq.0 <- suppressMessages( chisq.test(Mobs) )
        
        # Generate bootstrap X2 observations
        chisq.vals <- c()
        for (i in 1:reps) {
          y.intermediate <- apply(pRR, 1, function(x) {rmultinom(nrow(y), 1, x)[2,]})
          result.int <- RRirt(y=y.intermediate, p00=p00, p11=p11, noncompliance=TRUE)$estimates$est
          result.int <- head(result.int, M)
          prev.int <- cbind( result.int, 1-result.int)
          pRR.int <- prev.int %*% P
          pRR.int <- pRR.int[,c(2,1)]
          
          # Calculating bootstrap X2
          Mboot <- rbind( c(apply(y, 2, table)), c(apply(pRR.int, 1, function(x) {round( x*nrow(y), 0)})))
          rownames(Mboot) <- c("Expected","Observed")
          chi.test <- suppressMessages( chisq.test(Mboot) )
          chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
        }
      } else if (!is.null(x) & !noncompliance) {
        if (is.null(dim(x))) {
          x <- as.data.frame(x)
        }
        
        if(is.null(colnames(x))){
          if (ncol(x) > 2) {
            colnames(x) <- paste0("X",1:ncol(x))
          } else {
            colnames(x) <- c("X1")
          }
        }
        
        ## Dealing with numeric & non-numeric predictor data
        x        <- unclass(x)
        x        <- data.frame(x, stringsAsFactors = TRUE)
        x        <- model.matrix(~.+0, x)
        
        # Calculate reference X2
        result <- RRirt(y=y, x=x, p00=p00, p11=p11, bootin = TRUE)$estimates$est
        alpha <- head(result, M)
        beta <- tail(result, -M)
        xb <- x%*%beta
        yRR <- NULL
        for (i in 1:M) {
          prev <- cbind( c(lg(xb-alpha[i])), c(1-lg(xb-alpha[i])) )
          pRR <- prev%*%P
          yRR <- cbind(yRR, apply(pRR, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
        }
        Mobs <- rbind( c( apply(y, 2, table)) , c(apply(yRR, 2, table) ) )
        rownames(Mobs) <- c("Expected","Observed")
        chisq.0 <- suppressMessages( chisq.test(Mobs) )
        
        # Generate bootstrap X2 observations
        chisq.vals <- c()
        for (i in 1:reps) {
          y.intermediate <- NULL
          for (i in 1:M) {
            prev.int <- cbind( c(lg(xb-alpha[i])), c(1-lg(xb-alpha[i])) )
            pRR.int <- prev.int%*%P
            y.intermediate <- cbind(y.intermediate, apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
          }
            # Calculate IRT model
          result.int <- RRirt(y=y.intermediate, x=x, p00=p00, p11=p11, bootin=TRUE)$estimates$est
          alpha.int <- head(result.int, M)
          beta.int <- tail(result.int, -M)
          xb.int <- x%*%beta.int
          y.boot <- NULL
          for (i in 1:M) {
            prev.int <- cbind( c(lg(xb.int-alpha.int[i])), c(1-lg(xb.int-alpha.int[i])) )
            pRR.int <- prev.int%*%P
            y.boot <- cbind(y.boot, apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
          }
          Mboot <- rbind( c( apply(y, 2, table)) , c(apply(y.boot, 2, table) ) )
          rownames(Mboot) <- c("Expected","Observed")
          chi.test <- suppressMessages( chisq.test(Mboot) )
          chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
        }
      } else if (!is.null(x) & noncompliance) {
        if (is.null(dim(x))) {
          x <- as.data.frame(x)
        }
        
        if(is.null(colnames(x))){
          if (ncol(x) > 2) {
            colnames(x) <- paste0("X",1:ncol(x))
          } else {
            colnames(x) <- c("X1")
          }
        }
        
        ## Dealing with numeric & non-numeric predictor data
        x        <- unclass(x)
        x        <- data.frame(x, stringsAsFactors = TRUE)
        x        <- model.matrix(~.+0, x)
        
        # Calculate reference X2
        result <- RRirt(y=y, x=x, p00=p00, p11=p11, noncompliance = TRUE, bootin = TRUE)$estimates$est
        alpha <- head(result, M)
        beta.phi <- tail(result, -M)
        xb <- cbind(x,1)%*%beta.phi
        yRR <- NULL
        for (i in 1:M) {
          prev <- cbind( c(lg(xb-alpha[i])), c(1-lg(xb-alpha[i])) )
          pRR <- prev%*%P
          yRR <- cbind(yRR, apply(pRR, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
        }
        Mobs <- rbind( c( apply(y, 2, table)) , c(apply(yRR, 2, table) ) )
        rownames(Mobs) <- c("Expected","Observed")
        chisq.0 <- suppressMessages( chisq.test(Mobs) )
        
        # Generate bootstrap X2 observations
        chisq.vals <- c()
        for (i in 1:reps) {
          y.intermediate <- NULL
          for (i in 1:M) {
            prev.int <- cbind( c(lg(xb-alpha[i])), c(1-lg(xb-alpha[i])) )
            pRR.int <- prev.int%*%P
            y.intermediate <- cbind(y.intermediate, apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
          }
          # Calculate IRT model
          result.int <- RRirt(y=y.intermediate, x=x, p00=p00, p11=p11, noncompliance = TRUE, bootin = TRUE)$estimates$est
          alpha.int <- head(result.int, M)
          beta.phi.int <- tail(result.int, -M)
          xb.int <- cbind(x,1)%*%beta.phi.int
          y.boot <- NULL
          for (i in 1:M) {
            prev.int <- cbind( c(lg(xb.int-alpha.int[i])), c(1-lg(xb.int-alpha.int[i])) )
            pRR.int <- prev.int%*%P
            y.boot <- cbind(y.boot, apply(pRR.int, 1, function(x) {rmultinom(1, 1, x)[1,]}) )
          }
          Mboot <- rbind( c( apply(y, 2, table)) , c(apply(y.boot, 2, table) ) )
          rownames(Mboot) <- c("Expected","Observed")
          chi.test <- suppressMessages( chisq.test(Mboot) )
          chisq.vals <- append(chisq.vals, chi.test$statistic, after=length(chisq.vals))
        }
      }
    }
  
    # Output
    p.value <- sum(chisq.vals >= chisq.0$statistic)/reps
      # Based on: http://bit.ly/2CcX1Qz (\item{t0}; \item{t.star}) and http://bit.ly/2ojWId8 (Line 770)
    if (p.value==0) {
      p.display <- paste("<", 1/reps)
    } else {
      p.display = round(p.value, digits=4)
    }
    
    # Return output
    return(list(gof.test.result = data.frame(model=model, chiSquare0=chisq.0$statistic, df=chisq.0$parameter, chiSquarePval=chisq.0$p.value, boot.p.value=p.display),
                X2boot = chisq.vals))
}

#=============================================================================
# ADDITIONAL SUPPORT FUNCTIONS
#=============================================================================

# Logistic
lg <- function(x) {
  exp(x)/(1+exp(x))
}

# Paired logistic
lgp <- function(x) {
  c(lg(x),1-lg(x))
}

# Flipped paired logistic
lgpIRT <- function(z) {c(1-exp(z)/(1+exp(z)), exp(z)/(1+exp(z)))}

# Power analysis variance
var.pwr <- function(pi=NULL, P=NULL, n=NULL){
  pr  <- c(pi, 1-pi)
  pobs.RR <- P%*%pr
  return(pobs.RR[1]*(1-pobs.RR[1])/(n*(2*P[2,1]-1)^2))
}

# Hessian message
hess <- function() {message("Hessian not invertible")}

# Blank variable name message
blank.name<- function(x) {message(paste0("No names defined, assigning dummy variables for input ",x,"."))}

#=============================================================================
# DATA
#=============================================================================
#' @name huntRRdata
#' @aliases huntRRdata
#' @docType data
#' @title Randomized response technique (RRT) questionnaire on illegal bird hunting in Southwest China
#'
#' @description Individuals in Southwest China were asked whether or not they had hunted barbet (\emph{Psilopogon} spp.), bulbul (\emph{Pycnonotus} spp.), and/or partridge (\emph{Arborophila} spp) using the forced response RRT question design. This dataset presents one covariate: respondent age (centered and scaled).
#'
#' @format A data frame containing a sample of 189 observations.
#' Each question yielded a yes (\code{1}) or no (\code{0}) response. The RRT questions (response variables) were as follows:
#' \itemize{
#' \item \code{Barbet}: In the past year, have you hunted barbet?
#' \item \code{Bulbul}: In the past year, have you hunted bulbul?
#' \item \code{Partridge}: In the past year, have you hunted partridge?
#' }
#' In addition, the column \code{Age} presents respondent age as a predictor covariate (standardized to z-scores).
#'
#' @usage data(huntRRdata)
#'
#' @examples
#' ## Loading the data
#' data(huntRRdata)
#' p00 <- p11 <- 5/6
#'
#' ## Performing univariate RRT inference
#'    # How many individuals hunted barbets?
#' RRunivariate(y=huntRRdata['Barbet'], p00=p00, p11=p11)
#'    # Was age predictive of hunting barbets?
#' RRunivariate(y=huntRRdata['Barbet'], x=huntRRdata['Age'], p00=p00, p11=p11)
#'
#' ## Using RRbootstrap to determine confidence intervals for "Intercept" and Age
#'    # NB: Increase reps for more robust analyses
#' RRbootstrap(model="univariate", y=huntRRdata['Barbet'], x=huntRRdata['Age'],
#'             p00=p00, p11=p11, reps=50)
#'
#' ## Performing multivariate inference
#'     # Sum score model
#' RRsumscore(y=huntRRdata[,1:3], p00=p00, p11=p11)
#'     # Item response theory model
#' RRirt(y=huntRRdata[,1:3],p00=p00, p11=p11)
#'
#' @keywords datasets
#'
#' @references Chang, C.H., Williams, S.J., Cruyff, M.J.L.F., Zhang, M., Wilcove, D.S., Levin, S.A., and Quan, R.-C. (2017).
#' \emph{Predictors of illegal hunting in Southwest China.}
#' Unpublished Data.
#'
#' @seealso \code{RRdata} for RRT data on sensitive behaviors among university students.
"huntRRdata"

#' @name RRdata
#' @aliases RRdata
#' @docType data
#' @title Randomized Response Technique (RRT) Survey for a set of sensitive behaviors among university students
#'
#' @description  Undergraduate students were asked about a variety of sensitive topics including cheating, fighting, and bullying using the unrelated question (Horvitz) design.
#'
#' @format A data frame containing a sample of 710 observations.
#' Each question yielded a yes (\code{1}) or no (\code{0}) response. The RRT question topics are:
#' \itemize{
#'  \item copied: Have you ever copied in an exam?
#'  \item fought: Have you ever fought with a teacher?
#'  \item bullied: Have you been bullied?
#'  \item bullying: Have you ever bullied someone?
#'  \item drug: Have you ever taken drugs on campus?
#'  \item sex: Have you had sex on the premises of the university?
#' }
#'
#' @usage data(RRdata)
#'
#' @examples
#' ## Loading the data
#' data(RRdata)
#'
#' ## Calculating p00 and p11 for each question (see: ?RRunivariate)
#' college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
#' college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
#' college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
#'
#' ## Performing univariate RRT inference
#'    # How many students copied exams? (Question 1)
#' RRunivariate(RRdata['copied'], p00=college.p00[1], p11=college.p11[1])
#'    # How many students took drugs on campus? (Question 5)
#' RRunivariate(RRdata['drug'], p00=college.p00[5], p11=college.p11[5])
#'
#' ## Performing multivariate inference - note that the same values of p00 and p11 are required
#' question.indices <- c(1,6)
#'     # Sum score model
#' RRsumscore(RRdata[,question.indices], p00=college.p00[1], p11=college.p11[1])
#'     # Item response theory model
#' RRirt(RRdata[,question.indices],p00=college.p00[1], p11=college.p11[1])
#'
#' @keywords datasets
#'
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#'
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#'  Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#'
#' @references Cobo Rodriguez, B., del Mar Rueda Garcia, M., Arcos Cebrian, A. (2015).
#'  \emph{RRTCS: Randomized Response Techniques for Complex Surveys.}
#'  \url{https://CRAN.R-project.org/package=RRTCS}
#'
#' @seealso \code{\link[RRTCS]{HorvitzDataRealSurvey}} for the original RRdata source,
#' \code{huntRRdata} for RRT data on illegal hunting in Southwest China.
"RRdata"
