#' @title Item response theory model for RRT data
#'
#' @description \code{RRirt} estimates (1) the prevalence of a set of sensitive traits, and (2) performs logistic regression for multivariate RRT data (2 or more sensitive traits) when predictor variables are provided, and can specify question noncompliance (one-sided lying to avoid admitting yes).
#'
#' @param y Matrix or data frame of observed responses for two or more RRT items
#' @param x Predictor variable values (matrix or data frame).
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)
#' @param noncompliance Should the model specify question noncompliance (one-sided lying to avoid saying "yes")? \code{FALSE} (default) or \code{TRUE}.
#' @param n.round Number of digits for rounding estimates (default is \code{n.round=3}).
#' @param bootin Default: \code{bootin=FALSE} - do not change. For RRbootstrap internal functioning.
#'
#' @examples
#' ## Item response theory estimation
#' # Load or simulate multivariate RRT data (see ?RRsimulate for questions)
#' pi.1s <- c(0.45, 0.2, 0.1)
#' p00 <- p11 <- 5/6
#' IRTdata <- RRsimulate(pi.1 = pi.1s, p00=p00, p11=p11,
#'                       x.betas=c(1, 3), n=200, model="irt")
#'
#' # IRT prevalence estimation
#' RRirt(y=IRTdata[,1:3], p00=5/6, p11=5/6)
#'
#' # IRT logistic regression
#' RRirt(y=IRTdata[,1:3], x=IRTdata[,4:5], p00=5/6, p11=5/6)
#'
#' ## Item response theory example: hunting in Southwest China
#'     # Load the data and prepare the design parameters
#' data(huntRRdata) # See ?huntRRdata for information on this dataset
#' p00 <- p11 <- 5/6
#'     # IRT prevalence analysis
#' RRirt(y=huntRRdata[,c(1:3)], p00=p00, p11=p11)
#'     # IRT prevalence corrected for phi, or the rate of evasive response bias (noncompliance=TRUE)
#' RRirt(y=huntRRdata[,c(1:3)], p00=p00, p11=p11, noncompliance=TRUE)
#'     # Multivariate regression
#' RRirt(y=huntRRdata[,c(3,1,2)], x=huntRRdata['Age'], p00=p00, p11=p11)
#'
#' ## Item response theory example: university student behavior
#'     # Loading the data
#' data(RRdata) # See ?RRdata for more information on the sensitive questionnaire.
#'
#'     # Calculating p00 and p11 for each question
#' college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
#' college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
#' college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
#' question.indices <- c(1,6) # Note that the p00 and p11 values MUST be the same
#'                            # for all questions used in IRT or sumscore analyses
#'
#'     # Performing IRT call
#' RRirt(RRdata[,question.indices],p00=college.p00[1], p11=college.p11[1])
#'
#' @return When \code{x=NULL} (default setting), \code{RRirt} yields the prevalence of sensitive traits; otherwise, RRirt performs item randomized response regression. In both cases, \code{RRirt} returns a list:
#' \itemize{
#'    \item \code{$estimates}: a data frame with the parameter estimate (\code{est}), standard error (\code{se}), T-statistic (\code{tval}), lower bounds of the 95\% confidence interval (\code{min95}), and upper bounds of the 95\% confidence interval (\code{max95}).
#'    \item \code{$fitmeasures}: -2*log likelihood (\code{minus2LL}), degrees of freedom (df), Akaike Information Criterion (\code{AIC}), and AIC corrected for (small) sample sizes (\code{AICc}).
#' }
#' When the Hessian is singular, standard errors are inestimable, as are the t-statistic, p-value, and confidence limits; in this case, the model will return \code{NaN} for those entries.
#' @seealso \code{\link{RRbootstrap}} for bootstrap estimates when the Hessian is singular.
#' @export
RRirt <- function(y, x=NULL, p00, p11, noncompliance=F, n.round=3, bootin=FALSE){

  ## Checking input
  pi.check <- findInterval(c(p00,p11), c(0,1))
  pi.check <- ifelse(pi.check != 1, 0, 1)
  if ( sum(pi.check) != 2 ) {
    stop("All pi.1, p00, and p11 values must be within the interval (0,1).")
  }

  if (is.null(dim(y))) {
    stop("The input y must contain at least two or more question items.")
  } else if (ncol(y) < 2) {
    stop("The input y must contain at least two or more question items.")
  }

  y.check <- apply(y, 2, function(y) {!all(y %in% 0:1)})
  if(sum(y.check) > 0) {
    stop("The input y must contain responses coded as 0 or 1 (not possessing or possessing the sensitive trait).")
  }

  ## Preparing data
  y <- data.matrix(y)
  if(is.null(colnames(y))){
    names.y <- paste0("Y", 1:ncol(y))
    suppressMessages(blank.name("y"))
  }else{
    names.y <- colnames(y)
  }

  if(!is.null(x)){

    if (!isTRUE(bootin)) {

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
    }
    names.x <- colnames(x)
  }

  n  <- nrow(y)
  M  <- ncol(y)
  P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)
  Q  <- P%x%P
  if(M > 2){
    for(i in 3:M){
      Q     <- Q%x%P
    }
  }

  dt <- as.data.frame(table(as.data.frame(y)))

  lr    <- list()
  for(i in 1:M){lr[[i]] <- 0:1}
  r     <- rev(expand.grid(lr))
  for(j in 1:M){dt <- dt[order(r[, j]), ]}
  nobs  <- dt[, ncol(dt)]

  pr.r    <- apply(r, 1, function(x) paste(x, collapse=""))
  pr.y    <- apply(y, 1, function(x) paste(x, collapse=""))
  p.id    <- match(pr.y, pr.r)
  profile <- sort(unique(p.id))


  if(is.null(x) & !noncompliance){

    parnames <- paste0("P(", names.y,")")
    npar     <- length(parnames)

    logl <- function(b){
      pmat  <- matrix(lgpIRT(-b), M, 2)
      Probs <- pmat[1, ]%x%pmat[2, ]
      if(M > 2){
        for(i in 3:M){
          Probs <- Probs%x%pmat[i, ]
        }
      }
      logl  <- t(nobs)%*%log(Q%*%Probs)

      return(-logl)
    }

    fit <- optim(rep(0, M), logl, control=list(maxit=1000))

    est      <- lg(-fit$par)
    minus2LL <- 2*fit$value
    npar     <- length(fit$par)
    AIC      <- minus2LL + 2*npar
    AICc     <- AIC + 2*npar*(npar+1)/(n-npar-1)

    logl.p <- function(b){
      pmat   <- cbind(1-b, b)
      Probs <- pmat[1, ]%x%pmat[2, ]
      if(M > 2){
        for(i in 3:M){
          Probs <- Probs%x%pmat[i, ]
        }
      }

      return(-t(nobs)%*%log(Q%*%Probs))
    }

    vcv <- tryCatch(solve(optim(est, logl.p, hessian=T, control=list(maxit=1))$hessian),
                    error = function(e) return("singular"))

    if(class(vcv)=="matrix"){

      se    <- suppressWarnings( sqrt(diag(vcv)) )
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), n-npar)
      min95 <- ifelse(est-1.96*se < 0, 0, est-1.96*se)
      max95 <- ifelse(est+1.96*se > 1, 1, est+1.96*se)

    }else{

      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN

    }


  }else if(is.null(x) & noncompliance){

    parnames <- c(paste0("P(", names.y,")"), "phi")
    npar     <- length(parnames)


    logl <- function(b){
      pmat  <- matrix(lgpIRT(-b[-length(b)]), M, 2)
      Probs <- pmat[1, ]%x%pmat[2, ]
      if(M > 2){
        for(i in 3:M){
          Probs <- Probs%x%pmat[i, ]
        }
      }
      sp     <- lg(b[length(b)])
      R      <- (1-sp)*Q
      R[1, ] <- R[1, ] + sp

      logl  <- t(nobs)%*%log(R%*%Probs)

      return(-logl)
    }

    fit <- optim(c(rep(0, M), -3), logl, control=list(maxit=1000))

    est      <- c(lg(-fit$par[1:M]), lg(fit$par[M+1]))
    minus2LL <- 2*fit$value
    npar     <- length(fit$par)
    AIC      <- minus2LL + 2*(npar)
    AICc     <- AIC + 2*npar*(npar+1)/(n-npar-1)

    logl.p <- function(b){
      pmat   <- cbind(1-b, b)
      Probs <- pmat[1, ]%x%pmat[2, ]
      if(M > 2){
        for(i in 3:M){
          Probs <- Probs%x%pmat[i, ]
        }
      }
      sp     <- b[length(b)]
      R      <- (1-sp)*Q
      R[1, ] <- R[1, ] + sp

      return(-t(nobs)%*%log(R%*%Probs))
    }

    vcv <- tryCatch(solve(optim(est, logl.p, hessian=T, control=list(maxit=1))$hessian),
                    error   = function(e) return("singular"),
                    warning = function(w) return("infinite"))


    if(class(vcv)=="matrix"){

      se    <- suppressWarnings( sqrt(diag(vcv)) )
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), n-npar)
      min95 <- ifelse(est - 1.96*se < 0, 0, est - 1.96*se)
      max95 <- ifelse(est + 1.96*se > 1, 1, est + 1.96*se)

    }else{

      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN

    }


  }else if(!is.null(x) & !noncompliance){

    parnames <- c(paste("gamma", names.y, sep="."), names.x)
    npar     <- length(parnames)

    logl <- function(b){

      alpha <- b[1:M]
      beta  <- b[(M+1):length(b)]
      xb    <- x%*%beta

      p2 <- 2^(0:M)
      pv <- 1
      for(j in 1:M){
        Pr <- matrix(lgpIRT(xb-alpha[j]), ncol=2)
        pv  <- pv*Pr[, rep(1:2, times=p2[j], each=p2[M+1-j])]
      }

      logl    <- 0
      for(j in 1:length(profile)){
        logl <- logl + sum(log(Q[j, ]%*%t(pv[p.id==j, ])))
      }

      return(-sum(logl))
    }

    fit <- optim(rep(0, npar), logl)

    est      <- fit$par
    minus2LL <- 2*fit$value
    npar     <- length(fit$par)
    AIC      <- minus2LL + 2*(npar+1)
    AICc     <- AIC + 2*npar*(npar+1)/(n-npar-1)

    vcv <-  tryCatch(solve(optimHess(est, logl)),
                     error   = function(e) return("error"),
                     warning = function(w) return("infinite"))

    if(class(vcv)=="matrix"){

      se    <- suppressWarnings( sqrt(diag(vcv)) )
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), df=n-length(est)-1)
      min95 <- est-1.96*se
      max95 <- est+1.96*se

    }else{

      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN

    }


  }else if(!is.null(x) & noncompliance){

    parnames <- c(paste("gamma", names.y, sep="."), names.x, "phi")
    npar     <- length(parnames)

    logl <- function(b){

      alpha <- b[1:M]
      beta  <- b[(M+1):(length(b)-1)]
      phi   <- lg(b[length(b)])
      R     <- (1-phi)*Q
      R[1,] <- R[1, ] + phi
      x     <- as.matrix(x)
      xb    <- x%*%beta

      p2 <- 2^(0:M)
      pv <- 1
      for(j in 1:M){
        Pr <- matrix(lgpIRT(xb-alpha[j]), ncol=2)
        pv  <- pv*Pr[, rep(1:2, times=p2[j], each=p2[M+1-j])]
      }

      profile <- sort(unique(p.id))
      logl    <- 0
      for(j in 1:length(profile)){
        logl <- logl + sum(log(R[j, ]%*%t(pv[p.id==j, ])))
      }

      return(-sum(logl))
    }

    fit <- optim(c(rep(0, npar-1), -3), logl, control=list(maxit=1000))

    est       <- fit$par
    est[npar] <- lg(fit$par[npar])

    minus2LL <- 2*fit$value
    npar     <- length(fit$par)
    AIC      <- minus2LL + 2*(npar+1)
    AICc     <- AIC + 2*npar*(npar+1)/(n-npar-1)

    vcv <-  tryCatch(solve(optimHess(fit$par, logl)),
                     error   = function(e) return("singular"),
                     warning = function(w) return("infinite"))

    if(class(vcv)=="matrix"){

      se       <- suppressWarnings( sqrt(diag(vcv)) )
      se[npar] <- (exp(fit$par[npar])/(1+exp(fit$par[npar]))^2)*se[npar]
      tval <- est/se
      pval <- 2*pt(-abs(tval), df=n-npar)
      min95 <- est-1.96*se
      max95 <- est+1.96*se

    }else{

      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN

    }

  }

  # Global return call
  return(list(
    estimates   = round(data.frame(est       = est,
                                   se        = se,
                                   tval      = tval,
                                   pval      = pval,
                                   min95     = min95,
                                   max95     = max95,
                                   row.names = parnames), n.round),
    fitmeasures = round(data.frame(minus2LL  = 2*fit$value,
                                   df        = n - npar,
                                   AIC       = 2*(fit$value + npar),
                                   AICc      = 2*(fit$value  + 2*npar*(npar+1)/(n-npar-1)),
                                   row.names = ""), n.round)))
}
