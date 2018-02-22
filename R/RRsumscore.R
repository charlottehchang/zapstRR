#' @title Sum score model for RRT data
#'
#' @description \code{RRsumscore} estimates (1) the prevalence of 0..K sensitive traits, (2) ordinal regression when a matrix or data frame of predictor variables is provided.
#'
#' @param y Matrix or data frame of observed responses for K RRT items
#' @param x Predictor variable values (matrix or data frame).
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)
#' @param n.round Number of digits for rounding estimates (default is \code{n.round=3}).
#' @param bootin Default: \code{bootin=FALSE} - do not change. For RRbootstrap internal functioning.
#'
#' @examples
#' ## Sum score modeling with simulated data (see ?RRsimulate for questions)
#' p00 <- p11 <- 5/6
#' SSdata <- RRsimulate(model="sumscore",pi.1 = c(0.15, 0.2, 0.4),
#'                      p00=p00, p11=p11, x.betas=c(2), n=375)
#' head(SSdata)
#' # Sum score prevalence
#' RRsumscore(y=SSdata[,1:3], p00=p00, p11=p11, n.round=2)
#' # Ordinal regression modeling
#' RRsumscore(y=SSdata[,1:3], x=SSdata[,4], p00=p00, p11=p11)
#'
#' ## Sum score example: the number of species hunted in Southwest China
#'     # Load the data and prepare the design parameters
#' data(huntRRdata) # See ?huntRRdata for information on this dataset
#' p00 <- p11 <- 5/6
#'     # Sum score prevalence
#' RRsumscore(y=huntRRdata[,1:3], p00=p00, p11=p11)
#'     # Ordinal regression: the effect of age on the number of species hunted
#' RRsumscore(y=huntRRdata[,1:3], x=huntRRdata['Age'], p00=p00, p11=p11)
#'
#' ## Sum score example: University student behaviors
#'     # Loading the data
#' data(RRdata)  # See ?RRdata for more information on the sensitive questionnaire.
#'
#'     # Calculating p00 and p11 for each question
#' college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
#' college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
#' college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
#' question.indices <- c(1,6) # Note that the p00 and p11 values MUST be the same
#'                            # for all questions used in IRT or sumscore analyses
#'
#'     # Performing sum score modeling
#' RRsumscore(RRdata[,question.indices], p00=college.p00[1], p11=college.p11[1])
#'
#' @return When \code{x=NULL} (default setting), \code{RRsumscore} yields the prevalence of RRT sum scores; otherwise, RRsumscore performs ordinal regression. In both cases, \code{RRsumscore} returns a list:
#' \itemize{
#'    \item \code{$estimates}: a data frame with the parameter estimate (\code{est}), standard error (\code{se}), T-statistic (\code{tval}), lower bounds of the 95\% confidence interval (\code{min95}), and upper bounds of the 95\% confidence interval (\code{max95}).
#'    \item \code{$fitmeasures}: -2*log likelihood (\code{minus2LL}), degrees of freedom (df), Akaike Information Criterion (\code{AIC}), and AIC corrected for (small) sample sizes (\code{AICc}).
#' }
#' When the Hessian is singular, standard errors are inestimable, as are the t-statistic, p-value, and confidence limits; in this case, the model will return \code{NaN} for those entries.
#' @seealso \code{\link{RRbootstrap}} for bootstrap estimates when the Hessian is singular.
#' @export
RRsumscore <- function(y, x=NULL, p00, p11, n.round=3,bootin=FALSE){

  ## Data cleaning/preparation
  if (is.null(dim(y))) {
    stop("The input y must contain at least two or more question items.")
  } else if (ncol(y) < 2) {
    stop("The input y must contain at least two or more question items.")
  }

  y.check <- apply(y, 2, function(y) {!all(y %in% 0:1)})
  if(sum(y.check) > 0) {
    stop("The input y must contain responses coded as 0 or 1 (not possessing or possessing the sensitive trait).")
  }

  if(is.null(names(y))){
    names.y <- paste0("Y",1:ncol(y))
    suppressMessages(blank.name("y"))
  }else{
    names.y <- names(y)
  }

  pi.check <- findInterval(c(p00, p11), c(0,1))
  pi.check <- ifelse(pi.check != 1, 0, 1)
  if ( sum(pi.check) != 2) {
    stop("Both p00 and p11 must be within the interval (0,1).")
  }

  y  <- data.matrix(y)
  n  <- nrow(y)
  M  <- ncol(y)

  ## Generating transition matrix
  P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)
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

    npar <- M

    parnames <-  paste0("P(sumscore=", 0:ncol(y),")")

    logl <- function(b){

      a <- NULL
      for(j in 1:M){
        a <- c(a, sum(exp(b[1:j])))
      }

      cp  <- matrix(lg(a), M, 1)
      p   <- rbind(cp, 1) - rbind(0, cp)
      -t(ny)%*%log(Q%*%p)
    }

    fit <- optim(rep(-2, npar), logl, control=list(maxit=1000))

    alpha <- NULL
    for(j in 1:M){
      alpha <- c(alpha, sum(exp(fit$par[1:j])))
    }
    cum.pr <- matrix(lg(alpha), M, 1)
    est    <- rbind(cum.pr, 1) - rbind(0, cum.pr)

    logl.p <- function(b){
      p <- c(1-sum(b), b)
      -t(ny)%*%log(Q%*%p)
    }

    vcv <- tryCatch(solve(optim(est[-1], logl.p, hessian=T, control=list(maxit=1))$hessian),
                    error   = function(e)"singular",
                    warning = function(w)"infnite")

    if(class(vcv)=="matrix"){

      se    <- sqrt(c(sum(vcv), diag(vcv)))
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), n-length(est)-1)
      min95 <- ifelse(est - 1.96*se < 0, 0, est - 1.96*se)
      max95 <- ifelse(est + 1.96*se > 1, 1, est + 1.96*se)

    }else{

      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN

    }


  }else if(!is.null(x)){

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
        suppressMessages(blank.name("x"))
      }
      ## Dealing with numeric & non-numeric predictor data
      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~.+0, x)
    }

    parnames <- c(paste0("alpha", 0:(M-1)), colnames(x))
    npar <- M+ncol(x)

    logl <- function(b){

      alpha <- NULL
      for(j in 1:M){
        alpha <- c(alpha, sum(exp(b[1:j])))
      }
      alpha     <- matrix(alpha, n, M, byrow=T)
      beta      <- b[(M+1):length(b)]
      xb        <- x%*%beta
      cum.probs <- lg(alpha  - c(xb))
      probs     <- cbind(cum.probs, 1) - cbind(0, cum.probs)

      logl  <- 0
      for(j in 0:M){
        logl <- logl - sum(log(Q[j+1, ]%*%t(probs[sy==j, ])))
      }
      return(logl)
    }

    fit <- optim(c(rep(-2, M), rep(0, ncol(x))), logl, control=list(maxit=5000))

    alpha <- NULL
    for(j in 1:M){
      alpha     <- c(alpha, sum(exp(fit$par[1:j])))
    }

    est      <- c(alpha, fit$par[-(1:M)])

    logl.a <- function(b){

      a     <- matrix(b[1:M], n, M, byrow=T)
      beta  <- b[(M+1):length(b)]
      xb    <- x%*%beta
      cp    <- lg(a  - c(xb))
      p      <- cbind(cp, 1) - cbind(0, cp)

      logl  <- 0
      for(j in 0:M){
        logl <- logl - sum(log(Q[j+1, ]%*%t(p[sy==j, ])))
      }
      return(logl)
    }

    vcv <- tryCatch(solve(optim(est, logl.a, hessian=T, control=list(maxit=1))$hessian),
                    error   = function(e)"singular",
                    warning = function(w)"infinite")

    if(class(vcv)=="matrix"){

      se    <- suppressWarnings( sqrt(diag(vcv)) )
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), n-length(est)-1)
      min95 <- est - 1.96*se
      max95 <- est + 1.96*se

    }else{
      suppressMessages(hess())
      se <- tval <- pval <- min95 <- max95 <- NaN
    }
  }

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
