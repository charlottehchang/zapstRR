#' @title Single-item (univariate) RRT modeling for prevalence and logistic regression
#'
#' @description \code{RRunivariate} finds (1) the prevalence of a sensitive trait given a specific RRT question design, and (2) can perform logistic regression for individual RRT items.
#'
#' @param y Observed responses (vector or 1-column matrix/data frame)
#' @param x Predictor variable values (data frame). An intercept is automatically added.
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)
#' @param n.round Number of digits for rounding estimates (default is \code{n.round=3})
#' @param bootin Default: \code{bootin=FALSE} - do not change. For RRbootstrap internal functioning.
#'
#' @examples
#' ## Univariate analysis of simulated data (see ?RRsimulate for questions)
#' RRTdata <- RRsimulate(pi.1=0.3, p00=5/6, p11=5/6, n=300,
#'                       model="univariate", x.betas=c(-1,3))
#' head(RRTdata) # see first 6 rows of simulated data
#'    # Univariate estimation of sensitive trait prevalence
#' RRunivariate(RRTdata[,1], p00=5/6, p11=5/6)
#'    # Univariate logistic regression
#' RRunivariate(y=RRTdata[,1], x=RRTdata[,2:3], p00=5/6, p11=5/6)
#'
#' ## Univariate analysis of bulbul hunting in Southwest China
#'     # Load the data and prepare the design parameters
#' data(huntRRdata) # See ?huntRRdata for information on this dataset
#' p00 <- p11 <- 5/6
#'     # Prevalence of hunting bulbuls
#' RRunivariate(y=huntRRdata['Bulbul'], p00=p00, p11=p11)
#'     # Effect of age on the probability of hunting bulbuls
#' RRunivariate(y=huntRRdata['Bulbul'], x=huntRRdata['Age'],p00=p00, p11=p11)
#'
#' ## Univariate analysis of sensitive behaviors among university students
#'    # Loading the data
#' data(RRdata) # See ?RRdata for more information on the sensitive questionnaire.
#'
#'    # Calculating p00 and p11 for each question
#' college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
#' college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
#' college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
#'
#'    # Performing univariate RRT inference
#'      # How many students copied exams? (Question 1)
#' RRunivariate(y=RRdata['copied'], p00=college.p00[1], p11=college.p11[1])
#'      # How many students took drugs on campus? (Question 5)
#' RRunivariate(RRdata['drug'], p00=college.p00[5], p11=college.p11[5])
#'
#' @return \code{RRunivariate} returns a list:
#' \itemize{
#'    \item \code{$estimates}: a data frame with the parameter estimate (\code{est}), standard error (\code{se}), T-statistic (\code{tval}), lower bounds of the 95\% confidence interval (\code{min95}), and upper bounds of the 95\% confidence interval (\code{max95}).
#'    \item \code{$fitmeasures}: -2*log likelihood (\code{minus2LL}), degrees of freedom (df), Akaike Information Criterion (\code{AIC}), and AIC corrected for (small) sample sizes (\code{AICc}).
#' }
#' When the Hessian is singular, standard errors are inestimable, as are the t-statistic, p-value, and confidence limits; in this case, the model will return \code{NaN} for those entries.
#' @seealso \code{\link{RRbootstrap}} for bootstrap estimates when the Hessian is singular.
#' @export
RRunivariate <- function(y, x=NULL, p00, p11, n.round=3, bootin=FALSE){

  if (!is.data.frame(y) & !is.matrix(y)) {
    if(!all(y %in% 0:1)) {
      stop("The input y must contain responses coded as 0 or 1 (not possessing or possessing the sensitive trait).")
    }
  } else {
    if(ncol(y) > 1){
      stop("The univariate RRT model can only take single items; y must only have one column of data.")
    }
  }

  if(is.null(names(y))){
    names.y <- "Y"
    suppressMessages(blank.name("y"))
  }else{
    names.y <- names(y)
  }

  y  <- data.matrix(y)
  n  <- nrow(y)
  P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)

  if(is.null(x)){

    if(is.null(names(y))){
      parnames <-  "P(Y)"
    }else{
      parnames <- paste0("P(", colnames(y),")")
    }

    npar     <- length(parnames)

    pobs  <- c(sum(y==0), sum(y==1))/n
    est   <- (solve(P)%*%pobs)[2]
    se    <- sqrt(pobs[1]*pobs[2]/(n*(P[2,2]-P[2,1])^2))
    tval  <- est/se
    pval  <- 2*pt(-abs(tval), n-length(est)-1)
    min95 <- ifelse(est - 1.96*se < 0, 0, est - 1.96*se)
    max95 <- ifelse(est + 1.96*se > 1, 1, est + 1.96*se)

    fit   <-  data.frame(value = -n*t(pobs)%*%log(P%*%c(1-est, est)))


  } else if(!is.null(x)){

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

      x        <- unclass(x)
      x        <- data.frame(x, stringsAsFactors = TRUE)
      x        <- model.matrix(~., x) # Include intercept
    }

    parnames <-  colnames(x)

    npar     <- length(parnames)

    logl <- function(b){
      p     <- t(cbind(1-lg(x%*%b), lg(x%*%b)))
      logl  <- log((P[1, ] %*% p)^(1-c(y)))+log((P[2, ] %*% p)^c(y))
      return(-sum(logl))
    }

    fit <- optim(rep(0, npar), logl)

    vcv <-  tryCatch(solve(optimHess(fit$par, logl)),
                     error   = function(e){return("singular")},
                     warning = function(e){return("infinite")})

    if(class(vcv)!="matrix"){

      suppressMessages(hess())
      est <- round(fit$par, 3)
      se <- tval <- pval <- min95 <- max95 <- NaN

    }else{

      est   <- fit$par
      se    <- suppressWarnings( sqrt(diag(vcv)) )
      tval  <- est/se
      pval  <- 2*pt(-abs(tval), n-length(est)-1)
      min95 <- est - 1.96*se
      max95 <- est + 1.96*se

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
