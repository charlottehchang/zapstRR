% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRirt.R
\name{RRirt}
\alias{RRirt}
\title{Item response theory model for RRT data}
\usage{
RRirt(y, x = NULL, p00, p11, noncompliance = F, n.round = 3,
  bootin = FALSE)
}
\arguments{
\item{y}{Matrix or data frame of observed responses for two or more RRT items}

\item{x}{Predictor variable values (matrix or data frame).}

\item{p00}{The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)}

\item{p11}{The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)}

\item{noncompliance}{Should the model specify question noncompliance (one-sided lying to avoid saying "yes")? \code{FALSE} (default) or \code{TRUE}.}

\item{n.round}{Number of digits for rounding estimates (default is \code{n.round=3}).}

\item{bootin}{Default: \code{bootin=FALSE} - do not change. For RRbootstrap internal functioning.}
}
\value{
When \code{x=NULL} (default setting), \code{RRirt} yields the prevalence of sensitive traits; otherwise, RRirt performs item randomized response regression. In both cases, \code{RRirt} returns a list:
\itemize{
   \item \code{$estimates}: a data frame with the parameter estimate (\code{est}), standard error (\code{se}), T-statistic (\code{tval}), lower bounds of the 95\% confidence interval (\code{min95}), and upper bounds of the 95\% confidence interval (\code{max95}).
   \item \code{$fitmeasures}: -2*log likelihood (\code{minus2LL}), degrees of freedom (df), Akaike Information Criterion (\code{AIC}), and AIC corrected for (small) sample sizes (\code{AICc}).
}
When the Hessian is singular, standard errors are inestimable, as are the t-statistic, p-value, and confidence limits; in this case, the model will return \code{NaN} for those entries.
}
\description{
\code{RRirt} estimates (1) the prevalence of a set of sensitive traits, and (2) performs logistic regression for multivariate RRT data (2 or more sensitive traits) when predictor variables are provided, and can specify question noncompliance (one-sided lying to avoid admitting yes).
}
\examples{
## Item response theory estimation
# Load or simulate multivariate RRT data (see ?RRsimulate for questions)
pi.1s <- c(0.45, 0.2, 0.1)
p00 <- p11 <- 5/6
IRTdata <- RRsimulate(pi.1 = pi.1s, p00=p00, p11=p11,
                      x.betas=c(1, 3), n=200, model="irt")

# IRT prevalence estimation
RRirt(y=IRTdata[,1:3], p00=5/6, p11=5/6)

# IRT logistic regression
RRirt(y=IRTdata[,1:3], x=IRTdata[,4:5], p00=5/6, p11=5/6)

## Item response theory example: hunting in Southwest China
    # Load the data and prepare the design parameters
data(huntRRdata) # See ?huntRRdata for information on this dataset
p00 <- p11 <- 5/6
    # IRT prevalence analysis
RRirt(y=huntRRdata[,c(1:3)], p00=p00, p11=p11)
    # IRT prevalence corrected for phi, or the rate of evasive response bias (noncompliance=TRUE)
RRirt(y=huntRRdata[,c(1:3)], p00=p00, p11=p11, noncompliance=TRUE)
    # Multivariate regression
RRirt(y=huntRRdata[,c(3,1,2)], x=huntRRdata['Age'], p00=p00, p11=p11)

## Item response theory example: university student behavior
    # Loading the data
data(RRdata) # See ?RRdata for more information on the sensitive questionnaire.

    # Calculating p00 and p11 for each question
college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
question.indices <- c(1,6) # Note that the p00 and p11 values MUST be the same
                           # for all questions used in IRT or sumscore analyses

    # Performing IRT call
RRirt(RRdata[,question.indices],p00=college.p00[1], p11=college.p11[1])

}
\seealso{
\code{\link{RRbootstrap}} for bootstrap estimates when the Hessian is singular.
}
