% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRsumscore.R
\name{RRsumscore}
\alias{RRsumscore}
\title{Sum score model for RRT data}
\usage{
RRsumscore(y, x = NULL, p00, p11, n.round = 3, bootin = FALSE)
}
\arguments{
\item{y}{Matrix or data frame of observed responses for K RRT items}

\item{x}{Predictor variable values (matrix or data frame).}

\item{p00}{The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)}

\item{p11}{The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)}

\item{n.round}{Number of digits for rounding estimates (default is \code{n.round=3}).}

\item{bootin}{Default: \code{bootin=FALSE} - do not change. For RRbootstrap internal functioning.}
}
\value{
When \code{x=NULL} (default setting), \code{RRsumscore} yields the prevalence of RRT sum scores; otherwise, RRsumscore performs ordinal regression. In both cases, \code{RRsumscore} returns a list:
\itemize{
   \item \code{$estimates}: a data frame with the parameter estimate (\code{est}), standard error (\code{se}), T-statistic (\code{tval}), lower bounds of the 95\% confidence interval (\code{min95}), and upper bounds of the 95\% confidence interval (\code{max95}).
   \item \code{$fitmeasures}: -2*log likelihood (\code{minus2LL}), degrees of freedom (df), Akaike Information Criterion (\code{AIC}), and AIC corrected for (small) sample sizes (\code{AICc}).
}
When the Hessian is singular, standard errors are inestimable, as are the t-statistic, p-value, and confidence limits; in this case, the model will return \code{NaN} for those entries.
}
\description{
\code{RRsumscore} estimates (1) the prevalence of 0..K sensitive traits, (2) ordinal regression when a matrix or data frame of predictor variables is provided.
}
\examples{
## Sum score modeling with simulated data (see ?RRsimulate for questions)
p00 <- p11 <- 5/6
SSdata <- RRsimulate(model="sumscore",pi.1 = c(0.15, 0.2, 0.4),
                     p00=p00, p11=p11, x.betas=c(2), n=375)
head(SSdata)
# Sum score prevalence
RRsumscore(y=SSdata[,1:3], p00=p00, p11=p11, n.round=2)
# Ordinal regression modeling
RRsumscore(y=SSdata[,1:3], x=SSdata[,4], p00=p00, p11=p11)

## Sum score example: the number of species hunted in Southwest China
    # Load the data and prepare the design parameters
data(huntRRdata) # See ?huntRRdata for information on this dataset
p00 <- p11 <- 5/6
    # Sum score prevalence
RRsumscore(y=huntRRdata[,1:3], p00=p00, p11=p11)
    # Ordinal regression: the effect of age on the number of species hunted
RRsumscore(y=huntRRdata[,1:3], x=huntRRdata['Age'], p00=p00, p11=p11)

## Sum score example: University student behaviors
    # Loading the data
data(RRdata)  # See ?RRdata for more information on the sensitive questionnaire.

    # Calculating p00 and p11 for each question
college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order
question.indices <- c(1,6) # Note that the p00 and p11 values MUST be the same
                           # for all questions used in IRT or sumscore analyses

    # Performing sum score modeling
RRsumscore(RRdata[,question.indices], p00=college.p00[1], p11=college.p11[1])

}
\seealso{
\code{\link{RRbootstrap}} for bootstrap estimates when the Hessian is singular.
}
