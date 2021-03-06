% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerRRT.R
\name{powerRRT}
\alias{powerRRT}
\title{Power analysis for individual RRT question designs}
\usage{
powerRRT(pi.null = NULL, pi.alt = NULL, p00 = NULL, p11 = NULL,
  n = NULL, conf.level = 0.95, alternative = c("greater", "less",
  "two.sided"))
}
\arguments{
\item{pi.null}{Null expectation for sensitive trait prevalence}

\item{pi.alt}{True sensitive trait prevalence under the alternative hypothesis}

\item{p00}{The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)}

\item{p11}{The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)}

\item{n}{Sample size}

\item{conf.level}{Confidence level of the interval (default = 0.95, value must be within (0,1))}

\item{alternative}{The alternative hypothesis, must be a character string matching one of the following: \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.}
}
\value{
\code{powerRRT} returns the statistical power for a given RRT survey design and null hypothesis.
}
\description{
\code{powerRRT} performs power analysis for randomized response survey designs.
}
\examples{
## An example power analysis for a specific RRT design
# Here, we assume that the null is that the sensitive trait is not present.
# The alternative is that it is present.
  # Simulated vector of possible sample sizes (number of respondents)
nsim  <- seq(50,  1000, by = 10)
  # Vector of possible sensitive trait prevalences
prev  <- c(.02, .05, .1, .2, 0.3)
  # Assuming a null hypothesis of no sensitive trait possession in population
pi.null <- 0
  # Question design features for randomized response
p00 <- p11 <- 5/6
RRpwr <- matrix(0, length(nsim), length(prev))
for (j in 1:length(prev)) {
    pi.alt <- prev[j]
    for (i in 1:length(nsim)) {
         n <- nsim[i]
         RRpwr[i,j] <- powerRRT(pi.null, pi.alt, p00, p11, n, alternative="greater")
     }
}
     ## Plotting
colors <- rainbow(length(prev))
par(mar=c(5,4,2,2),mfrow=c(1,1))
plot(nsim, RRpwr[, 1], ylim=c(0,1), type="l",
     xlab="Sample size (n)", ylab="Power", col=colors[1], lwd=1.25)
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = F)
for (j in 2:length(prev)) {
     lines(nsim, RRpwr[,j], col=colors[j])
}
legend("bottomright",legend=prev,
       title=expression(paste("Prevalence: (",pi[alt],")")),
       lty=rep(1,length(prev)),bty="n",horiz=TRUE,
       col=colors, title.adj=1, lwd=1.25)
}
