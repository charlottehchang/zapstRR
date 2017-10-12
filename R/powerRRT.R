###====================================================================================
### Power analysis
###====================================================================================

#' @export
#' @title Power analysis for individual RRT question designs
#'
#' @description \code{powerRRT} performs power analysis for randomized response survey designs.
#' @param pi.null Null expectation for sensitive trait prevalence
#' @param pi.alt True sensitive trait prevalence under the alternative hypothesis
#' @param p00 The conditional classification probability \eqn{y=0|z=0} (observed no conditional on true no)
#' @param p11 The conditional classification probability \eqn{y=1|z=1} (observed yes conditional on true yes)
#' @param n Sample size
#' @param conf.level Confidence level of the interval (default = 0.95, value must be within (0,1))
#' @param alternative  The alternative hypothesis, must be a character string matching one of the following: \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#'
#' @examples
#' ## An example power analysis for a specific RRT design
#' # Here, we assume that the null is that the sensitive trait is not present.
#' # The alternative is that it is present.
#'   # Simulated vector of possible sample sizes (number of respondents)
#' nsim  <- seq(50,  1000, by = 10)
#'   # Vector of possible sensitive trait prevalences
#' prev  <- c(.02, .05, .1, .2, 0.3)
#'   # Assuming a null hypothesis of no sensitive trait possession in population
#' pi.null <- 0
#'   # Question design features for randomized response
#' p00 <- p11 <- 5/6
#' RRpwr <- matrix(0, length(nsim), length(prev))
#' for (j in 1:length(prev)) {
#'     pi.alt <- prev[j]
#'     for (i in 1:length(nsim)) {
#'          n <- nsim[i]
#'          RRpwr[i,j] <- powerRRT(pi.null, pi.alt, p00, p11, n, alternative="greater")
#'      }
#' }
#'      ## Plotting
#' colors <- rainbow(length(prev))
#' par(mar=c(5,4,2,2),mfrow=c(1,1))
#' plot(nsim, RRpwr[, 1], ylim=c(0,1), type="l",
#'      xlab="Sample size (n)", ylab="Power", col=colors[1], lwd=1.25)
#' grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = F)
#' for (j in 2:length(prev)) {
#'      lines(nsim, RRpwr[,j], col=colors[j])
#' }
#' legend("bottomright",legend=prev,
#'        title=expression(paste("Prevalence: (",pi[alt],")")),
#'        lty=rep(1,length(prev)),bty="n",horiz=TRUE,
#'        col=colors, title.adj=1, lwd=1.25)
#' @return \code{powerRRT} returns the statistical power for a given RRT survey design and null hypothesis.
#' @export
powerRRT <- function(pi.null=NULL, pi.alt=NULL, p00=NULL, p11=NULL, n=NULL, conf.level=0.95, alternative=c("greater","less","two.sided")) {
  ## Error handling
  if (!is.numeric(p00) | !is.numeric(p11) | !is.numeric(pi.null) | !is.numeric(pi.alt)) {
    stop("You must provide pi.null, pi.alt, p00, and p11. All pi.null, pi.alt, p00, and p11 values must be numeric and within the interval (0, 1).")
  }

  pi.check <- findInterval(c(p00, p11, pi.null, pi.alt), c(0,1))
  pi.check <- ifelse(pi.check != 1, 0, 1)
  if ( sum(pi.check) != 4) {
    stop("All pi.null, pi.alt, p00, and p11 values must be within the interval (0,1).")
  }

  if (is.null(n)) {
    n <- 150
    warning("Nothing provided for n, setting n to 150")
  }

  if (alternative != "greater" & alternative != "less" & alternative != "two.sided") {
    message("Please provide one of the following: 'greater', 'less', or 'two.sided' for alternative; setting to 'greater'.")
    alternative <- "greater"
  }

  ## Generate RRT classification process
  P  <- matrix(c(p00, 1-p00, 1-p11, p11), 2, 2)

  ## Calculating variance
  var0    <- var.pwr(pi=pi.null, P=P, n=n)
  var.alt <- var.pwr(pi=pi.alt, P=P, n=n)

  ## Calculating power
  z.score <- qnorm(conf.level)
  if (alternative=="greater") {
    pwr.score <- pnorm( (pi.alt - (pi.null+z.score*sqrt(var0)))/sqrt(var.alt))
  } else if (alternative=="less") {
    pwr.score <- pnorm( (pi.null - (pi.alt+z.score*sqrt(var.alt)))/sqrt(var0))
  } else if (alternative=="two.sided") {
    z.score <- qnorm(1-((1-conf.level)/2))
    upper <- pnorm( (pi.alt - (pi.null+z.score*sqrt(var0)))/sqrt(var.alt))
    lower <- pnorm( (pi.null - (pi.alt+z.score*sqrt(var.alt)))/sqrt(var0))
    pwr.score <- upper + lower
  } else {
    stop("Please specify 'greater', 'less', or 'two.sided' for alternative.")
  }
  return(pwr.score)
}
