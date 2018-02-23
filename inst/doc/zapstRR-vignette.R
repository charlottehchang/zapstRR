## ----uni_bulbul----------------------------------------------------------
###==================================================================================
## Preparation
###==================================================================================

## Load the package
library(zapstRR)

## Load the data
data("huntRRdata")

## Study design parameters
p00 <- p11 <- 5/6 # p00: Probability of observing no conditional on a true no; p11: Probability of observing yes conditional on a true yes.

###==================================================================================
## Univariate prevalence analysis: Bulbul hunting
###==================================================================================

RRunivariate(y=huntRRdata['Bulbul'], p00=p00, p11=p11)$estimates

## ----bulbul_glm----------------------------------------------------------
###==================================================================================
## Univariate logistic regression: Bulbul hunting
###==================================================================================
RRunivariate(y=huntRRdata['Bulbul'], x=huntRRdata['Age'], p00=p00, p11=p11)$estimates

## ----barbet_glm_boot-----------------------------------------------------
###==================================================================================
## Univariate logistic regression bootstrap: Barbet hunting
###==================================================================================
# Barbet has low prevalence, which combined with the limited sample size (n=189),
# precludes SE identification in RRunivariate.
RRunivariate(y=huntRRdata['Barbet'], x=huntRRdata['Age'], p00=p00, p11=p11)$estimates

# Compare the results against the bootstrap estimate:
RRbootstrap(model='univariate',y=huntRRdata[,'Barbet'], x=huntRRdata['Age'], p00=p00, p11=p11, reps=50)

## ----bulbul-GLM-gof------------------------------------------------------
###==================================================================================
## Goodness of fit test example: Bulbul hunting
###==================================================================================
reps <- 50 # This should be increased for substantive analyses; kept at this low number for convenience
gof.bulbul <- RRgof(model="univariate", y=huntRRdata['Bulbul'], x=huntRRdata['Age'], p00=p00, p11=p11, reps=reps)
gof.bulbul$gof.test.result

# For more information use the command below. (Delete the leading hash symbol, #)
# ?RRgof 

## ----sumscore_prev-------------------------------------------------------
###==================================================================================
## Sum score prevalence: The number of species hunted
###==================================================================================

RRsumscore(y=huntRRdata[,1:3], p00=p00, p11=p11)$estimates

## ----sumscore_university-------------------------------------------------
###==================================================================================
## Sum score prevalence: University student behaviors
###==================================================================================

## Loading the data
data(RRdata)

## Calculating p00 and p11 for each question (see: ?RRunivariate)
college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order

## Performing multivariate inference - note that the same values of p00 and p11 are required
question.indices <- c(1,6)
    # Sum score model
RRsumscore(RRdata[,question.indices], p00=college.p00[1], p11=college.p11[1])

## ----ordinal_reg---------------------------------------------------------
###==================================================================================
## Ordinal regression: Age and the number of species hunted
###==================================================================================

RRsumscore(y=huntRRdata[,c(2,1,3)], x=huntRRdata['Age'], p00=p00, p11=p11)

## ----ordinal_bootstrap---------------------------------------------------
###==================================================================================
## Ordinal regression bootstrap: Age and the number of species hunted
###==================================================================================

RRbootstrap(model="sumscore", y=huntRRdata[,c(2,1,3)], x=huntRRdata['Age'], p00=p00, p11=p11, reps=50)

## ----irt_prev------------------------------------------------------------
###==================================================================================
## IRT prevalence: rate of hunting each bird species
###==================================================================================

IRT.MLEs <- RRirt(y=huntRRdata[,c(1:3)], p00=p00, p11=p11)$estimates
IRT.MLEs[,c(1,2,5,6)] <- IRT.MLEs[,c(1,2,5,6)]*100
signif(IRT.MLEs, 3)

## ----uni_prev3-----------------------------------------------------------
###==================================================================================
## Univariate prevalence of hunting
###==================================================================================

uniMLEs <- rbind(RRunivariate(y=huntRRdata['Barbet'], p00=p00, p11=p11)$estimates, # Barbet
                 RRunivariate(y=huntRRdata['Bulbul'], p00=p00, p11=p11)$estimates, # Bulbul
                RRunivariate(y=huntRRdata['Partridge'], p00=p00, p11=p11)$estimates) # Partridge
rownames(uniMLEs) <- c("Barbet","Bulbul","Partridge")
uniMLEs[,c(1,2,5,6)] <- uniMLEs[,c(1,2,5,6)]*100
signif(uniMLEs,3)

## ----irt_glm-------------------------------------------------------------
###==================================================================================
## IRT regression: The effect of age on hunting across multiple bird species
###==================================================================================

RRirt(y=huntRRdata[,c(3,1,2)], x=huntRRdata['Age'], p00=p00, p11=p11)

## ----irt_phi-------------------------------------------------------------
###===========================================================================
### NONCOMPLIANCE
###===========================================================================

## Prevalence only
RRirt(y=huntRRdata[,1:3], p00=p00, p11=p11, noncompliance = TRUE)

## IRT regression
RRirt(y=huntRRdata[,1:3], x=huntRRdata['Age'],p00=p00, p11=p11, noncompliance = TRUE)

## Bootstrap
RRbootstrap(model="irt", y=huntRRdata[,1:3], x=huntRRdata['Age'],p00=p00, p11=p11, noncompliance = TRUE, reps=50)

## ----power---------------------------------------------------------------
###===========================================================================
### POWER ANALYSIS
###===========================================================================

## An example power analysis for a specific RRT design
op <- par()
# Here, we assume that the null is that the sensitive trait is not present. The alternative is that it is present.
nsim  <- seq(50,  1000, by = 10)  # Simulated vector of possible sample sizes (number of respondents)
prev  <- c(.02, .05, .1, .2, 0.3) # Vector of possible sensitive trait prevalences
pi.null <- 0                      # Assuming a null hypothesis of no sensitive trait possession in population
p00 <- p11 <- 5/6                 # Question design features for randomized response
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
plot(nsim, RRpwr[, 1], ylim=c(0,1), type="l", xlab="Sample size (n)", ylab="Power", col=colors[1], lwd=1.25)
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = F)
for (j in 2:length(prev)) {
     lines(nsim, RRpwr[,j], col=colors[j])
}
legend("bottomright",legend=prev,title=expression(paste("Prevalence: ",pi[alt],"")),lty=rep(1,length(prev)),bty="n",horiz=TRUE, col=colors, title.adj=1, lwd=1.25)

suppressWarnings( par(op) ) # reset plotting parameters

## ----RRsimulate----------------------------------------------------------
###===========================================================================
### NONCOMPLIANCE
###===========================================================================

## UNIVARIATE SIMULATION - assume a logistic regression covariate of 2 for y ~ x
UNIdataLR <- RRsimulate(model="univariate",pi.1 = 0.35, p00=4/5, p11=4/5, x.betas=c(2), n=300)
head(UNIdataLR)

## Analysis
RRunivariate(y=UNIdataLR[,1], x=UNIdataLR[,2], p00=4/5, p11=4/5)

