% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/suppoRRT.R
\docType{data}
\name{RRdata}
\alias{RRdata}
\title{Randomized Response Technique (RRT) Survey for a set of sensitive behaviors among university students}
\format{A data frame containing a sample of 710 observations.
Each question yielded a yes (\code{1}) or no (\code{0}) response. The RRT question topics are:
\itemize{
 \item copied: Have you ever copied in an exam?
 \item fought: Have you ever fought with a teacher?
 \item bullied: Have you been bullied?
 \item bullying: Have you ever bullied someone?
 \item drug: Have you ever taken drugs on campus?
 \item sex: Have you had sex on the premises of the university?
}}
\usage{
data(RRdata)
}
\description{
Undergraduate students were asked about a variety of sensitive topics including cheating, fighting, and bullying using the unrelated question (Horvitz) design.
}
\examples{
## Loading the data
data(RRdata)

## Calculating p00 and p11 for each question (see: ?RRunivariate)
college.py <- c(1/12, 1/10, 20/30, 1/10, 10/30, 1/12)
college.p00 <- 1-0.5*college.py # p00 for each question prompt in order
college.p11 <- 1-0.5*(1-college.py) # p11 for each question prompt in order

## Performing univariate RRT inference
   # How many students copied exams? (Question 1)
RRunivariate(RRdata['copied'], p00=college.p00[1], p11=college.p11[1])
   # How many students took drugs on campus? (Question 5)
RRunivariate(RRdata['drug'], p00=college.p00[5], p11=college.p11[5])

## Performing multivariate inference - note that the same values of p00 and p11 are required
question.indices <- c(1,6)
    # Sum score model
RRsumscore(RRdata[,question.indices], p00=college.p00[1], p11=college.p11[1])
    # Item response theory model
RRirt(RRdata[,question.indices],p00=college.p00[1], p11=college.p11[1])

}
\references{
Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
\emph{The unrelated question RR model: Theoretical framework.}
Journal of the American Statistical Association, 64, 520-539.

Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
\emph{The unrelated question RR model.}
 Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.

Cobo Rodriguez, B., del Mar Rueda Garcia, M., Arcos Cebrian, A. (2015).
 \emph{RRTCS: Randomized Response Techniques for Complex Surveys.}
 \url{https://CRAN.R-project.org/package=RRTCS}
}
\seealso{
\code{\link[RRTCS]{HorvitzDataRealSurvey}} for the original RRdata source,
\code{huntRRdata} for RRT data on illegal hunting in Southwest China.
}
\keyword{datasets}
