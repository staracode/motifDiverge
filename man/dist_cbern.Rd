\name{dist_cbern}
\alias{dcbern}
\alias{pcbern}
\alias{rcbern}
\title{ 
	Distribution of the difference between the number of sucesses in two correlated sequences of Bernoulli trials.
}
\description{
Probability mass function, cumulative distribution function and random generation of the difference of the number of successes between two correlated sequences of Bernoulli trials, marginalized over the number of sucesses.
}
\usage{
dcbern(k, n1, n2, p10, p01, p11, g.approx=TRUE)
pcbern(k, n1, n2, p10, p01, p11, lower.tail=TRUE, g.approx=TRUE)
rcbern(n, n1, n2, p10, p01, p11)
}
\arguments{
  \item{k}{
	Integer: Difference of successes between the trials.
}
  \item{n}{
        Integer: The number of observations.
}
  \item{n1}{
	Integer: Number of trials in the first sequence of Bernoulli trials.
}
  \item{n2}{
	Integer: Number of trials in the second sequence of Bernoulli trials.
}
  \item{p10}{
	Numeric: Joint probability of success in the first trial and no-success in the second trial.
}
  \item{p01}{
	Numeric: Joint probability of no-success in the first trial and success in the second trial.
}
  \item{p11}{
        Numeric: Joint probability of success in the first trial and success in the second trial.
}
 \item{g.approx}{
	Should a Gaussian approximation be employed when appropriate?
}
 \item{lower.tail}{
        Logical: Should the lower tail P[K <= k] be calculated?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{dcbern} returns the probability mass \eqn{P(k)}{P(k)}, 
\code{pcbern} returns \eqn{P(x\leq k)}{P(x<=k)} and 
\code{rcbern} returns random observations of \eqn{k}{k}. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	k.obs   = rcbern(10   ,100,110,.01*.99,.01*.99,.0001) 
	d.k.obs = dcbern(k.obs,100,110,.01*.99,.01*.99,.0001) 
	p.k.obs = pcbern(k.obs,100,110,.01*.99,.01*.99,.0001) 
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
