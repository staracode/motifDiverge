\name{Logg odds score distribution}
\alias{scoreDist}
\alias{scoreCutFromErr}
\alias{errFromScoreCut}
\title{ 
	Distribution of the log-odds scores of motifs related to TFBSs.
}
\description{
Distribution of log-odds scores for a TFBS motif, given a position specific scoring matrix, 
a position specific probability matrix and a background model. Methods for calculating 
corrsponding score cutoffs for Type I, Type II and balanced errors. 
}
\usage{
scoreDist(pssm,pspm=rep(1/4,4),gran=0.01)
scoreCutFromErr(err,pssm,pspm,bg=rep(1/4,4),type="type1",gran=0.01)
errFromScoreCut(cut,pssm,pspm,bg=rep(1/4,4),type="type1",gran=0.01)
}
\arguments{
  \item{pssm}{
	Position specific scoring matrix.
}
  \item{pspm}{
        Position specific probability matirx. Columns must sum to one, dimensions must match \code{pssm}. 
}
  \item{bg}{
	Background frequencies. 4-vector, must sum to one.
}
  \item{err}{
	Error of the type \code{type}, \eqn{(0\leq \mbox{err}\leq 1)}{(0 <= err <= 1)}. 
}
  \item{cut}{
	Log-odds score cutoff.
}
  \item{type}{
	Denotes the type of error. Possible values are "type1", "type2" and "balanced".
}
  \item{gran}{
        The granularity with which the scores in \code{pssm} are approximated.
}
}
\details{
	For \code{scoreDist} \code{pspm} can be a 4-vector of background frequencies. 
	For \code{scoreCutFromErr} the argumetn \code{err} is ignored when \code{type} is "balanced".
	For \code{errFromScoreCut} "balanced" is not a valid value for \code{type}
}
\value{
\code{scoreDist} returns a matrix with two colums: The first contains log-odds score values, the second the according probabilities. 
\code{scoreCutFromErr} returns a cutoff value according to the value of \code{err}. 
\code{errFromScoreCut} returns an error value according to the value of \code{cut}. 
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
	#- generate random pspm and pssm
	pspm = matrix(runif(4*8),nrow=4)
	pspm = t(t(pspm)/colSums(pspm))
	pssm = log2(pspm) - log2(1/4)

	#- get the score distributions
	sdbg = (scoreDist(pssm))
	sdfg = (scoreDist(pssm=pssm,pspm=pspm))

	#- visualize
	plot(sdfg,col="red",pch=19,xlab="Log odds score",ylab="Probability")
	points(sdbg,col="blue",pch=19) ; grid()

	#- calculate the different cutoffs
	t1cut = scoreCutFromErr(.05,pssm,pspm)
	t2cut = scoreCutFromErr(.05,pssm,pspm,type="type2")
	blcut = scoreCutFromErr(.05,pssm,pspm,type="balanced")
	
	#- add them to plot
	abline(v=t1cut,col="blue")
	abline(v=t2cut,col="red")
	abline(v=blcut,col="green")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
