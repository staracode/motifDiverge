\name{Scoring matrices for TFBS motifs}
\alias{pscmToPspm}
\alias{pspmToPssm}
\alias{pscmToPssm}
\title{ 
	Convert position specific count, frequency or probability matrices into scoring matices.
}
\description{
TFBS motifs are often described as position specific count (pscm), probability (pspm) scoring (pssm) matices. 
The latter are typically derived from the former.
}
\usage{
pscmToPspm(pscm, bg=NULL, digits=4, reg=TRUE)
pspmToPssm(pspm, bg=NULL, digits=4, return.pspm=FALSE)
}
\arguments{
  \item{pscm}{
	Position specific count matrix.
}
  \item{pspm}{
        Position specific probability matirx.
}
  \item{bg}{
        The background distribution.
}
  \item{digits}{
        How many digitis should the pspm/pssm be based on. 
}
 \item{reg}{
        Wheter a bg-distributed pseudocount should be added.
}
 \item{return.pspm}{
	Should the \code{pspm} be returned?
}
}
\details{
\code{pscmToPspm} returns adds a pseudocount (if \code{reg} is \code{TRUE}) to \code{pscm} and returns a \code{pspm} with \code{digits} 
positions after the decimal point.
\code{pspmToPssm} takes a \code{pspm} and returns the log-odds score of \code{pspm} compared to \code{bg}. To avoid singularities 
\code{pspm} is regularized if it contains entries smaller than 10^(-\code{digits}).
\code{pscmToPssm} Is just a wrapper.
}
\value{
\code{pscmToPspm} returns a pspm matrix. 
\code{pspmToPssm} returns a pssm matrix, or a list with two matrices if \code{return.pspm} is \code{TRUE}. 
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
	#- generate random pscm
	pscm = matrix(sample(1:20,4*8,rep=TRUE),nrow=4)
	pspm = pscmToPspm(pscm)
	pssm = pspmToPssm(pspm)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
