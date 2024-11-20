\name{Parameters for correlated Binomial trials}
\alias{cbernEstimateModelPars}
\title{ 
	Model parameters for differences in motif hit-counts.
}
\description{
	Given a set of sequence-pairs motif hits are annotated. Then parameters for
	a statistical model for motif hit-count differences are estimated.

}
\usage{
	cbernEstimateModelPars(seqs.x, seqs.y, pssm, pspm, cut.fw, cut.rc, indep=TRUE, useCounts=TRUE, modList=NULL, NSAMP=1E+5)
}
\arguments{
  \item{seqs.x}{ 
	DNAStringSet: Sequences to be compared.
}
  \item{seqs.y}{
	DNAStringSet: Sequences to be compared. 
}
  \item{pssm}{
        The pssm for motif hit annotation.
}
 \item{pspm}{
        The pspm for motif hit annotation.
}
\item{cut.fw}{
        The cutoff for motif hit annotation on the forward strand.
}
\item{cut.rc}{
        The cutoff for motif hit annotation on the reverse strand.
}
\item{indep}{
        Are the sequence pairs independent or homologous?
}
  \item{useCounts}{
        Return count-based or model-based estimates?
}
 \item{modList}{
        If model-based estimates are to be calculated this needs to be a list of tree models suitable to 
	generate the phyloHMM.
}
 \item{NSAMP}{
	Number of samples used for paramter estimates in case of model-based parameter estimation.
}
}
\details{
\code{modList} needs to be a list of tree models used by the \code{rphast} package to build a phyloHMM.
}
\value{
A matrix (one row for each sequence pair) containing model estimates and motif hit counts.
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
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
