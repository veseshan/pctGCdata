\name{hg38}
\alias{hg38}
\title{GC percentages for human genome build hg38}
\description{
  GC percentages in 1 kiolobase sliding windows along chromosomes 1-22,
  X, Y & MT for human genome build hg38.
}
\usage{data(hg38)}
\format{
  A list containing 25 vectors corresponding to chromosomes 1-22, X, Y &
  MT in that order. For example \code{hg38[[7]]} gives the GC percentage
  for chromosome 7.
}
\details{
  GC percentage is the proportion of nucleotides that are GC in a 1000
  base window. The windows begin with bases 1-1000 and slide 100 bases
  at a time. See the function \code{getGCpct} to get GC percentages for
  a subset of genomic positions.
}
\keyword{datasets}
