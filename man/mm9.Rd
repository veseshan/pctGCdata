\name{mm9}
\alias{mm9}
\title{GC percentages for mouse genome build mm9}
\description{
  GC percentages in 1 kiolobase sliding windows along chromosomes 1-19,
  X, Y & MT for mouse genome build mm9.
}
\usage{data(mm9)}
\format{
  A list containing 25 vectors corresponding to chromosomes 1-19, X, Y &
  MT in that order. For example \code{mm9[[7]]} gives the GC percentage
  for chromosome 7.
}
\details{
  GC percentage is the proportion of nucleotides that are GC in a 1000
  base window. The windows begin with bases 1-1000 and slide 100 bases
  at a time. See the function \code{getGCpct} to get GC percentages for
  a subset of genomic positions.
}
\keyword{datasets}
