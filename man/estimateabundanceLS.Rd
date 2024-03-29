% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateabundanceLS.R
\name{estimateabundanceLS}
\alias{estimateabundanceLS}
\title{Estimate abundance maps: spectral unmixing}
\usage{
estimateabundanceLS(x,endmembers, method)
}
\arguments{
\item{x}{RasterStack, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral),}

\item{endmembers}{Matrix with dimensions L(bands) x q (endmembers)}

\item{method}{Character. Specifies the least-squares method to use for estimating the abundance maps. By default Fully constrained least-squares method "fcls".
Other options are Unconstrained least-squares method "ucls" and 'Nonnegative constrained least-squares method "ncls".}
}
\value{
Endmember abundance map as RasterStack or matrix. Each layer displays a different endmember.
}
\description{
Estimates abundance maps  from an image and a set of endmember spectral signatures using least square method.
}
\references{
Nirmal Keshava  & John Mustard, Spectral Unmixing, Signal
Processing Magazine, IEEE. 19. 44 - 57. 10.1109/79.974727, 2002. S.M. Kay, Fundamentals of Statistical Signal Processing: Estimation
Theory, Englewood Cliffs, NJ: Prentice Hall, 1993.
}
\author{
Christian Rossi christian.rossi1990@gmail.com
}
