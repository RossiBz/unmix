

#' Noise-whitened Harsanyi–Farrand–Chang method
#'
#' Finds the number of endmembers present in a multispectral or hyperspectral image by using the noise-whitened Harsanyi–Farrand–Chang (NWHFC) method.
#'
#' @param x RasterStack, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral),
#' @param noisewhiten Logical. If TRUE, performs noise-whitening before computing the number of endmembers.
#' FALSE consist in the Harsanyi–Farrand–Chang (HFC) method. By default TRUE.
#' @param far  Positive integer.  Probability of false alarm . By default 10^-3
#'
#' @return Number of endmembers in the image, returned as a positive integer.
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @import pracma
#' @import Rfast
#' @import stats
#'
#' @usage hfcvd(x,noisewithen, far)
#'
#' @references C.-I. Chang and Q. Du, "Estimation of number of spectrally distinct
#' signal sources in hyperspectral imagery," IEEE Transactions on Geoscience and Remote Sensing, vol. 43, no. 3, mar 2004.
#
#'
#'
#'
#' @export
#'





hfcvd <- function(x, noisewhiten, far)

{
  if (class(x)[1] != "RasterStack" && class(x)[1] != "matrix") {
    stop("x is not a matrix or raster")
  }

  if (missing(noisewhiten))
  {
    noisewhiten <- TRUE #default
  }

  if (missing(far))
  {
    far <- 10 ^ -3
  }


  if (class(far)[1] == "numeric" && far > 0 && far <= 1) {
    far <- far
  } else{
    stop("far is not a positive scalar in the range (0,1]")
  }



  if (class(noisewhiten)[1] == "logical") {
    noisewhiten <- noisewhiten
  } else{
    stop("noisewhiten is not logical")
  }



  if (class(x)[1] == "RasterStack")
  {
    #transfrom raster in matrix remove NA
    Ma <- stats::na.omit(as.matrix(x))
  } else{
    Ma <- stats::na.omit(Rfast::transpose(x))
  }



  numBands <- dim(Ma)[2]
  N <- dim(Ma)[1]


  if (noisewhiten)
  {
    # Noise Whiten the data

    # Covariance matrix
    covMatrix <- Rfast::cova(Ma)

    # Inverse covariance matix
    invCovMatrix <- pracma::inv(covMatrix)


    # Diagonal elements of inverse covariance matrix
    zeta <- diag(invCovMatrix)

    # Noise covariance matrix
    noiseCovMatrix <- diag(1 / zeta)
    # Noise whitening
    I <- Rfast::mat.mult(Ma, pracma::sqrtm(noiseCovMatrix)$Binv)
  } else{
    I <- Ma
  }



  # calculate eigenvalues of covariance and correlation between bands
  lambda_cov <-
    sort(eigen(Rfast::cova(I), only.values = TRUE)$values, decreasing = TRUE)
  lambda_corr <-
    sort(eigen(Rfast::Crossprod(I, I) / N, only.values = TRUE)$values, decreasing =
           TRUE)


  # Difference between eigen values
  diff <- lambda_corr - lambda_cov

  # Estimate noise variance
  variance <- (2 / N) * (lambda_corr ^ 2 + lambda_cov ^ 2)

  # Binary hypothesis test
  tau <- sqrt(2 * variance) * pracma::erfcinv(2 * far)

  numEndmembers <- sum(diff > tau)

  return(numEndmembers)

}
