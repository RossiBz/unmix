

#' Pixel purity index
#'
#' Extracts endmember signatures from multispectral or hyperspectral image by using the pixel purity index (PPI) algorithm.
#'
#' @param x RasterStack, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral),
#' @param q Positive integer. number of endmembers in the scene
#' @param numSkewers  Positive integer. Number of randomly generated unit vectors to be used for projection. By default 10^4
#' @param reductionmethod  Character. Select a dimensionality reduction method. By default maximum noise fraction "MNF". Other options are
#' principal component analysis "PCA" or "None".
#'
#' @return Estimated mixing matrix (endmembers signatures)
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @import pracma
#' @import Rfast
#' @import stats
#'
#' @usage ppi(x, q, numSkewers, reductionmethod)
#'
#' @references J.W Boardman, F.A. Kruse and R.O. Green, "Mapping target signatures via partial unmixing of AVIRIS data.",
#' Technical Report, California, USA, 1995.
#'
#'
#'
#' @export
#'




ppi <- function(x, q, numSkewers, reductionmethod) {
  if (class(x)[1] != "RasterStack" && class(x)[1] != "matrix") {
    stop("x is not a matrix or raster")
  }

  if (class(q)[1] == "numeric" && q > 0) {
    q <- q
  } else{
    stop("q is not a positive integer")
  }

  if (missing(numSkewers))
  {
    numSkewers <- 10 ^ 4 #default
  }


  if (class(numSkewers)[1] == "numeric" &&
      numSkewers > 0) {
    numSkewers <- numSkewers
  } else{
    stop("numSkewers is not a positive integer")
  }



  if (class(x)[1] == "RasterStack")
  {
    #transfrom raster in matrix remove NA
    Ma <- stats::na.omit(as.matrix(x))
  } else{
    Ma <- stats::na.omit(transpose(x))
  }


  if (missing(reductionmethod))
  {
    reductionmethod <- "MNF" #default
  }



  if (reductionmethod == "PCA")
  {
    #Principal component ------------------
    R <- Rfast::transpose(stats::prcomp(Ma, rank. = q)$x)

  } else if (reductionmethod == "MNF") {
    #Minimum Noise Fraction ------------------

    # Find the covariance of the noise.

    V <- svd(Rfast::cova(base::diff(Ma)))


    # Consider only nonzero eigen values
    getNonZeros <- sum(!(V$d == 0))


    U <- V$d[1:getNonZeros]


    if (isempty(U)) {
      # Noise whitening matrix when all the eigen values are zero.
      noiseWhitening <- Rfast::transpose(V$u)
    } else{
      # Noise whitening matrix
      noiseWhitening <- Rfast::transpose(V$u)

      # Consider only non zero eigen values.
      nonzeroRotations <- base::solve(t(sqrt(diag(U))),
                                noiseWhitening[1:getNonZeros, 1:getNonZeros])
      noiseWhitening[1:getNonZeros, 1:getNonZeros] <-
        nonzeroRotations
    }

    # Get the noise adjusted data.


    cube <-  Rfast::Tcrossprod(noiseWhitening, Ma - colmeans(Ma)[col(Ma)])

    # Second rotation by principal components.
    R <-
      Rfast::transpose(prcomp(
        Rfast::transpose(cube),
        rank. = q,
        center = FALSE
      )$x)

  } else if (reductionmethod == "None") {
    #no reduction method ----------------------

    R <- Rfast::transpose(Ma)

  } else {
    stop("reductionmethod is not valid")
  }




  # Run ppi -----------------------------

  # rows are bands
  # columns are signals
  p <- dim(R)[1]
  N <- dim(R)[2]

  # Generate random unit vectors.
  skewers <- pracma::randn(numSkewers, p)


  votes <- rep(0, N)

  for (kk in 1:numSkewers)
  {
    # skewers[,kk] is already a row
    tmp <- abs(skewers[kk,] %*% R)
    idx <- which.max(tmp)
    votes[idx] <- votes[idx] + 1
  }

  # Find the maximum projected votes to decide the endmember signatures.
  max_idx <- order(votes, decreasing = TRUE)

  # the last right idx..s at the max_idx list are
  # those with the max votes
  endmember_idx <- max_idx[1:q]
  endmember_spectra <- Rfast::transpose(Ma[endmember_idx,])

  return(endmember_spectra)

}
