


#' HySime
#'
#' Hyperspectral signal subspace estimation to find the number of endmembers present in a multispectral or hyperspectral image
#'
#' @param x RasterStack, RasterBrick, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral).
#' @param noisetype Character. Infers the noise in a image, by assuming that the reflectance at a given band is well modelled by
#' a linear regression on the remaining bands. By default the noise type is "additive". Other option is "poisson".
#'
#' @return Number of endmembers or virtual dimensionality of the image, returned as a positive integer.
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @rawNamespace import(pracma, except = c(Norm,squareform,Rank))
#' @import Rfast
#' @rawNamespace import(stats, except = c(predict,density,weighted.mean,aggregate,quantile,update))
#'
#' @usage hysime(x, noisetype)
#'
#' @references J. M. Bioucas-Dias and J. M. P. Nascimento, "Hyperspectral Subspace Identification,"
#' in IEEE Transactions on Geoscience and Remote Sensing, vol. 46, no. 8, pp. 2435-2445, Aug. 2008.
#'
#'
#' @export
#'


hysime <- function(x, noisetype)

{
  if if (class(x)[1] != "RasterStack" && class(x)[1] != "matrix" && class(x)[1] !="RasterBrick") {
    stop("x is not a matrix or raster")
  }

  if (missing(noisetype))
  {
    noisetype <- "additive"
  }


  if (class(x)[1] == "RasterStack" | class(x)[1] =="RasterBrick")
  {
    #transfrom raster in matrix remove NA
    Ma <- Rfast::transpose(stats::na.omit(as.matrix(x)))
  } else{
    Ma <-x[ , colSums(is.na(x)) == 0]
  }




  # rows are bands
  # columns are signals
  p <- dim(Ma)[1]
  N <- dim(Ma)[2]


  #noise estimation function

  estAdditiveNoise <- function(r) {
    p <- dim(r)[1]
    N <- dim(r)[2]
    small <- 1e-6

    w <- matrix(ncol = N, nrow = p, 0)
    RR <- Rfast::Tcrossprod(r, r)
    RRi <- pracma::inv(RR + small * pracma::eye(p))

    for (i in 1:p)
    {
      XX <- RRi - (RRi[, i] %*% t(RRi[i,])) / RRi[i, i]
      RRa <- RR[, i]
      RRa[i] <- 0
      beta <- XX %*% RRa
      beta[i] <- 0
      w[i,] <- r[i,] - t(beta) %*% r

    }

    Rw <- diag(diag(Rfast::Tcrossprod(w, w) / N))

    return(list(w, Rw))

  }

  #noise estimation

  if (noisetype == "poisson")
  {
    u_Ru <- list(NULL)

    sqy <- sqrt(Ma * (Ma > 0))          # prevent negative values
    u_Ru_p <- estAdditiveNoise(sqy)   # noise estimates
    x <- (sqy - u_Ru_p[[1]]) ^ 2                #signal estimates
    u_Ru[[1]] <- sqrt(x) * u_Ru_p[[1]] * 2
    u_Ru[[2]] <- Rfast::Tcrossprod(u_Ru[[1]], u_Ru[[1]]) / N
  } else if (noisetype == "additive") {
    u_Ru <- estAdditiveNoise(Ma) # noise estimates

  } else{
    stop("noisetype is not valid")
  }


  x <- Ma - u_Ru[[1]]

  Ry <- Rfast::Tcrossprod(Ma, Ma) / N   # sample correlation matrix
  Rx <- Rfast::Tcrossprod(x, x) / N  # signal correlation matrix estimates

  #Computing the eigen vectors of the signal correlation matrix
  E_D <- svd(Rx)
  dx <- E_D$d

  #Estimating the number of endmembers
  Rn <- u_Ru[[2]] + sum(diag(Rx)) / p / 10 ^ 5 * pracma::eye(p)

  Py <- diag(t(E_D$u) %*% Ry %*% E_D$u)
  Pn <- diag(t(E_D$u) %*% Rn %*% E_D$u)

  cost_F <- -Py + 2 * Pn
  kf <- sum(cost_F < 0) #number of endmembers
  ind_asc <- order(cost_F)
  Ek <-
    E_D$u[, ind_asc[1:kf]] # Ek matrix which columns are the eigenvectors that span the signal subspace

  return(kf)

}
