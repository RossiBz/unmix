
#' Vertex Component Analysis
#'
#' Extracts endmember signatures from multispectral or hyperspectral image by using the Vertex Component Analysis method (VCA).
#'
#' @param x RasterStack, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral).
#' @param q Positive integer. number of endmembers in the scene
#'
#' @return Estimated mixing matrix (endmembers signatures)
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @import pracma
#' @import RSpectra
#' @import Rfast
#' @import stats
#'
#' @usage vca(x, q)
#'
#' @references José M. P. Nascimento and José M. B. Dias "Vertex Component Analysis: A Fast Algorithm to Unmix Hyperspectral Data"
#' submited to IEEE Trans. Geosci. Remote Sensing, vol. .., no. .., pp. .-., 2004
#'
#'
#' @export
#'



vca <- function(x, q) {
  if (class(x)[1] != "RasterStack" && class(x)[1] != "matrix") {
    stop("x is not a matrix or raster")
  }

  if (class(q)[1] == "numeric" && q > 0) {
    q <- q
  } else{
    stop("q is not a positive integer")
  }


  if (class(x)[1] == "RasterStack")
  {
    #transfrom raster in matrix remove NA
    R <- transpose(na.omit(as.matrix(x)))
  } else{
    R <-x[ , colSums(is.na(x)) == 0]
  }

  # rows are bands
  # columns are signals
  p <- dim(R)[1]
  N <- dim(R)[2]


  #SNR estimates --------------------------------------------------------------
  R_o <- R - rowmeans(R)[row(R)]       # data with zero-mean
  D <-
    svds(Tcrossprod(R_o, R_o) / N, q)    # computes the p-projection matrix
  x_p <-
    Crossprod(D$u, R_o)                 # project the zero-mean data onto p-subspace


  r_m <- rowmeans(R)
  P_y <- sum(R ^ 2) / N
  P_x <- sum(x_p ^ 2) / N + t(r_m) %*% r_m

  SNR <- 10 * log10((P_x - q / p * P_y) / (P_y - P_x))

  SNR_th <- 15 + 10 * log10(q)


  # Choosing Projective Projection or projection to p-1 subspace-----------------


  if (SNR < SNR_th)
  {
    d <- q - 1
    Ud <- D$u[, 1:d]

    Rp <-
      mat.mult(Ud, x_p[1:d, ]) + rowmeans(R)[row(R)]     # again in dimension L

    x_ <- x_p[1:d, ]
    c <- max(colSums(x_ ^ 2)) ^ 0.5
    y <- rbind(x_, c * ones(1, N))

  } else{
    d <- q
    D <-
      svds(Tcrossprod(R, R) / N, d)         # computes the p-projection matrix

    x_ <-  Crossprod(D$u, R)
    Rp <-
      mat.mult(D$u, x_[1:d, ])     # again in dimension L (note that x_p has no null mean)

    u <- matrix(colSums(x_ * rowmeans(x_)), nrow = 1)

    y <-  x_ / (u[rep(seq_len(nrow(u)), each = d),])




  }



  # VCA algorithm-------------------------------


  indice <- matrix(0, 1, q)
  A <- matrix(0, q, q)
  A[q, 1] <- 1

  for (i in 1:q)
  {
    w <- rand(q, 1)
    f <- w - A %*% pinv(A) %*% w
    f <- f / sqrt(sum(f ^ 2))

    v <- Crossprod(f, y)

    indice[i] <- which.max(abs(v))

    A[, i] <- y[, indice[i]]
  }

  Ae <- Rp[, indice]

  return(Ae)

}
