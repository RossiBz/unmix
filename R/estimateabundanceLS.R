


#' Estimate abundance maps: spectral unmixing
#'
#' Estimates abundance maps  from an image and a set of endmember spectral signatures using least square method.
#'
#' @param x RasterStack, RasterBrick, matrix with dimensions L(bands) x N(pixels). Remote sensing imagery (usually hyperspectral),
#' @param endmembers Matrix with dimensions L(bands) x q (endmembers)
#' @param method  Character. Specifies the least-squares method to use for estimating the abundance maps. By default Fully constrained least-squares method "fcls".
#' Other options are Unconstrained least-squares method "ucls" and 'Nonnegative constrained least-squares method "ncls".
#'
#' @return Endmember abundance map as RasterStack or matrix. Each layer displays a different endmember.
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @rawNamespace import(pracma, except = c(Norm,squareform,Rank))
#' @import Rfast
#'
#' @usage estimateabundanceLS(x,endmembers, method)
#'
#' @references Nirmal Keshava  & John Mustard, Spectral Unmixing, Signal
#' Processing Magazine, IEEE. 19. 44 - 57. 10.1109/79.974727, 2002. S.M. Kay, Fundamentals of Statistical Signal Processing: Estimation
#' Theory, Englewood Cliffs, NJ: Prentice Hall, 1993.
#
#'
#'
#'
#' @export
#'






estimateabundanceLS <- function(x, endmembers, method)
{
   if (class(x)[1] == "RasterStack" &&  class(x)[1] =="RasterBrick")
   {
      Ma <- as.matrix(x)
   } else if (class(x)[1] == "matrix") {
      Ma <- Rfast::transpose(x)
   } else {
      stop("x is not a matrix or raster")
   }


   if (class(endmembers)[1] != "matrix") {
      stop("endmembers need to be a matrix")
   }

   if (missing(method))
   {
      method <- "fcls"
   }
  
  #divide by 10000 to avoid overflow errors due to large numbers 
  
  if (mean(Ma,na.rm=TRUE)>10)
  {
  MA <- MA/10000
  endmembers <- endmembers/10000
  }

   # columns are bands
   # rows are signals
   N <- dim(Ma)[1]
   p <- dim(Ma)[2]

   #number of endmembers
  
   q <- dim(endmembers)[2]

   if(q >= p)
   {
   stop("Error: The numbers of endmembers equals o is greater than the number of bands")
   }
   abundanceMapTemp <- matrix(nrow = q, ncol = N, 0)

   if (method == "ucls")
   {
      sigInv <- pracma::pinv(endmembers)
      for (i in 1:N)
      {
         if (anyNA(Ma[i, ]))
         {
            abundanceMapTemp[, i] <- NA
         } else {
            abundanceMapTemp[, i] <- sigInv %*% Ma[i, ]
         }
      }

   } else if (method == "ncls") {
      for (i in 1:N)
      {
         if (anyNA(Ma[i, ]))
         {
            abundanceMapTemp[, i] <- NA
         } else {
            abundanceMapTemp[, i] <-
               pracma::lsqlincon(endmembers,
                         Ma[i,],
                         lb = 0,
                         ub = 1)
         }

      }


   } else if (method == "fcls") {
      for (i in 1:N)
      {
        
         if (anyNA(Ma[i,]))
         {
            abundanceMapTemp[, i] <- NA
         } else {
            abundanceMapTemp[, i] <-
               pracma::lsqlincon(
                  endmembers,
                  Ma[i, ],
                  lb = 0,
                  ub = 1,
                  beq = 1,
                  Aeq = matrix(ncol = q, nrow = 1, 1)
               )
         }
      }
   } else{
      stop("method is not valid")
   }



   if (class(x)[1] == "RasterStack" | class(x)[1] =="RasterBrick")
   {
      x <- x[[1:q]]

      values(x) <-  Rfast::transpose(abundanceMapTemp)

      names(x) <- rep("endmember", q)

   } else {
      x <- abundanceMapTemp
   }



   return(x)

}
