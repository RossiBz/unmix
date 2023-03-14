

#' Subpixel spectral diversity
#'
#' Calculates different diversity indices for each pixel based on the endmember abundance map.
#'
#' @param x RasterStack, RasterBrick.  Each layer displays the abundance map of an endmember.
#'
#' @return Endmember diversity per pixel as RasterBrick. The four layers include Richness, Shannon-Weaver index, Simpson index and Evenness
#' calculated based on the endmember abundance map
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @rawNamespace import(vegan, except = procrustes)
#'
#' @usage endDiv(x)
#'
#' @references Rossi, C., and Gholizadeh, H. (2023). Uncovering the hidden diversity: Leveraging sub-pixel spectral diversity to estimate plant diversity from space.
#' submitted to Remote Sensing of Environment
#
#'
#'
#'
#' @export
#'


endDiv <- function(x)

{
  if (class(x)[1] != "RasterStack" & class(x)[1] !="RasterBrick")
  {
    stop("x is not RasterStack or RasterBrick")
  }



  ric_endmembers <- raster::calc(
    x,
    fun = function(y) {
      #taxonomic diversity
      Richness_endm <- vegan::specnumber(y)
      Shannon_weaver_endm <- vegan::diversity(y)
      Simpson_endm <- vegan::diversity(y, index = "invsimpson")
      Eveness_endm <- Shannon_weaver_endm / log(Richness_endm)

      cbind(Richness_endm,
            Shannon_weaver_endm,
            Simpson_endm,
            Eveness_endm)
    }
  )

  return(ric_endmembers)

}
