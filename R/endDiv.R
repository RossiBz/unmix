

#' Subpixel spectral diversity
#'
#' Calculates different diversity indices for each pixel based on the endmember abundance map. 
#'
#' @param x RasterStack.  Each layer displays the abundance map of an endmember.
#'
#' @return Endmember diversity per pixel as RasterStack. The four layers include Richness, Shannon-Weaver index, Simpson index and Evenness 
#' calculated based on the endmember abundance map
#'
#' @author Christian Rossi christian.rossi1990@gmail.com
#'
#' @import raster
#' @import vegan
#'
#' @usage endDiv(x)
#'
#' @references Rossi, C., and Gholizadeh, H. (2022). Subpixel spectral diversity: Using endmember diversity to capture plant diversity.
#' submitted to Methods in Ecology and Evolution.
# 
#'
#'
#'
#' @export
#'  


endDiv <- function(x)
  
{
  if (class(x)[1] != "RasterStack")
  {
    stop("x is not RasterStack")
  }
  
  
  
  ric_endmembers <- raster::calc(
    x,
    fun = function(y) {
      #taxonomic diversity
      Richness_endm <- specnumber(y)
      Shannon_weaver_endm <- diversity(y)
      Simpson_endm <- diversity(y, index = "simpson")
      Eveness_endm <- Shannon_weaver_endm / log(Richness_endm)
      
      cbind(Richness_endm,
            Shannon_weaver_endm,
            Simpson_endm,
            Eveness_endm)
    }
  )
  
  return(ric_endmembers)
  
}