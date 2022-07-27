# unmix


## An R package for unmixing of multi- and hyperspectral images

The packages include state of the art endmember extraction algortihms like the Pixel Purity Index (PPI) and the Vertex Component Analysis (VCA). 
In addition it allows to find the number of endmembers present in a hyperspectral data cube by using Harsanyi–Farrand–Chang (NWHFC) or HySime method.
Furthermore, spectral unmixing per pixel can be performed,  which is the decomposition of the spectra of each pixel into a given set of endmember spectra.
Finally, the endmember diversity, i.e. subpixel spectral diveristy, per pixel can be calculated based on different diversity indices (i.e., Richness, Shannon-Weaver, Simpson and Evenness) 
as proposed by Rossi and Gholizadeh (2022).


## Main features

* Virtual dimensionality HySime and NWHFC
* PPI and VCA
* Linear unmixing with three options (ucls — Unconstrained least-squares method, fcls — Fully constrained least-squares method, 
ncls — Nonnegative constrained least-squares method)
* Subpixel spectral diversity from the endmember diversity


## Installation

To load (using `devtools`):

```Rscript
library(devtools)
install_github("RossiBz/unmix")
```


## Authors

* Christian Rossi

## License

Licensed under the GNU General Public License, Version 3.0: https://www.gnu.org/licenses/gpl-3.0.html

## Reference

Rossi, C., and Gholizadeh, H. (2022). Subpixel spectral diversity: Using endmember diversity to capture plant diversity.
submitted to Methods in Ecology and Evolution.

