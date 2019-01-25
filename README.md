# For Use with Supernova Spectropolarimetry (FUSS)

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![DOI](https://zenodo.org/badge/79924113.svg)](https://zenodo.org/badge/latestdoi/79924113)

Our code For Use with Supernova Spectropolarimetry was used to reduce and analyse the data for our papers:
 - [The evolution of the 3D shape of the broad-lined Type Ic SN 2014ad](http://adsabs.harvard.edu/abs/2017MNRAS.469.1897S)
 - [Probing the rotational velocity of Galactic WO stars with spectropolarimetry](http://adsabs.harvard.edu/abs/2018MNRAS.479.4535S)
 - [The 3D shape of Type IIb SN 2011hs](http://adsabs.harvard.edu/abs/2019arXiv190107562S)

If you want to use our code, or part of our code, we would be grateful if you could cite the aforementioned publication.

For documentation, see our wiki, or the pdf in the docs/FUSS folder. 

## Installation
### Python
`pip install fusspy`

### IRAF scripts

Create tasks in your login.cl file for the .cl scripts (not the python scripts!) located in FUSS/iraf_scripts

e.g. `task    $ascii_f = [path_to_iraf_scripts]$toascii_flux.cl`

You will also have to edit the path to the python scripts in `arc_apertures.cl`, `create_id.cl`, `toascii_flux.cl` and `toascii_pol.cl` manually to fit the location of your iraf_scripts folder.

## Licensed

This project is Copyright (c) H. F. Stevance and licensed under the terms of the MIT license. See the licenses folder for more information.

