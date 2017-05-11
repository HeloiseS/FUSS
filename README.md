# For Use with Supernova Spectropolarimetry (FUSS)

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![DOI](https://zenodo.org/badge/79924113.svg)](https://zenodo.org/badge/latestdoi/79924113)

Our code For Use with Supernova Spectropolarimetry was used to reduce and analyse the data for our paper on SN 2014ad, see https://arxiv.org/pdf/1704.06270.pdf.

If you want to use our code, or part of our code, we would be grateful if you could cite the aforementioned publication.

For documentation, see our wiki, or the pdf in the docs/FUSS folder. 

## Installation
### Python
`pip install fussy`

### IRAF scripts

Create tasks in your login.cl file for the .cl scripts located in FUSS/iraf_scripts

e.g. `task    $ascii_f = [path_to_iraf_scripts]$toascii_flux.cl`

## Licensed

This project is Copyright (c) H. F. Stevance and licensed under the terms of the MIT license. See the licenses folder for more information.
