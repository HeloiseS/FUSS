"""
9 - May - 2017 / H. F. Stevance fstevance1@sheffield.ac.uk

Available modules:
 - datred
 - interactive_graph
 - isp 
 - polmisc (loaded here for backwards compatibility)
 - polplot
 - stat 
"""
from ._astropy_init import *

if not _ASTROPY_SETUP_:
    from polmisc import *


