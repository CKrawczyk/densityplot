"""
Adds some new cubehelix color maps to matplotlib.cm
This file must be imported *before* matplotlib for it to work.
After importing, these color maps can be used in exactly the
same way as matplotlib's internal color maps.

Each new color map is a variant on the cubhelix color map
with no rotation in color space, each one passing thorugh
a different color. The _l versions have more color variation
for smaller values than larger values.

Usage:
    from densityplot import *
    from pylab import *

    new_colormap = cm.cubehelix_purple
    new_colormap_reverse = cm.cubehelix_purple_r
"""
from matplotlib import _cm as c

c.datad['cubehelix_light'] = c.cubehelix(gamma=0.6)
c.datad['cubehelix_purple'] = c.cubehelix(gamma=0.6, s=0.5, r=0.0, h=1.0)
c.datad['cubehelix_red'] = c.cubehelix(gamma=0.6, s=1.0, r=0.0, h=1.0)
c.datad['cubehelix_yellow'] = c.cubehelix(gamma=0.6, s=1.5, r=0.0, h=1.0)
c.datad['cubehelix_green'] = c.cubehelix(gamma=0.6, s=2.0, r=0.0, h=1.0)
c.datad['cubehelix_teal'] = c.cubehelix(gamma=0.6, s=2.5, r=0.0, h=1.0)
c.datad['cubehelix_blue'] = c.cubehelix(gamma=0.6, s=3.0, r=0.0, h=1.0)
c.datad['cubehelix_purple_l'] = c.cubehelix(gamma=0.4, s=0.5, r=0.0, h=1.5)
c.datad['cubehelix_red_l'] = c.cubehelix(gamma=0.4, s=1.0, r=0.0, h=1.5)
c.datad['cubehelix_yellow_l'] = c.cubehelix(gamma=0.4, s=1.5, r=0.0, h=1.5)
c.datad['cubehelix_green_l'] = c.cubehelix(gamma=0.4, s=2.0, r=0.0, h=1.5)
c.datad['cubehelix_teal_l'] = c.cubehelix(gamma=0.4, s=2.5, r=0.0, h=1.5)
c.datad['cubehelix_blue_l'] = c.cubehelix(gamma=0.4, s=3.0, r=0.0, h=1.5)
