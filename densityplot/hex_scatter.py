import matplotlib as mpl
import pylab as pl
import numpy as np

# Plot a hexbin object with outliers plotted as points:
# Useage: hex_scatter(x = [numpy array], y = [numpy array], min_cnt = [2 | integer], levels = [2 | integer], std = [False | True], smoothing = [3 | integer], hkwargs = {},skwargs = {},ckwargs = {})
#
# Parameters
#----------
# x : array_like
#   the x coordinates for the input data
# y : array_like
#   the x coordinates for the input data
# min_cnt: integer (optional)
#   minimum number of hex bin counts to draw contour, less than or equal to this will be points
# levels: integer or array (optional)
#   the level curves to draw, if integer will draw that number of curves Default = 2
#   if std = True, level contours will be drawn at 1, 2, 3, 4, and 5 sigma up to the number of levels.
#   if array, levels will be drawn at those values.
# std : bool (optional)
#   make level contours of the fraction total number of points within the level curve Default = False
#   if levels is an integer level contours will be drawn at 1, 2, 3, 4, and 5 sigma up to the number of levels.
# smoothing: integer (optional)
# See: http://matplotlib.org/1.3.0/api/tri_api.html?highlight = uniformtrirefiner#matplotlib.tri.UniformTriRefiner
#   Default = 3, values from 1 to 4 are recommended.
# additional keyword arguments are passed in hkwargs = {},skwargs = {},ckwargs = {}
# hkwargs:
#   a dictionary of keywords passed to hexbin
# skwargs:
#   a dictionary of keywords passed to plot (for scatter points)
# ckwargs:
#   a dictionary of keywords passed to tricontour (if std is set)
#
# Returns
#-------
# f : pylab.hexbin object
#   hexbin plotted by the routine
# P : pylab.plot object
#   outlier points plotted by the routine
# T : pylab.tricontour object
#   T = None, is std = False
#   the contour created by the routine
#
#
#
#
# Plot a contour with outliers plotted as points.
#
# Useage: hex_contour(x = [numpy array], y = [numpy array], min_cnt = [2 | integer],levels = [2 | integer], std = [False | True], smoothing = [3 | integer], hkwargs = {},skwargs = {},ckwargs = {})
#
# Parameters
#----------
# x : array_like
#   the x coordinates for the input data
# y : array_like
#   the x coordinates for the input data
# min_cnt: integer (optional)
#   minimum number of hex bin counts to draw contour, less than or equal to this will be points
# levels: integer or array (optional)
#   the level curves to draw, if integer will draw that number of curves Default = 2
#   if std = True, level contours will be drawn at 1, 2, 3, 4, and 5 sigma up to the number of levels.
#   if array, levels will be drawn at those values.
# std : bool (optional)
#   make level contours of the fraction total number of points within the level curve Default = False
#   if levels is an integer level contours will be drawn at 1, 2, 3, 4, and 5 sigma up to the number of levels.
# smoothing: integer (optional)
# See: http://matplotlib.org/1.3.0/api/tri_api.html?highlight = uniformtrirefiner#matplotlib.tri.UniformTriRefiner
#   Default = 3, values from 1 to 4 are recommended.
# additional keyword arguments are passed in hkwargs = {},skwargs = {},ckwargs = {}
# hkwargs:
#   a dictionary of keywords passed to hexbin
# skwargs:
#   a dictionary of keywords passed to plot (for scatter points)
# ckwargs:
#   a dictionary of keywords passed to tricontour (if std is set)
#
# Returns
#-------
# T : pylab.tricontour object
#   the contour created by the routine
# P : pylab.plot object
#   outlier points plotted by the routine
#
#
#
# Common keywords used in dictionaries for both functions
#
# hkwargs (See http://matplotlib.org/1.3.1/api/pyplot_api.html?highlight = hexbin#matplotlib.pyplot.hexbin)
# gridsize: [ 100 | integer ]
#   The number of hexagons in the x-direction, default is 100.
#   The corresponding number of hexagons in the y-direction is chosen such that the hexagons are approximately regular.
#   Alternatively, gridsize can be a tuple with two elements specifying the number of hexagons in the x-direction and the y-direction.
#
# bins: [ None | 'log' | integer | sequence ]
#   If None, no binning is applied; the color of each hexagon directly corresponds to its count value.
#   If 'log', use a logarithmic scale for the color map. Internally,  is used to determine the hexagon color.
#   If an integer, divide the counts in the specified number of bins, and color the hexagons accordingly.
#   If a sequence of values, the values of the lower bound of the bins to be used.
#
# xscale: [ 'linear' | 'log' ]
#   Use a linear or log10 scale on the horizontal axis.
#
# scale: [ 'linear' | 'log' ]
#   Use a linear or log10 scale on the vertical axis.
#
# mincnt: [ None | a positive integer ]
#   If not None, only display cells with more than mincnt number of points in the cell.
#
# marginals: [ True | False ]
#   if marginals is True, plot the marginal density as colormapped rectagles along the bottom of the x-axis and left of the y-axis.
#
# extent: [ None | scalars (left, right, bottom, top) ]
#   The limits of the bins. The default assigns the limits based on gridsize, x, y, xscale, and yscale.
#   Changing the default extent does not work with log scales.
#
# cmap: [ None | Colormap ]
#   a matplotlib.colors.Colormap instance. If None, defaults to rc image.cmap.
#
# norm: [ None | Normalize ]
#   matplotlib.colors.Normalize instance is used to scale luminance data to 0,1.
#
# vmin / vmax: scalar
#   vmin and vmax are used in conjunction with norm to normalize luminance data. If either are None, the min and max of the color array C is used. Note if you pass a norm instance, your settings for vmin and vmax will be ignored.
#
# alpha: scalar between 0 and 1, or None
#   the alpha value for the patches
#
# linewidths: [ None | scalar ]
#   If None, defaults to rc lines.linewidth. Note that this is a tuple, and if you set the linewidths argument you must set it as a sequence of floats, as required by RegularPolyCollection.
#   Other keyword arguments controlling the Collection properties:
#
# edgecolors: [ None | 'none' | mpl color | color sequence ]
#   If 'none', draws the edges in the same color as the fill color. This is the default, as it avoids unsightly unpainted pixels between the hexagons.
#   If None, draws the outlines in the default color.
#   If a matplotlib color arg or sequence of rgba tuples, draws the outlines in the specified color.
#
#
#
# skwargs (See: http://matplotlib.org/api/pyplot_api.html?highlight = plot#matplotlib.pyplot.plot)
#
# alpha: float (0.0 transparent through 1.0 opaque)
#
# color or c: any matplotlib color
#
# label:string or anything printable with '%s' conversion.
#
# marker: ['.' | ',' | 'o'] (See http://matplotlib.org/api/markers_api.html#module-matplotlib.markers for more)
#
# markeredgecolor or mec: any matplotlib color
#
# markeredgewidth or mew: float value in points
#
# markerfacecolor or mfc: any matplotlib color
#
# markerfacecoloralt or mfcalt: any matplotlib color
#
# markersize or ms: [2 | float]
#
#
#
# ckwargs (See: http://matplotlib.org/1.3.1/api/pyplot_api.html?highlight = tricontour#matplotlib.pyplot.tricontour)
#
# colors: [ None | string | (mpl_colors) ]
#   If None, the colormap specified by cmap will be used.
#   If a string, like 'r' or 'red', all levels will be plotted in this color.
#   If a tuple of matplotlib color args (string, float, rgb, etc), different levels will be plotted in different colors in the order specified.
#
# alpha: float
#   The alpha blending value
#
# cmap: [ None | Colormap ]
#   A cm Colormap instance or None. If cmap is None and colors is None, a default Colormap is used.
#
# norm: [ None | Normalize ]
#   A matplotlib.colors.Normalize instance for scaling data values to colors. If norm is None and colors is None, the default linear scaling is used.
#
#levels [level0, level1, ..., leveln]
# A list of floating point numbers indicating the level curves to draw; eg
# to draw just the zero contour pass levels = [0]

try:
    from matplotlib.tri import Triangulation, UniformTriRefiner
except:
    print 'Note: You need matplotlib 1.3.0 or higher to use hex_contour or standard deviation contours on hex_scatter!'


def hex_scatter(x, y, min_cnt=2, std=False, levels=2, level_at_edge=False, smoothing=3,
                hkwargs={}, skwargs={}, ckwargs={}):
    """Plot a hexbin object with outliers plotted as points.

    Parameters
    ----------
    x : array_like
        the x coordinates for the input data
    y : array_like
        the y coordinates for the input data

    Keywords
    --------
    min_cnt: integer (optional)
        Minimum number of hex bin counts to draw contour. Bins
        containing less than or equal to this will be scatter points.
        Default: min_cnt = 2
    std : bool (optional)
        Draw contours of the fraction of points contained within each
        level. The levels are set using the levels keword (see below).
        Default: std = False
    levels: integer or array (optional)
        Only available if std = True, and indicates either the number of
        contour levels to draw (integer), or an array of levels to draw.
        If integer, levels are drawn at [1, 2, 3, 4, 5][:levels] simga,
        i.e. levels=2 will draw contours at 1 and 2 sigma. If array,
        the values indicate the fraction of points contained within each
        contour.  Default: levels = 2
    level_at_edge: bool (optional)
        Only available if std = True, and indicates if the the scatter
        points should start at the last contour level given.
        Dfault: cont_at_edge = False
    smoothing: integer (optional)
        Only available if std = True, and is used to subdivide (i.e.
        smooth) the contouring grid. Default: smoothing = 3, values from
        1 to 4 are recommended.
        See subdiv on: http://matplotlib.org/1.3.0/api/tri_api.html?highlight=uniformtrirefiner#matplotlib.tri.UniformTriRefiner
        

    Passed Keywords
    ---------------
    hkwargs:
        a dictionary of keywords passed to hexbin
    skwargs:
        a dictionary of keywords passed to plot (for scatter points)
    ckwargs:
        a dictionary of keywords passed to tricontour (if std is set)

    Returns
    -------
    f: pylab.hexbin object
        Object returned by hexbin
    P: pylab.plot object
        Object returned by plot
    T: pylab.tricontour object
        if std = True, object returned by tricontour
        if std = False, None is returned
    """

    # check if you are adding curves on top of the hex bin, if so reverse the
    # default colorbar. (lower contours have higher levels)
    if std:  
        colorbar = hkwargs.pop('cmap', mpl.cm.winter_r)
    else:
        colorbar = hkwargs.pop('cmap', mpl.cm.winter)
    # set the default scales on x and y axis to linear
    hkwargs.setdefault('xscale', 'linear')
    hkwargs.setdefault('yscale', 'linear')
    # hexbin the data to a mincnt of 0 (needed for sig_cont to look good),
    # visible = False so it does not draw yet
    f = pl.hexbin(x, y, mincnt=0, visible=False, cmap=colorbar, **hkwargs)
    # no draw commands are sent until visible = True
    c = f.get_array()  # get the number of points in each bin
    idx = (c <= min_cnt)  # mask for all bins with less than the min_cnt
    if std:  # check if you are adding curves on top of the hex bin
        # if you just indicate the number of contours you want
        if isinstance(levels,int):
            # dfault to [1, 2, 3, 4, 5][:levels] sigma contorus
            ckwargs['levels'] = [0.682689492, 0.954499736, 0.997300204,
                                 0.99993666, 0.999999426697][0:levels]
        else:
            ckwargs['levels'] = levels
        # convert the hexbin counts to fraction of the total number of points
        c1 = convert_to_stdev(c)
        cmin_std = max(ckwargs['levels'])
        # mask for all bins with less than the min_cnt 
        idx = (c<=min_cnt)
        # mask for all bins less than cmin_std
        if level_at_edge:
            idx = idx | (c1 > cmin_std)
        f.set_array(c1[~idx])  # keep values above cutoff
        vmin=hkwargs.get('vmin',0.0)
        vmax=hkwargs.get('vmax',1.0)
        f.set_clim(vmin=vmin, vmax=vmax)  # set the limits on the colorbar
    else:  # if just plotting the hexbins
        f.set_array(c[~idx])  # keep values above cutoff

    H = f.get_paths()  # get the hexagon Path objects
    xout = pl.array([])  # array to hold x scatter points
    yout = pl.array([])  # array to hold y scatter points

    # either this is a very boring hexbin or you have the newer version of
    # matplotlib and binning on a linear scale (log scale still uses the old
    # hexbin convention)
    if (len(H) == 1) \
        and (hkwargs['xscale'] == 'linear') \
        and (hkwargs['yscale'] == 'linear'):
        ext = H[0].get_extents()  # the extents for this one hex
        h = f.get_offsets()  # the centers for each hex
        if std:  # check if you are adding curves on top of the hexbin
            # plot contours
            T = triangle_contour(h.T[0], h.T[1], c1, smoothing, ckwargs)
        out = h[idx]  # only the hexagons with less than min_cnt
        # keep only the hexagons with values larger than min_cnt
        f.set_offsets(h[~idx])
        idx = idx & (c > 0)
        # loop through all the hexagons with less than the min_cnt
        for o in out:
            # move data to current bin's center
            xnow = x - o[0]  
            ynow = y - o[1]
            # get the index of the points in extent of the hexagon
            jdx = (xnow >= ext.xmin) & (xnow <= ext.xmax) & \
                  (ynow >= ext.ymin) & (ynow <= ext.ymax)  
            # loop through each point in extent of the hexagon
            for p in zip(xnow[jdx], ynow[jdx]):
                # if the point is actual in the hexagon add it to the list of
                # outliers
                if H[0].contains_point(p):
                    xout = np.append(xout, p[0] + o[0])
                    yout = np.append(yout, p[1] + o[1])
    # you have the old version of matplotlib or non linear scaling on an axis
    else:
        if std:  # check if you are adding curves on top of the hex bin
            # get the centers of each hexagon
            xs = np.array([(0.5 * (i.get_extents().xmax - i.get_extents().xmin) \
                            + i.get_extents().xmin) for i in H])
            ys = np.array([(0.5 * (i.get_extents().ymax - i.get_extents().ymin) \
                            + i.get_extents().ymin) for i in H])
            # plot contours
            T = triangle_contour(xs, ys, c1, smoothing, ckwargs)
        out = pl.array(H)[idx]  # hexagons to be cut
        for o in out:  # loop over Paths
            H.remove(o)  # remove hexagons below cutoff
            ext = o.get_extents()  # get the x,y extent of the hexagon
            # mask for data points only in the hexagon extent
            jdx = (x >= ext.xmin) & (x <= ext.xmax) & \
                  (y >= ext.ymin) & (y <= ext.ymax) 
            for p in zip(x[jdx], y[jdx]):  # loop over these points
                if o.contains_point(p):  # check if point is in the hexagon
                    xout = np.append(xout, p[0])  # if so append the x value
                    yout = np.append(yout, p[1])  # if so append the y value
    f.set_visible(True)  # draw remaining hexagons
    # set default type and size for outlier points
    skwargs.setdefault('marker', '.')
    skwargs.setdefault('ms', 2)
    skwargs.setdefault('ls', 'None')
    P = pl.plot(xout, yout, **skwargs)  # plot the outlier points
    if std:
        return f, P, T
    else:
        return f, P, None


def triangle_contour(x_center, y_center, values, smoothing, ckwargs={}):
    """A helper funciton used to draw contorus with pl.tricontour"""
    # make Triangulation object using the centers of each of the hexbins
    triag = Triangulation(x_center, y_center)
    refiner = UniformTriRefiner(triag) # refines the mesh of triangle
    # returns refines triangle field of triangles and interpolated
    # contour values by dividing each triangle into 4**subdiv triangles
    tri_refi, c_refi = refiner.refine_field(values, subdiv=smoothing)
    T = pl.tricontour(tri_refi, c_refi, **ckwargs)
    return T
    
def hex_contour(x, y, min_cnt=2, std=False, levels=2, smoothing=3,
        hkwargs={}, skwargs={}, ckwargs={}):
    """Plot density contours with outliers plotted as scatter points.
                
    Parameters
    ----------
    x : array_like
        the x coordinates for the input data
    y : array_like
        the y coordinates for the input data

    Keywords
    --------
    min_cnt: integer (optional)
        Minimum number of hex bin counts to draw contour. Bins
        containing less than or equal to this will be scatter points.
        Default: min_cnt = 2
    std : bool (optional)
        Draw contours of the fraction of points contained within each
        level. The levels are set using the levels keword (see below).
        Default: std = False
    levels: integer or array (optional)
        Indicates either the number of contour levels to draw (integer),
        or an array of levels to draw. If integer and std = True, levels
        are drawn at [1, 2, 3, 4, 5][:levels] simga, otherwise the
        levels will be evenly spaced in density from min_cnt upto the
        maximum bin count. If array and std = True, the values indicate
        the fraction of points contained within each contour, otherwise
        they indicate the number of points in each bin.
        Default: levels = 2
    smoothing: integer (optional)
        Is used to subdivide (i.e. smooth) the contouring grid.
        Default: smoothing = 3, values from 1 to 4 are recommended.
        See subdiv on: http://matplotlib.org/1.3.0/api/tri_api.html?highlight=uniformtrirefiner#matplotlib.tri.UniformTriRefiner
        

    Passed Keywords
    ---------------
    hkwargs:
        a dictionary of keywords passed to hexbin
    skwargs:
        a dictionary of keywords passed to plot (for scatter points)
    ckwargs:
        a dictionary of keywords passed to tricontour (if std is set)

    Returns
    -------
    f: pylab.hexbin object
        Object returned by hexbin
    P: pylab.plot object
        Object returned by plot
    T: pylab.tricontour object
        Object returned by tricontour
    """
    # set the default scales on x and y axis to linear
    hkwargs.setdefault('xscale', 'linear')  
    hkwargs.setdefault('yscale', 'linear')
    # hexbin the data to a mincnt of 0 (needed for sig_cont to look good),
    # visible = False so it does not draw
    f = pl.hexbin(x, y, mincnt=1, visible=False, **hkwargs)
    # NOTE: may have to change the mincnt if the C keyword is set...
    c = f.get_array()  # get the number of points in each bin
    H = f.get_paths()  # get the hexagon Path objects
    # if binning is linear hexbin uses the new convention
    if (hkwargs['xscale'] == 'linear') and (hkwargs['yscale'] == 'linear'):
        h = f.get_offsets()  # the centers for each hex
        xs = h.T[0]
        ys = h.T[1]
    # if bining is non-linear it uses the old convention
    else:
        # get the centers of each hexagon
        xs = np.array([(0.5 * (i.get_extents().xmax - i.get_extents().xmin) + \
                        i.get_extents().xmin) for i in H])  
        ys = np.array([(0.5 * (i.get_extents().ymax - i.get_extents().ymin) + \
                        i.get_extents().ymin) for i in H])
    if std:  # check if you are making standard deviation curves
        # convert the hexbin counts to fraction of the total number of points
        c1 = convert_to_stdev(c)
        # if you just indicate the number of contours you want
        if isinstance(levels, int):
            ckwargs['levels'] = [0.682689492, 0.954499736, 0.997300204,
                                 0.99993666, 0.999999426697][:levels]
        else:
            ckwargs['levels'] = levels
        # set where the last contour goes and where outlier points begin
        cmin_std = max(ckwargs['levels'])
        # mask for all bins with less than the min_cnt or less than cmin_std
        idx = ((c1 > cmin_std) | (c <= min_cnt)) & (c > 0)
        # make Triangulation object using the centers of each of the hexbins
        # plot contours
        T = triangle_contour(xs, ys, c1, smoothing, ckwargs)
    else:  # if you are making normal contours
        if isinstance(levels, int):
            # default to contours from min_cnt to just below the maximum value
            # in even steps.
            ckwargs['levels'] = np.linspace(min_cnt, c.max(), levels + 1)[:-1]
        else:
            ckwargs['levels'] = levels
        idx = (c <= min_cnt) & (c > 0)  # mask for all bins with less than the min_cnt
        # plot contours
        T = triangle_contour(xs, ys, c, smoothing, ckwargs)
    xout = pl.array([])  # array to hold x scatter points
    yout = pl.array([])  # array to hold y scatter points
    # if binning is linear hexbin uses the new convention
    if (hkwargs['xscale'] == 'linear') and (hkwargs['yscale'] == 'linear'):
        out = h[idx]  # only the hexagons with less than min_cnt
        H = H[0]  # there is only hexagon one path now centered at 0,0
        ext = H.get_extents()  # the extents of this one hexagon
        for o in out:
            # move data to current bin's center
            xnow = x - o[0]  
            ynow = y - o[1]
            # get the index of the points in extent of the hexagon
            jdx = (xnow >= ext.xmin) & (xnow <= ext.xmax) & \
                  (ynow >= ext.ymin) & (ynow <= ext.ymax)  
            # loop through each point in extent of the hexagon
            for p in zip(xnow[jdx], ynow[jdx]):
                # if the point is actual in the hexagon add it to the list of
                # outliers
                if H.contains_point(p):
                    xout = np.append(xout, p[0] + o[0])
                    yout = np.append(yout, p[1] + o[1])
    # if binning is non-linear it uses the old convention
    else:
        out = pl.array(H)[idx]  # hexagons to be cut
        for o in out:  # loop over Paths
            H.remove(o)  # remove hexagons below cutoff
            ext = o.get_extents()  # get the x,y extent of the hexagon
            # mask for data points only in the hexagon extent
            jdx = (x >= ext.xmin) & (x <= ext.xmax) & \
                  (y >= ext.ymin) & (y <= ext.ymax) 
            for p in zip(x[jdx], y[jdx]):  # loop over these points
                if o.contains_point(p):  # check if point is in the hexagon
                    xout = np.append(xout, p[0])  # if so append the x value
                    yout = np.append(yout, p[1])  # if so append the y value
    # set default type and size for outlier points
    skwargs.setdefault('marker', '.')
    skwargs.setdefault('ms', 2)
    skwargs.setdefault('ls', 'None')
    P = pl.plot(xout, yout, **skwargs)  # plot the outlier points
    return T, P


def convert_to_stdev(grid):
    """
    Based on http://pypi.python.org/simple/astroML/
    Given a grid of values, convert them to cumulative standard deviation.
    This is useful for drawing contours with stand deviation values as the levels.
    """
    shape = grid.shape
    grid = grid.ravel()
    
    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(grid)[::-1]
    i_unsort = np.argsort(i_sort)
    grid_cumsum = grid[i_sort].cumsum()
    grid_cumsum /=  grid_cumsum[-1]
    return grid_cumsum[i_unsort].reshape(shape)

if __name__ == '__main__':
    print """
    What example would you like to see?
    1: hex_scatter and hex_contour for a 2D normal distribution
    2: hex_contour for stelar color-color plot (takes ~3 mins)
    3: both
    """
    example_number = raw_input("example number: ")

    try:
        example_number = int(example_number)
        assert example_number in [1,2,3], "You must enter 1, 2, or 3!"
    except:
        raise Exception("You must enter 1, 2, or 3!")

    if (example_number == 1) or (example_number ==3):
        pl.np.random.seed(0)
        n = 100000
        x = pl.np.random.standard_normal(n)
        y = 2 + 3 * x + 4 * pl.np.random.standard_normal(n)

        pl.figure(1)
        pl.subplot(121)
        hex_scatter(x, y, min_cnt=2, levels=3, std=True, smoothing=1,
            hkwargs={'gridsize': 20}, skwargs={'color': 'k'})
        pl.subplot(122)
        hex_contour(x, y, min_cnt=2, levels=3, std=True, smoothing=1,
            hkwargs={'gridsize': 20}, skwargs={'color': 'k'})
    if (example_number == 2) or (example_number ==3):
        pl.figure(2)
        d = np.loadtxt('test.out', dtype='float')
        hex_contour(d[0], d[1], levels=[0.95, 0.9, 0.8, 0.6, 0.4, 0.2],
                    std=True, min_cnt=1, smoothing=3,
                    hkwargs={'gridsize': 100, 'xscale': 'linear'},
                    skwargs={'color': 'red', 'marker': ','},
                    ckwargs={'colors': 'red', 'linewidths': 0.7})
    pl.show()
