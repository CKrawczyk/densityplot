import matplotlib as mpl
import pylab as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def hex_difference(xy1,xy2,show_all=False,hkwargs={},color_bar=True,fignum=1):
    """A function that plots the difference between two hexbin plots.

    Parameters
    ----------
    xy1 : A tuple of (x,y) corrdianates for the first hexbin. A tuple of (x,y,C)
          can also be passed in where C is the value for each (x,y) point.
    xy2 : A tuple of (x,y) corrdianates for the second hexbin. A tuple of (x,y,C)
          can also be passed in where C is the value for each (x,y) point.

    NOTE : the 'C' functinality is untested and may not work as expected.

    Keywords
    --------
    show_all  : bool (optional)
                If True all intermediate hexbin plots are returned.
                Default: show_all=False
    color_bar : bool (optional)
                If True a colorbar is placed on the plot(s)
                Default: colorbar=True
    fignum    : int (optional)
                The number to give the resulting figure(s).  If 
                show_all=True, the intermediate plots will be
                fignum+1 and fignum+2 while the difference will
                be fignum.
                default: fignum=1
               
    Passed Keywords
    ---------------
    hkwargs:
        a dictionary of keywords passed to hexbin
        (see http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hexbin
        for additional keywords that can be set)

    Returns
    -------
    d  : pylab.hexbin object
         Object returned by the difference hexbin
    h1 : pylab.hexbin object
         Object returned by hexbin of first data set
         NOTE: only returned if show_all=True
    h2 : pylab.hexbin object
         Object returned by hexbin of second data set
         NOTE: only returned if show_all=True
    c  : matplotlib.colorbar.Colorbar instance
         NOTE: only returned if color_bar=True

    Usage
    -----
    import numpy as np
    n=100000
    x1=np.random.standard_normal(n) #random x points
    y1=2+3*x1+4*np.random.standard_normal(n) #random y points
    
    x2=np.random.standard_normal(n) #random x points
    y2=2-3*x2+4*np.random.standard_normal(n) #random y points
    
    hex_difference((x1,y1),(x2,y2),show_all=True,color_bar=True,hkwargs={'gridsize':100,'extent':[-4.5,4.5,-25,25],'vmin':-180,'vmax':180})
    pl.show()

    """
    if show_all: #if shoing all hexbins then draw them as you go (you can't change the drawing axis object after creation)
        pl.figure(fignum+1)
        hex1=pl.hexbin(*xy1,**hkwargs)
        if color_bar:
            pl.colorbar()
        pl.figure(fignum+2)
        hex2=pl.hexbin(*xy2,**hkwargs)
        if color_bar:
            pl.colorbar()
    else: #make but don't draw the 2 hexbins
        hex1=pl.hexbin(*xy1,visible=False,**hkwargs) #make the hexbins, visible it False to avoid drawing them to a plot
        hex2=pl.hexbin(*xy2,visible=False,**hkwargs)
    pl.figure(fignum)
    hex_dif=pl.hexbin(*xy1,visible=False,**hkwargs) #this will have the counts overwritten (so don't draw yet)
    c1=hex1.get_array() #the counts for hex1
    c2=hex2.get_array() #the counts for hex2
    c_dif=c1-c2 #difference between plots
    gdx=~((c1==0)&(c2==0)) #the bins to draw (removes where both hists had no counts)
    #NOTE: if the 'C' values are set checking against 0 is NOT a good idea...
    hex_dif.set_array(c_dif[gdx]) #set the defferences into the hex_dif object

    h=hex_dif.get_paths() #get the hexagon Path object(s)
    if len(h)>1: #you have an old version of matplotlib, use this bit of code
        rem_me=pl.array(h)[~gdx] #bins to remove
        for r in rem_me:
            h.remove(r) #remove blank bins
    else: #either you have a boaring hexbin or a newer version of matplotlib
        h=hex_dif.get_offsets()
        hex_dif.set_offsets(h[gdx])
    hex_dif.set_visible(True) #this draws the new hex_dif
    ret=[hex_dif]
    if show_all:
        ret.append(hex1)
        ret.append(hex2)
    if color_bar:
        ains=inset_axes(pl.gca(),width='80%',height='5%',loc=9)
        #TODO: externalize colorbar keywords
        c=pl.colorbar(hex_dif,cax=ains,orientation='horizontal')
        ret.append(c)
    return tuple(ret)

if __name__=='__main__':
    import numpy as np
    n=100000
    x1=np.random.standard_normal(n) #random x points
    y1=2+3*x1+4*np.random.standard_normal(n) #random y points
    
    x2=np.random.standard_normal(n) #random x points
    y2=2-3*x2+4*np.random.standard_normal(n) #random y points
    
    hex_difference((x1,y1),(x2,y2),show_all=True,color_bar=True,hkwargs={'gridsize':100,'extent':[-4.5,4.5,-25,25],'vmin':-180,'vmax':180})
    pl.show()
