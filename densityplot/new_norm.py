from matplotlib.colors import Normalize,LinearSegmentedColormap,cbook,ma
from matplotlib.cm import datad
import numpy as np

"""These functions allow for the creation of 'divergent' colorbars
out of *any* two existing matplotlib colorbars.  They are
stitched together such that the split between the two colorbars
always happens at zero.

Usage
-----
from densityplot import *
from pylab import *
    
#create some data to plot
x = arange(0, pi, 0.1)
y = arange(0, 2*pi, 0.1)
X, Y = meshgrid(x,y)
Z = cos(X) * sin(Y) * 10

#stitch together two colorbars
#(the second colorbar is automatically reversed)
dub_cm=mk_dub_color('cubehelix_purple','cubehelix_green')
#set the normalization such that the split is at zero
n=MidNorm(vmax=10,vmin=-5)

#use this colorbar and normiazation to plot the image
imshow(Z, interpolation='nearest', cmap=dub_cm, norm=n, vmax=10, vmin=-15)
colorbar()
show()
"""

class MidNorm(Normalize):
    """A subclass on Normalize to map all pisitive values to
    the range [.5,1] and all negitve valeus to the range
    [0,.5).  This means 0 will always be mapped to the
    *middle* of a colorbar.

    Usage:
        from desnityplot import *
        #make a normialization that maps [-5,10] to
        #[0,1] such that 0 is mapped to 0.5
        n=MidNorm(vmax=10,vmin=-5)
    """
    def __init__(self,vmin=None,vmax=None,clip=False):
        Normalize.__init__(self,vmin,vmax,clip)

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > 0:
            raise ValueError("minvalue must be less than 0")
        if vmax < 0:
            raise ValueError("maxvalue must be more than 0")            
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            resdat[resdat>0] /= vmax
            resdat[resdat<0] /= -vmin
            resdat=resdat/2.+0.5
            result = np.ma.array(resdat, mask=result.mask, copy=False)
        if is_scalar:
            result = result[0]
        return result
    
    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax
        if cbook.iterable(value):
            val = ma.asarray(value)
            val=2*(val-0.5) 
            val[val>0]*=vmax
            val[val<0]*=-vmin
            return val
        else:
            if val<0.5: 
                return  2*val*(-vmin)
            else:
                return val*vmax

def mk_dub_color(cmHigh,cmLow):
    """
    This function takes two matplotlib colomaps and makes a single
    colorbar out of them.
    Parameters
    ----------
    cmHigh : name of mpl colormap (as string) to be used for
             high values (normalized values > 0.5)
             *note* this colorbar will be reversed
    cmLow :  name of mpl colormap (as string) to be used for
             low values (normalized valeus < 0.5)

    Return
    ------
    dub_color : matplotlib colormap
             
    Usage
    -----
    from densityplot import *
    from pylab import *
    dub_cm=mk_dub_color('cubehelix_purple','cubehelix_green')
    """
    cH=datad[cmHigh] #get the color dict for cmHigh
    cL=datad[cmLow] #get the color dict for cmLow
    norm_high=Normalize(vmin=0.5,vmax=1.0)
    norm_low=Normalize(vmin=0.0,vmax=0.5)
    def dub_color_get(cH,cL,key):
        def color(x):
            hdx=(x>=0.5)
            ldx=(x<0.5)
            out=np.zeros(len(x))
            if hdx.any():
                xn=norm_high(x[hdx])
                xn=abs(xn-1) #invert this color map so it goes white to black
                out[hdx]=cH[key](xn)
            if ldx.any():
                xn=norm_low(x[ldx])
                out[ldx]=cL[key](xn)
            return out    
        return color
    dub_color_dict={'red':dub_color_get(cH,cL,'red'),
                    'green':dub_color_get(cH,cL,'green'),
                    'blue':dub_color_get(cH,cL,'blue')}
    return LinearSegmentedColormap('dub_color',dub_color_dict,256)

if __name__=='__main__':
    from custom_cm import *
    from pylab import *
    
    #create some data to plot
    x = arange(0, pi, 0.1)
    y = arange(0, 2*pi, 0.1)
    X, Y = meshgrid(x,y)
    Z = cos(X) * sin(Y) * 10
    
    #stitch together two colorbars
    #(the second colorbar is automatically reversed)
    dub_cm=mk_dub_color('cubehelix_purple','cubehelix_green')
    #set the normalization such that the split is at zero
    n=MidNorm(vmax=10,vmin=-5)

    #use this colorbar and normiazation to plot the image
    imshow(Z, interpolation='nearest', cmap=dub_cm, norm=n, vmax=10, vmin=-15)
    colorbar()
    show()
