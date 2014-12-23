from matplotlib.patches import Ellipse
import scipy.stats as st
import scipy.linalg as sl
from scipy import dot,zeros,array,rad2deg,diag,newaxis,arctan2,sqrt

class Ellipse_proj():
    """A class to find the n-dim ellipsoid containing a percentage p of the data points
    and porjects the ellipsoid down to lower dimensions. This is useful for finding
    confidence regions in n-dim data and then projecting it onto a series of 2d plots.

    NOTE: this makes the assumption of a n-dim normal distribution
    """
    def __init__(self,data,cov_matrix=False,loc=None):
        """Parameters
        ----------
        data : array of data, shape=(number points,number dim)
               If cov_matrix is True then data is the covariance matrix (see below)

        Keywords
        --------
        cov_matrix : bool (optional)
                     If True data is treated as a covariance matrix with shape=(number dim, number dim)
        loc        : the mean of the data if a covarinace matrix is given, shape=(number dim)
        """
        if cov_matrix:
            self.dim=data.shape[0]
            self.n=None
            self.data_t=None
            self.mu=loc
            self.evec,eval,V=sl.svd(data,full_matrices=False)
            self.sigma=sqrt(eval)
            self.Sigma=diag(1./self.sigma)
            self.B=dot(self.evec,self.Sigma)
            self.Binv=sl.inv(self.B)
        else:
            self.n,self.dim=data.shape #the shape of input data
            self.mu=data.mean(axis=0) #the mean of the data
            self.data_t=data-self.mu #remove the mean
            self.evec,eval,V=sl.svd(self.data_t.T,full_matrices=False) #get the eigenvectors (axes of the ellipsoid)
            data_p=dot(self.data_t,self.evec) #project the data onto the eigenvectors
            self.sigma=data_p.std(axis=0) #get the spread of the distribution (the axis ratos for the ellipsoid)
            self.Sigma=diag(1./self.sigma) #the eigenvalue matrix for the ellipsoid equation
            self.B=dot(self.evec,self.Sigma) #used in the ellipsoid equation
            self.Binv=sl.inv(self.B) #also useful to have around
    def proj(self,proj_vars):
        """This function will project the ellipsoid down onto a lower dimension space
        Parameters
        ----------
        proj_vars : array that is 1 for the projection dimension, and 0 otherwise.
                    i.e. array([0,0,1,0,1]) will project a 5d ellipsoid onto the plane
                    span by the 3rd and 5th variable.

        Return
        ------
        mu : The center of the ellipse projected into the lower dim space
        U  : The eigenvector for the projected ellipsoid axes (only returned if dim>1)
        S  : The non-scaled extent of the ellipsoid projected into the lower dim space
             NOTE: All scaling is done later to captuer 'p' percent of points inside
             the ellipsoid
        
        """
        pdim=proj_vars.sum() #the dim being projected into
        T=zeros([self.dim,pdim]) #the projection matrix
        #TODO: let the user specify the T matrix so they can do general projections
        for idx,i in enumerate(proj_vars.nonzero()[0]): #make a simple projection matrix
                T[i,idx]=1.
        mu_proj=dot(T.T,self.mu[:,newaxis]).flatten() #project the mean into the new space
        if pdim==1: #project ellipsoid down to a line
            S=sl.norm(dot(self.Binv,T),ord=2) #and this is why we want the B matrix...
            return mu_proj[0],S #return the mean and spread (NOTE: S is not scaled yet, that happens later)
        elif pdim==self.dim: #if there is no projection done then return the values calculated in the __init__
            return mu_proj,self.evec,self.sigma
        else: #project ellipsoid down to dim given by non-zero elements in proj_vars
            U,S,V=sl.svd(dot(T.T,self.Binv.T),full_matrices=False) #need to get the eigenvectors/values of the projection
            return mu_proj,U,S #return the mean, eigenvectors, and spread (NOTE: S is not scaled yet, that happens later)
        
    def inside(self,p):
        """This function will return a mask array that is True for data points in the ellipsoid and False otherwise
        Parameters
        ----------
        p : the percent of points contained in the ellipsoid, either a single value of a list of values
            i.e. 0.68 or [0.68,0.955].

        Return
        ------
        inside : A bool mask that is True for data 'p' percent points inside the ellipsoid.
                 If p is a list then the output will be an array of masks with
                 shape=(number of points,number of p values).
        """
        if self.data_t is None:
            return None
        try: #if a list get the length
            l=len(p)
        except: #if not then make it a list of length 1
            l=1
            p=[p]
        invp=st.chi.ppf(p,self.dim) #use the chi distribution to get the scaling factor for S (see, now we scale it)
        inside=zeros([self.n,l],dtype=bool) #list to hold the masks
        for jdx,j in enumerate(invp): #loop of each p values
            B=self.B/j #scale the B matrix
            for idx,i in enumerate(self.data_t): #loop over each data point
                A=dot(i,B) #dot it with B
                inside[idx,jdx]=(dot(A,A.T)<=1) #see if it is in the ellipsoid
        if l==1: #return 1d array if only 1 p value was given
            return inside.flatten()
        else: #return 2d array if multiple p values were given
            return inside
    def one_d_ellipse(self,proj_vars,p):
        """Return the 1d projection center and limts for given p values
        Parameters
        ----------
        proj_vars : array that is 1 for the projection dimension, and 0 other wise
                    i.e. array([0,0,1,0,0]) will project 5d ellipsoid onto the 3rd dim
        p         : the percent of points contained in the ellipsoid, either a single
                    value of a list of values i.e. 0.68 or [0.68,0.955], 

        Return
        ------
        mu  : the center of the ellipsoid projected down to 1D
        sig : the edges of and ellipsoid containing 'p' percent of points
              projected down to 1D. If p is a list then an array of limits will
              be returned, one for each p value.
        """
        mu,s=self.proj(proj_vars) #use the projection function
        return mu,s*st.chi.ppf(p,self.dim) #scale it using a chi distribution (see, now we scale it)
    def two_d_ellipse(self,proj_vars,p,**kwargs):
        """Return the 2d projection as a matplotlib Ellipse object for the given p values
        Parameters
        ----------
        proj_vars : array that is 1 for the projection dimension, and 0 other wise
                    i.e. array([0,0,1,0,1]) will project 5d ellipsoid onto the plane
                    span by the 3rd and 5th variable.
        p         : the percent of points contained in the ellipsoid, either a single
                     value of a list of values i.e. 0.68 or [0.68,0.955],

                     if p is a list then a list of Ellipse objects will be returned, one for each p value
        Keywords
        --------
        kwargs : keywords to be passed into the matplotlib Ellipse object

        Return
        ------
        ells : matplotlib Ellipse object
        """
        mu,u,s=self.proj(proj_vars) #get the mean, eigenvectors, and eigenvales for projected array
        try: #if a list get the length
            l=len(p)
        except: #if not then make it a list of length 1
            l=1
            p=[p]
        invp=st.chi.ppf(p,self.dim) #scale it using a chi distribution (see, now we scale it)
        angle=rad2deg(arctan2(u[0,1],u[0,0])) #angle the first eignevector makes with the x-axis
        ells=[] #list to hold the Ellipse objects
        for i in invp:
            ells.append(Ellipse(xy=mu,width=s[0]*i*2,height=s[1]*i*2,angle=angle,**kwargs))#make the Ellipse objects, the *2 is needed since Ellipse takes the full axis vector
        if l==1: #if only one p values was given return the Ellipse object (not as a list)
            return ells[0]
        else: #else return the list of Ellipse objects
            return ells

if __name__=='__main__':
    #lets test the code and make sure it works as expected
    from pylab import * #import this so we can plot
    N=10000 #number of data points to simulate
    Data=np.random.multivariate_normal([2.,-2.,0.],[[1,.25,-.5],[.25,.1,.6],[-.5,.6,5]],size=N)#random draw some data
    E=Ellipse_proj(Data) #make my Ellipse_proj object
    inside=E.inside([0.68,.955]) #get the list of points inside the 1 and 2 sigma confidence intervals
    print 'percent inside: ',inside[:,0].sum()/float(N),inside[:,1].sum()/float(N) #check that the correct fraction is being found

    #plot all 1d and 2d projections of the data with the 1 and 2 sigma confidence intervals
    
    ax1=subplot(3,3,1)
    a,b,c=ax1.hist(Data.T[0],25,normed=True,histtype='step') #make a histogram
    mu,sig=E.one_d_ellipse(array([1,0,0]),[0.68,.955]) #get the mean and sigma points
    ax1.vlines([mu,mu+sig[0],mu-sig[0],mu+sig[1],mu-sig[1]],0,a.max()) #plot the mean and confidence intervals as virtical lines

    ax2=subplot(3,3,4)
    ax2.plot(Data.T[0],Data.T[1],',') #plot 2d projection
    el=E.two_d_ellipse(array([1,1,0]),[0.68,.955],alpha=0.6,zorder=30) #get the matplotlib Ellipse objects for the 1 and 2 sigma levels
    el[1].set_facecolor('#ff0000') #change one ellipse color to red
    ax2.add_artist(el[1]) #add the object to the plot
    ax2.add_artist(el[0])

    ax3=subplot(3,3,5)
    a,b,c=ax3.hist(Data.T[1],25,normed=True,histtype='step')
    mu,sig=E.one_d_ellipse(array([0,1,0]),[0.68,.955])
    ax3.vlines([mu,mu+sig[0],mu-sig[0],mu+sig[1],mu-sig[1]],0,a.max())

    ax4=subplot(3,3,7)
    ax4.plot(Data.T[0],Data.T[2],',') #plot the data
    el=E.two_d_ellipse(array([1,0,1]),[0.68,0.955],alpha=0.6,zorder=30)
    el[1].set_facecolor('#ff0000')
    ax4.add_artist(el[1])
    ax4.add_artist(el[0])

    ax5=subplot(3,3,8)
    ax5.plot(Data.T[1],Data.T[2],',') #plot the data
    el=E.two_d_ellipse(array([0,1,1]),[0.68,0.955],alpha=0.6,zorder=30)
    el[1].set_facecolor('#ff0000')
    ax5.add_artist(el[1])
    ax5.add_artist(el[0])
    
    ax6=subplot(3,3,9)
    a,b,c=ax6.hist(Data.T[2],25,normed=True,histtype='step')
    mu,sig=E.one_d_ellipse(array([0,0,1]),[0.68,0.955])
    ax6.vlines([mu,mu+sig[0],mu-sig[0],mu+sig[1],mu-sig[1]],0,a.max())

    """
    ells=error_ellipse(Data,[0.05,0.68,0.955],alpha=0.2)
    ax=gca()
    colors=[[1,0,0],[0,1,0],[0,0,1]]
    for e,c in zip(ells,colors):
        ax.add_artist(e)
        e.set_facecolor(c)
    """
    show()
