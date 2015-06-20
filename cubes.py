import astropy
import astropy.units as u
import astropy.constants as c
import numpy as np
import pyfits
import matplotlib.pyplot as plt

class PhysicalQuantity:
    """Physical quantity : has a value, a unit, a description, and a type"""
    def __init__(self,
                 value=0.0, # The numerical value of the quantity
                 unit=u.dimensionless_unscaled, # The unit of the quantity
                 descr="", # A character string to describe the quantity, e.g. "Total proton density", "Dust temperature", "magnetic field"...
                 typ="scalar" # The tensorial nature of the quantity. Can be "scalar" or "vector", in which case the value must be an array of 2 values for a 2D vector, of 3 values for a 3D vector
                 ):
        self.v=value*unit
        self.descr=descr

class Axis:
    """Coordinate axis : a 1D array of physical quantities. WARNING : BETTER USE ONLY REGULAR AXES! AND USE FLOATS !"""
    def __init__(self,
                 value=np.empty(shape=(0),dtype=np.float64), # The values of the coordinate axis
                 unit=u.dimensionless_unscaled,
                 descr=""):
        self.v=value*unit
        self.unit=unit
        self.descr=descr

    def setcenters(self,centers):
        u=self.v.unit
        self.v=centers*u

    def Npts(self):
        return(self.v.shape[0])

    def getintermediatewalls(self):
        Npts=self.Npts()
        result=np.zeros(shape=(Npts-1),dtype=self.v.dtype)
        for i in range(Npts-1):
            val=0.5*(self.v[i]+self.v[i+1])
            result[i]=val.value
        return(result*self.v[0].unit)

    def getcellsizes(self):
        Npts=self.Npts()
        result=np.zeros(shape=(Npts),dtype=self.v.dtype)
        iw=self.getintermediatewalls()        
        for i in range(Npts-1):
            val=2.0*(iw[i].value-self.v[i].value)
            result[i]=val
        result[Npts-1]=2.0*(self.v[Npts-1].value-iw[Npts-2].value)
        return(result*self.v[0].unit)
            
    def getwalls(self):
        Npts=self.Npts()
        result=np.zeros(shape=(Npts+1),dtype=self.v.dtype)
        iw=self.getintermediatewalls()
        cs=self.getcellsizes()
        result[1:Npts]=iw[0:Npts-1].value
        result[0]=iw[0].value-cs[0].value
        result[Npts]=iw[Npts-2].value+cs[Npts-1].value
        return(result*self.v[0].unit)

    def setwalls(self,walls):
        Npts=walls.shape[0]-1
        centers=np.zeros(shape=(Npts),dtype=walls.dtype)
        for i in range(Npts):
            centers[i]=0.5*(walls[i]+walls[i+1])        
        self.setcenters(centers)

class Beam():
    """Beam : contains two physical quantities for major and minors axes of the (assumed Gaussian) beam, plus a position angle"""
    def __init__(self,bmaj,bmin,bpa):
        self.bmaj=bmaj
        self.bmin=bmin
        self.bpa=bpa

class Map():
    """Map : a 2D array of physical quantities, two axes, a type (scalar or vector), and a beam"""
    def __init__(self,
                 value=np.empty(shape=(0,0),dtype=np.float64),
                 unit=u.dimensionless_unscaled,
                 descr="",
                 axis1=Axis(descr="X",unit=u.dimensionless_unscaled),
                 axis2=Axis(descr="Y",unit=u.dimensionless_unscaled),
                 typ="scalar"):
        self.v=value*unit
        self.unit=unit
        self.descr=descr
        self.axis1=axis1
        self.axis2=axis2
        self.type=typ
        self.beam=None
    
    def readfromfits(self,fitsfile):
        # Read Header        
        h=(pyfits.open(fitsfile)[0]).header
        # Read data
        d=((pyfits.open(fitsfile)[0]).data).astype('>f4')         
        # Read data unit : COMMENTED OUT FOR THE MOMENT AS TEST MAPS HAVE NO BUNIT KEYWORD
        #self.unit=u.Unit(h['BUNIT'])        
        self.v=d*self.unit
        # Read axes units, automatic conversion from the strings read in CUNITi keywords. THIS MAY LEAD TO UNEXPECTED PROBLEMS...
        cunits=[None,None]
        for i in [1,2]:cunits[i-1]=u.Unit(h['CUNIT'+str(i)].lower())
        self.axis1=Axis(value=float(h['CRVAL1'])+float(h['CDELT1'])*(np.linspace(0.0,float(int(h['NAXIS1']))-1.0,num=int(h['NAXIS1']))-float(h['CRPIX1'])),unit=cunits[0],descr=h['CTYPE1'])
        self.axis2=Axis(value=float(h['CRVAL2'])+float(h['CDELT2'])*(np.linspace(0.0,float(int(h['NAXIS2']))-1.0,num=int(h['NAXIS2']))-float(h['CRPIX2'])),unit=cunits[1],descr=h['CTYPE2'])
        # Retrieve the beam. Note that this assumes the same unit for CUNIT1 and CUNIT2, and that BPA unit is degrees
        bmaj=float(h['BMAJ'])*cunits[0]
        bmin=float(h['BMIN'])*cunits[0]
        bpa=float(h['BPA'])*u.deg
        self.beam=Beam(bmaj,bmin,bpa)

    def writetofits(self,fitsfile):
        pass

    def show_stats(self):
        print("Statistics for {0}:\n Unit : {1}\n Min={2} ; Max={3} ; Mean={4} ; Median={5} ; StdDev={6} ".format(self.descr,self.unit,np.min(self.v.value),np.max(self.v.value),np.mean(self.v.value),np.median(self.v.value),np.std(self.v.value)))#

    def PDF(self,scaling="lin",nbins=200,normed=0,density=True):
        v_1D=np.reshape(self.v,np.prod(np.shape(self.v))).value
        if(scaling=="lin"):hist,bin_edges=np.histogram(v_1D, nbins, normed=normed,density=density)
        if(scaling=="log"):hist,bin_edges=np.histogram(np.log10(v_1D), nbins, normed=normed,density=density)        
        bincenters=np.zeros_like(hist)
        for i in range(len(hist)):bincenters[i]=0.5*(bin_edges[i]+bin_edges[i+1])
        return(bincenters,hist)

class Cube():
    """Cube : a 3D array of physical quantities, three axes, and a type (scalar or vector)"""
    def __init__(self,
                 value=np.empty(shape=(0,0,0),dtype=np.float64),
                 unit=u.dimensionless_unscaled,
                 descr="",
                 axis1=Axis(descr="X",unit=u.dimensionless_unscaled),
                 axis2=Axis(descr="Y",unit=u.dimensionless_unscaled),
                 axis3=Axis(descr="Z",unit=u.dimensionless_unscaled),
                 typ="scalar"):
        self.v=value*unit
        self.unit=unit
        self.descr=descr
        self.axis1=axis1
        self.axis2=axis2
        self.axis3=axis3
        self.type=typ
    
    def readfromfits(self,fitsfile):
        # Read Header        
        h=(pyfits.open(fitsfile)[0]).header
        #print h
        # Read data
        d=((pyfits.open(fitsfile)[0]).data).astype('>f4') 
        # Read data unit, automatic conversion from the string read in BUNIT keyword
        self.unit=u.Unit(h['BUNIT'])
        self.v=d*self.unit        
        # Read axes units, automatic conversion from the strings read in CUNITi keywords. THIS MAY LEAD TO UNEXPECTED PROBLEMS...
        cunits=[None,None,None]
        for i in [1,2,3]:cunits[i-1]=u.Unit(h['CUNIT'+str(i)])
        self.axis1=Axis(value=float(h['CRVAL1'])+float(h['CDELT1'])*(np.linspace(0.0,float(int(h['NAXIS1']))-1.0,num=int(h['NAXIS1']))-float(h['CRPIX1'])),unit=cunits[0],descr=h['CTYPE1'])
        self.axis2=Axis(value=float(h['CRVAL2'])+float(h['CDELT2'])*(np.linspace(0.0,float(int(h['NAXIS2']))-1.0,num=int(h['NAXIS2']))-float(h['CRPIX2'])),unit=cunits[1],descr=h['CTYPE2'])
        self.axis3=Axis(value=float(h['CRVAL3'])+float(h['CDELT3'])*(np.linspace(0.0,float(int(h['NAXIS3']))-1.0,num=int(h['NAXIS3']))-float(h['CRPIX3'])),unit=cunits[2],descr=h['CTYPE3'])

    def readfromascii(self,asciifile):
        pass

    def writetofits(self,fitsfile):
        pass

    def rotate(self,rotation_axis,angle):
        pass

    def project(self,projection_axis):
        pass

    def show_stats(self):
        print("Statistics for {0}:\n Unit : {1}\n Min={2} ; Max={3} ; Mean={4} ; Median={5} ; StdDev={6} ".format(self.descr,self.unit,np.min(self.v.value),np.max(self.v.value),np.mean(self.v.value),np.median(self.v.value),np.std(self.v.value)))#

    def getdims(self):
        return(self.v.shape)

    def show_axes(self):
        for a in [self.axis1,self.axis2,self.axis3]:
            print("axis {0} : {1} points, from {2} to {3}".format(a.descr,a.v.shape[0],a.v[0],a.v[a.v.shape[0]-1]))

    def PDF(self,scaling="lin",nbins=200,normed=0,density=True):
        v_1D=np.reshape(self.v,np.prod(np.shape(self.v))).value
        if(scaling=="lin"):hist,bin_edges=np.histogram(v_1D, nbins, normed=normed,density=density)
        if(scaling=="log"):hist,bin_edges=np.histogram(np.log10(v_1D), nbins, normed=normed,density=density)        
        bincenters=np.zeros_like(hist)
        for i in range(len(hist)):bincenters[i]=0.5*(bin_edges[i]+bin_edges[i+1])
        return(bincenters,hist)

if __name__ == "__main__":
    print 'toto'
    nH=Cube(descr="Total proton density")
    nH.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-dens-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    nH.show_stats()
    nH.show_axes()
#    nH.PDF(xscal="log",yscal="log",ymin=6e-5,ymax=10.,showmedian=True,showmean=True,show=False)
    Bx=Cube(descr="X component of the B field")
    Bx.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnX-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    Bx.show_stats()
    By=Cube(descr="Y component of the B field")
    By.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnY-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    Bz=Cube(descr="Z component of the B field")
    Bz.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnZ-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")



#    Bx.show_stats()
#    Bx.show_axes()
#    Bx.PDF(xscal="lin",yscal="log",showmedian=True,showmean=True)
#    M=Map(descr="I",value=np.array([[0,1],[2,3]]),unit=u.K)
#    print M.descr,M.axis1.descr,M.v
#    nH=Cube(descr="Total proton density",value=np.zeros((5,3,2),dtype=np.float64),unit=u.cm**-3,typ="scalar")
#    print nH.descr,nH.v.unit
#    p=PhysicalQuantity(value=1.4,unit=u.cm,descr="position")
#    print p.v.value,p.v.unit,p.descr
#    x=Axis(descr="X",unit=u.cm)
#    print x.v.value,x.v.unit,x.descr
#    x.setcenters(np.linspace(0.,100.,num=11))
#    print "centers",x.v
#    print "iw",x.getintermediatewalls()
#    print "cs",x.getcellsizes()
#    print "cw",x.getwalls()
#    print x.v.value
#
#    x.setwalls(np.array([-1.,0.,1.,2.,3.,4.]))
#    print "centers",x.v
#    print "iw",x.getintermediatewalls()
#    print "cs",x.getcellsizes()
#    print "cw",x.getwalls()
#    print x.v.value

#    x.setwalls(np.linspace(-2,3,num=10)*u.cm)
#    print x.v.value,x.v.unit,x.descr
#    I=Map(descr="I",value=np.array([[0,1],[2,3]]),unit=u.km/u.s)
#    print I.v.value,I.v.unit,I.descr
#    c=Cube(descr="c",unit=u.pc)
#    print c.v.value,c.v.unit,c.descr
#    print x.values
#    print x.unit
#    x.regular(10,-3.,3.)
#    print x.values
