import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
from matplotlib.patches import Ellipse
from cubes import *

# Plot a Map instance M
def fig_map(M,data_range=None,cmap=plt.get_cmap('hsv'),scale="lin",beaminfo=None,title="",figfile="",axisunit=None,xvar=None,yvar=None):
    # Create the figure canvas
    fig=plt.figure(figsize=(8.8,8.3))
    # if no specific data range is input, use the whole range
    if (data_range==None):data_range=[np.min(M.v.value),np.max(M.v.value)]
    # Build color scale
    D_norm = mpl.colors.Normalize(vmin=np.min(data_range), vmax=np.max(data_range))
    # Insert axes on figure canvas
    ax1 = fig.add_axes([0.03,0.15,0.85,0.8])
    # If no specific unit is wanted, use the one from axis1 - THIS AS ALWAYS ASSUMES THE SAME UNIT IS USED IN BOTH AXES    
    if (axisunit==None):axisunit=M.axis1.unit 
    # Build X and Y extent of the map
    extent=[M.axis1.getwalls()[0].to(axisunit).value,
            M.axis1.getwalls()[-1].to(axisunit).value,
            M.axis2.getwalls()[0].to(axisunit).value,
            M.axis2.getwalls()[-1].to(axisunit).value]
    # Plot the map (linear or logarithmic scalings)
    if (scale=="lin"):ax1.imshow(M.v.value, interpolation='nearest',origin='lower',cmap=cmap,extent=extent,norm=D_norm)
    if (scale=="log"):ax1.imshow(np.log(M.v.value), interpolation='nearest',origin='lower',cmap=cmap,extent=extent)
    if (scale=="log10"):ax1.imshow(np.log10(M.v.value), interpolation='nearest',origin='lower',cmap=cmap,extent=extent) 
    # If no specific names for the axes are given, use descriptors of axes in the Map instance
    if (xvar==None):xvar=M.axis1.descr
    if (yvar==None):xvar=M.axis2.descr
    # Add axes labels
    ax1.set_xlabel("$"+xvar+"$ "+"[{0:latex}]".format(axisunit),fontsize=20)
    ax1.set_ylabel("$"+yvar+"$ "+"[{0:latex}]".format(axisunit),fontsize=20)
    # Add title
    ax1.set_title(title+"\n",fontsize=20)
    # Add beam as black-rimmed, yellow-filled ellipse in lower left corner
    ax1.add_artist(Ellipse(xy=[0.9*extent[0]+0.1*extent[1],0.9*extent[2]+0.1*extent[3]], 
                           width=M.beam.bmin.to(axisunit).value, 
                           height=M.beam.bmaj.to(axisunit).value, 
                           angle=M.beam.bpa.to(u.deg).value,
                           fc='yellow',ec='black',lw=1))      
    # Specify the X and Y extent of the map
    ax1.axis(extent)
    # Change tick lengths and tick number sizes to more readable values
    for tick in ax1.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(1)
        tick.tick2line.set_markeredgewidth(1)
        tick.tick1line.set_markersize(18)
        tick.tick2line.set_markersize(18)
    for tick in ax1.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(1)
        tick.tick2line.set_markeredgewidth(1)
        tick.tick1line.set_markersize(18)
        tick.tick2line.set_markersize(18)
    # Retrieve locations and labels of ticks on the y axis
    locs, labels = plt.yticks()
    # Rotate these y axis labels by 90 degrees in accordance with Planck Style Guide
    plt.setp(labels, rotation=90)
    # Insert new set of axes to plot the color bar to the right of the map
    ax2 = fig.add_axes([0.85,0.15,0.03, 0.8])
    # Plot the color bar
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=D_norm,orientation='vertical')
    # Add a label to the color bar specifying the units of the map
    cb1.set_label(M.unit,fontsize=20)
    # Select whether we save the figure to a file (with specified name) or we just plot it on the screen
    if (figfile!=""):
        plt.savefig(figfile)
        plt.clf()
    else:
        plt.show()


# Plot a Map instance M, overlaid with contours using a (possibly different) Map instance C. Use lv to specify levels, cw to specify contour width and cc to specify contour colors
def fig_map_contours(M,C,lv=None,cw=None,cc=None,data_range=None,cmap=plt.get_cmap('hsv'),scale="lin",beaminfo=None,title="",figfile="",axisunit=None,xvar=None,yvar=None):
    # Create the figure canvas
    fig=plt.figure(figsize=(8.8,8.3))
    # if no specific data range is input, use the whole range
    if (data_range==None):data_range=[np.min(M.v.value),np.max(M.v.value)]
    # Build color scale
    D_norm = mpl.colors.Normalize(vmin=np.min(data_range), vmax=np.max(data_range))
    # Insert axes on figure canvas
    ax1 = fig.add_axes([0.03,0.15,0.85,0.8])
    # If no specific unit is wanted, use the one from axis1 - THIS AS ALWAYS ASSUMES THE SAME UNIT IS USED IN BOTH AXES    
    if (axisunit==None):axisunit=M.axis1.unit 
    # Build X and Y extent of the map
    extent=[M.axis1.getwalls()[0].to(axisunit).value,
            M.axis1.getwalls()[-1].to(axisunit).value,
            M.axis2.getwalls()[0].to(axisunit).value,
            M.axis2.getwalls()[-1].to(axisunit).value]
    # Plot the map (linear or logarithmic scalings)
    if (scale=="lin"):ax1.imshow(M.v.value, interpolation='nearest',origin='lower',cmap=cmap,extent=extent,norm=D_norm)
    if (scale=="log"):ax1.imshow(np.log(M.v.value), interpolation='nearest',origin='lower',cmap=cmap,extent=extent)
    if (scale=="log10"):ax1.imshow(np.log10(M.v.value), interpolation='nearest',origin='lower',cmap=cmap,extent=extent) 
    # Add contours : if levels are not specified, they are computed automatically. 
    if (lv==None):CS=ax1.contour(C.v.value,extent=extent,linewidths=cw,colors=cc)
    else:CS=ax1.contour(C.v.value,lv,extent=extent,linewidths=cw,colors=cc)
    # Add labels to the contours
    plt.clabel(CS,inline=1,fontsize=10)
    # If no specific names for the axes are given, use descriptors of axes in the Map instance
    if (xvar==None):xvar=M.axis1.descr
    if (yvar==None):xvar=M.axis2.descr
    # Add axes labels
    ax1.set_xlabel("$"+xvar+"$ "+"[{0:latex}]".format(axisunit),fontsize=20)
    ax1.set_ylabel("$"+yvar+"$ "+"[{0:latex}]".format(axisunit),fontsize=20)
    # Add title
    ax1.set_title(title+"\n",fontsize=20)
    # Add beam as black-rimmed, yellow-filled ellipse in lower left corner
    ax1.add_artist(Ellipse(xy=[0.9*extent[0]+0.1*extent[1],0.9*extent[2]+0.1*extent[3]], 
                           width=M.beam.bmin.to(axisunit).value, 
                           height=M.beam.bmaj.to(axisunit).value, 
                           angle=M.beam.bpa.to(u.deg).value,
                           fc='yellow',ec='black',lw=1))      
    # Specify the X and Y extent of the map
    ax1.axis(extent)
    # Change tick lengths and tick number sizes to more readable values
    for tick in ax1.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(1)
        tick.tick2line.set_markeredgewidth(1)
        tick.tick1line.set_markersize(18)
        tick.tick2line.set_markersize(18)
    for tick in ax1.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(1)
        tick.tick2line.set_markeredgewidth(1)
        tick.tick1line.set_markersize(18)
        tick.tick2line.set_markersize(18)
    # Retrieve locations and labels of ticks on the y axis
    locs, labels = plt.yticks()
    # Rotate these y axis labels by 90 degrees in accordance with Planck Style Guide
    plt.setp(labels, rotation=90)
    # Insert new set of axes to plot the color bar to the right of the map
    ax2 = fig.add_axes([0.85,0.15,0.03, 0.8])
    # Plot the color bar
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=D_norm,orientation='vertical')
    # Add a label to the color bar specifying the units of the map
    cb1.set_label(M.unit,fontsize=20)
    # Select whether we save the figure to a file (with specified name) or we just plot it on the screen
    if (figfile!=""):
        plt.savefig(figfile)
        plt.clf()
    else:
        plt.show()


## Plot a PDF from data held in a Cube or Map instance
#def plot_PDF(data,nbins=200,normed=0,density=True,figfile="",xscal="lin",yscal="lin",showmean=False,showmedian=False,ymin=None,ymax=None,xmin=None,xmax=None):
#    # Reshape the data values to 1D array
#    data_1D=np.reshape(data.v.value,np.prod(np.shape(data.v.value)))
#    # Build the histogram, i.e. the centers of bins and the hstogram values for each bin
#    bincenters,hist=data.PDF(scaling=xscal,nbins=nbins,normed=normed,density=density)   
#    # Create the figure canvas
#    fig = plt.figure(figsize=(17,10))
#    # Insert axes on figure canvas
#    ax = fig.add_subplot(111,position=[0.12,0.12,0.8,0.8])
#    # Plot the histogram, depending on the scaling required for x and y axes
#    if (xscal=="lin")and(yscal=="lin"):plt.plot(bincenters,hist,'k-',lw=4)
#    if (xscal=="lin")and(yscal=="log"):plt.semilogy(bincenters,hist,'k-',lw=4)
#    if (xscal=="log")and(yscal=="lin"):plt.semilogx(np.power(10.,bincenters),hist,'k-',lw=4)
#    if (xscal=="log")and(yscal=="log"):plt.loglog(np.power(10.,bincenters),hist,'k-',lw=4)
#    # If specified, show the mean of the data in red
#    if(showmean):plt.axvline(np.mean(data_1D),color='r',lw=4)
#    # If specified, show the median of the data in blue
#    if(showmedian):plt.axvline(np.median(data_1D),color='b',lw=4)
#    # If specific values are given for the minimum and maximum values to plot in y, use them
#    if (ymin!=None)and(ymax!=None):ax.set_ylim(ymin,ymax)
#    # Insert labels for x and y axes. The one on the x axis is based on data descriptor and unit, the one on the y axis is simply "Distribution function"
#    ax.set_xlabel(data.descr+" "+"[{0:latex_inline}]".format(data.unit),fontsize=35)
#    ax.set_ylabel("Distribution function"+"\n",fontsize=35)
#    # Change size of ticks on x and y axes for better readability
#    for tick in ax.xaxis.get_major_ticks():tick.label.set_fontsize(25) 
#    for tick in ax.yaxis.get_major_ticks():tick.label.set_fontsize(25) 
#    # Retrieve locations and labels of ticks on the y axis
#    locs, labels = plt.yticks()
#    # Rotate these y axis labels by 90 degrees in accordance with Planck Style Guide
#    plt.setp(labels, rotation=90)
#    # Select whether we save the figure to a file (with specified name) or we just plot it on the screen
#    if (figfile!=""):
#        plt.savefig(figfile)
#        plt.clf()
#    else:
#        plt.show()


# Plot a PDF from data held in a Cube or Map instance
def plot_PDF(input_data,nbins=200,normed=0,density=True,figfile="",xscal="lin",yscal="lin",showmean=False,showmedian=False,ymin=None,ymax=None,xmin=None,xmax=None,descr=None,colors=None):
    # Create the figure canvas
    fig = plt.figure(figsize=(17,10))
    # Insert axes on figure canvas
    ax = fig.add_subplot(111,position=[0.12,0.12,0.8,0.8])
    # Identify whether data is a list or a single instance of either Map or Cube
    if not(type(input_data) is list): input_data=[input_data]
    for data in input_data:
        # Reshape the data values to 1D array
        data_1D=np.reshape(data.v.value,np.prod(np.shape(data.v.value)))
        # Build the histogram, i.e. the centers of bins and the hstogram values for each bin
        bincenters,hist=data.PDF(scaling=xscal,nbins=nbins,normed=normed,density=density)   
        # Plot the histogram, depending on the scaling required for x and y axes
        if (xscal=="lin")and(yscal=="lin"):plt.plot(bincenters,hist,ls='-',lw=4,label=data.descr)
        if (xscal=="lin")and(yscal=="log"):plt.semilogy(bincenters,hist,ls='-',lw=4,label=data.descr)
        if (xscal=="log")and(yscal=="lin"):plt.semilogx(np.power(10.,bincenters),hist,ls='-',lw=4,label=data.descr)
        if (xscal=="log")and(yscal=="log"):plt.loglog(np.power(10.,bincenters),hist,ls='-',lw=4,label=data.descr)
        # If specified, show the mean of the data in red
        if(showmean):plt.axvline(np.mean(data_1D),color='r',lw=4)
        # If specified, show the median of the data in blue
        if(showmedian):plt.axvline(np.median(data_1D),color='b',lw=4)        
    # If specific values are given for the minimum and maximum values to plot in y, use them
    if (ymin!=None)and(ymax!=None):ax.set_ylim(ymin,ymax)
    # If no specific descriptor is given, use the one form the data
    if (descr==None):descr=data.descr
    # Insert labels for x and y axes. The one on the y axis is simply "Distribution function"
    ax.set_xlabel(descr+" "+"[{0:latex_inline}]".format(data.unit),fontsize=35)
    ax.set_ylabel("Distribution function"+"\n",fontsize=35)
    # Change size of ticks on x and y axes for better readability
    for tick in ax.xaxis.get_major_ticks():tick.label.set_fontsize(25) 
    for tick in ax.yaxis.get_major_ticks():tick.label.set_fontsize(25) 
    # Retrieve locations and labels of ticks on the y axis
    locs, labels = plt.yticks()
    # Rotate these y axis labels by 90 degrees in accordance with Planck Style Guide
    plt.setp(labels, rotation=90)
    plt.legend()
    # Select whether we save the figure to a file (with specified name) or we just plot it on the screen
    if (figfile!=""):
        plt.savefig(figfile)
        plt.clf()
    else:
        plt.show()


if __name__ == "__main__":
    print 'tut'
    I=Map(descr="Total intensity")
    I.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Output-Planck-beamed/I_y0_z_T_sm15.0.fits")
    I.show_stats()
    p=Map(descr="Polarization fraction")
    p.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Output-Planck-beamed/PI_y0_z_T_sm15.0.fits")
    p.show_stats()
#    fig_map(I,axisunit=u.deg,xvar="l",yvar="b")
#    fig_map_contours(I,p,lv=[0.05,0.1,0.15],cw=2,cc='k',axisunit=u.deg,xvar="l",yvar="b")
    nH=Cube(descr="Total proton density")
    nH.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-dens-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    nH.show_stats()
    nH.show_axes()
#    plot_PDF(nH,figfile="",xscal="log",yscal="log",ymin=6e-5,ymax=10.)
#    plot_PDF(I,nbins=200,normed=0,density=True,figfile="",xscal="log",yscal="log",showmean=True,showmedian=True,ymin=6e-5,ymax=10.,xmin=None,xmax=None)

    Bx=Cube(descr="X component of the B field")
    Bx.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnX-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    By=Cube(descr="Y component of the B field")
    By.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnY-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")
    Bz=Cube(descr="Z component of the B field")
    Bz.readfromfits("/Users/levrier/PLANCK/Simulations/Data/Clump-0/Original-data/SF-THY3D_grav_mag_bcl-036-magnZ-X0.52-Y0.78-Z0.82-18.0pc-Q9.fits")

    plot_PDF([Bx,By,Bz],yscal="log",descr="Component of the B field",ymin=2.0,ymax=8e5)
