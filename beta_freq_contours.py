#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import sys, os.path
import scipy as np
from scipy.constants import c, hbar, pi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=14)
matplotlib.rc('text.latex', preamble = \
        '\usepackage{amsmath}, \usepackage{palatino}, \usepackage{upgreek}')
#matplotlib.rc('font',**{'family':'serif','serif':['palatino, times']})  ## select fonts
matplotlib.rc('font',**{'family':'sans-serif', 'sans-serif':['Computer Modern Sans serif']})


## ==== User settings ====

## Select which structure properties will be plot and which not (each in separate file)
quantities = []
quantities = ['reflection', 'transmission', 'loss']
#quantities += ['absNimag']
#quantities += ['absNre']
#quantities += ['Nre']
#quantities += ['eps']
#quantities += ['mu']

## Settings for the parameter that was scanned through
#parameter_name = 'Ky'
#ylabel = u"$s$-polarisation incidence angle $\\beta_0$ [rad]"
#parameter_name = 'Kx'
#ylabel = u"$p$-polarisation incidence angle $\\beta_0$ [rad]"
#logarithmic = 0
#ylim = (0, np.pi/2); recalculate_to_angle = True
#minf, maxf  = 0, 1e12           # span of the horizontal axis to be plot
#frequnit = 1e12
#yunit = 1
#interp_anisotropy = .05       # value lower than 1. interpolates rather vertically
#xlim = 0,  1.

parameter_name = 'xradius'
ylabel = 'X-radius [$\\upmu$m]'
#parameter_name = 'yradius'
#ylabel = 'Y-radius [$\\upmu$m]'
#parameter_name = 'zradius'
#ylabel = 'Z-radius [$\\upmu$m]'
recalculate_to_angle = False
#
# Settings for the plot range
logarithmic = 1
frequnit = 1e12
yunit = 1e-6
minf, maxf  = .5e12, 3e12           # span of the horizontal axis to be plot
xlim        = (minf/frequnit, maxf/frequnit)
ylim = (5.01, 30.01)
interp_anisotropy = .01       # value lower than 1. interpolates rather vertically; optimize if plot desintegrates
## ==== / User settings ====


## Plot each of quantities in separate file
for quantity in quantities:

    ##Start figure + subplot
    plt.figure(figsize=(7,5))

    ## Load data from multiple files
    print "Got %d files to plot %s..." % (len(sys.argv), quantity)
    x, y, z, z2 = [np.array([]) for _ in range(4)] ## three empty arrays

    filenames = sys.argv[1:]
    #filenames.sort(key=lambda name: float(name.split('radius=')[1].split('_')[0]))     # sort (optional)

    for datafile_name in filenames: 
        ## Getting 1D data
        (freq, s11_ampli, s11p, s12_ampli, s12p, Nre, Nim, Zre, Zim, eps_r, eps_i, mu_r, mu_i) = \
                np.loadtxt(datafile_name, usecols=range(13), unpack=True)
        if quantity == 'reflection':        znew = s11_ampli
        elif quantity == 'transmission':    znew = s12_ampli
        elif quantity == 'loss':            znew = np.log(1 - s11_ampli**2 - s12_ampli**2)
        elif quantity == 'absNimag':        znew = np.clip(abs(Nim), 0, 10)
            #znew = np.log10(np.clip(abs(Nim), 0, 300 ))
        elif quantity == 'absNre':          
            znew = abs(np.arcsin(np.sin(np.real(Nre*freq*100e-6/c) * np.pi)) / np.pi)
            znew2 = np.clip(abs(Nim), 0, 10)
        elif quantity == 'Nre':             znew = Nre #np.real(Nre*freq*100e-6/c)
        elif quantity == 'eps':             znew = eps_r
        elif quantity == 'epsangle':        znew = np.angle(eps_r + 1j*eps_i)
        elif quantity == 'mu':              znew = mu_r
        else: print 'Error!, select a known quantity to plot!'

        ## Truncate the data ranges
        truncated = np.logical_and(freq>minf, freq<maxf)
        (freq, znew) = map(lambda x: x[truncated], [freq, znew])
        #(znew2) = map(lambda x: x[truncated], [znew2]) ## XXX

        ## Load header
        parameters  = {}
        columns     = []
        with open(datafile_name) as datafile:
            for line in datafile:
                if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
                key, value = line.replace(',', ' ').split()[-2:]
                try: value = float(value) ## Try to convert to float, if possible
                except: pass                ## otherwise keep as string
                parameters[key] = value
            for line in datafile:
                if ('column' in line.lower()): columns.append(line.strip().split(' ', 1)[-1])

        print parameters[parameter_name],

        if not recalculate_to_angle:
            ynew = parameters[parameter_name] * np.ones_like(freq) / yunit
        else:
            try:
                ynew = np.arcsin(parameters[parameter_name]/(2*pi*freq/c)) / yunit
            except:
                print freq, c
                print parameters[parameter_name], (2*pi*freq/c)
        x = np.append(x, freq/frequnit)
        y = np.append(y, ynew)
        z = np.append(z, znew)
        #z2 = np.append(z2, znew2) ## XXX

    if not ylim: ylim=(min(y), max(y))
    plt.ylabel(ylabel); 
    plt.xlabel(u"Frequency [THz]");

    xi = np.linspace(min(x), max(x), 500)
    yi = np.linspace(min(y), max(y), 300)
    # grid the data.
    from matplotlib.mlab import griddata
    zi = griddata(x, y*interp_anisotropy, z, xi, yi*interp_anisotropy, interp='linear')
    #zi2 = griddata(x, y*interp_anisotropy, z2, xi, yi*interp_anisotropy, interp='linear')
    # contour the gridded data, plotting dots at the nonuniform data points.

    if quantity in ('reflection', 'transmission'):        
                                cmap = cm.jet;      levels = np.arange(0.,1.05,.03)
    elif quantity=='loss':      cmap = cm.jet;      levels = np.arange(-10.,0.01,.2)
    elif quantity=='absNimag':  cmap = cm.jet;      levels = np.arange(0.,10.1,.2)
    elif quantity=='absNre':    cmap = cm.GnBu;     levels = np.arange(0.,.51,.05)
    elif quantity=='eps':       cmap = cm.RdBu_r;   levels = np.arange(-10.,10,.5)
    elif quantity=='epsangle':  cmap = cm.jet;      levels = np.arange(-np.pi,np.pi,.1)
    elif quantity=='mu':        cmap = cm.RdBu_r;   levels = np.arange(-10.,10,.5)
    else:                       cmap = cm.jet;      levels = np.arange(-2.,10,.1)

    # Standard contour plot
    contours = plt.contourf(xi,yi,zi, cmap=cmap, levels=levels)  
    for contour in contours.collections: contour.set_antialiased(False) ## fix aliasing for old Matplotlib

    # 3-shade contour plot for visualisation of photonic bands
    #contours = plt.contourf(xi, yi, zi, cmap=cmap, levels=[0., 0.03, 0.48, .5])  
    #for contour in contours.collections: contour.set_antialiased(False) ## fix aliasing for old Matplotlib

    ## Black for high permittivity
    #contours = plt.contourf(xi, yi, zi2, colors='k', alpha=1, levels=[2.5, 100])    ## XXX for eps100
    #contours = plt.contourf(xi, yi, zi2, colors='k', alpha=1, levels=[.8, 100])    ## XXX for eps012
    #for contour in contours.collections: contour.set_antialiased(False) ## fix aliasing for old Matplotlib


    #plt.clabel(plt.contour(xi, yi, zi, linewidths=0.5, levels=np.arange(0.,1.01,.25), colors='k'))
    #plt.clabel(plt.contour(xi, yi, zi, linewidths=0.5, levels=np.arange(.15), colors='k'))
    # plot data points.
    #plt.scatter(x, y, marker='o', color='black', s=.1, alpha=.2)
    #plt.xlim(0, max(freq)/frequnit)
    #plt.ylim(ylim)

    # The line at 15 um
    plt.plot([1e9/frequnit, 10e12/frequnit], [15e-6/yunit, 15e-6/yunit], c='k', lw=2)

    #xunit = 1e12
    #def reasonable_ticks(a): 
        #""" Define the grid and ticks a bit denser than by default """
        #x=np.trunc(np.log10(a)); y=a/10**x/10
        #return (10**x, 2*10**x,5*10**x)[int(3*y)]
    #xticks = np.arange(xlim[0], xlim[1], reasonable_ticks((xlim[1]-xlim[0])/10))
    #xnumbers = [("%.2f"%(f/xunit) if abs(f%reasonable_ticks((xlim[1]-xlim[0])/5))<(xunit/1000) else "") for f in xticks]
    #plt.xticks(xticks, xnumbers); 
    plt.minorticks_on();  
    #plt.grid(True)

    plt.xlim((1,3))
    plt.ylim(ylim)
    if logarithmic:
        plt.xscale('log')
        xticks = [500e-3, 1, 2]
        plt.xticks(xticks, xticks)
        plt.yscale('log')
        yticks = [5,10,20]
        plt.yticks(yticks, yticks)

    plt.grid()
    if quantity == 'absNre':
        plt.title('$\\varepsilon = 100$')
    else:
        plt.title('%s' % (quantity.capitalize()))

    a, b  = os.path.split(os.getcwd())
    print a,b
    a, upupdirname  = os.path.split(a)
    print a, upupdirname
    plt.savefig('%s_%s.png' % (quantity, upupdirname), bbox_inches='tight')
