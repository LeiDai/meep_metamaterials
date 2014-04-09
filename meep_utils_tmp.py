#!/usr/bin/env python
#coding:utf8 
"""
Here you can find various functions and classes that facilitate the work with python-meep.
I believe some of these functions ought to be implemented in the meep module.
Filip Dominec 2012-2013
"""
import numpy as np
import os, os.path, sys, subprocess, time
from scipy.constants import c, epsilon_0, mu_0

#import meep
import meep_mpi as meep

## XXX
import _meep_mpi


## Define the simulated models as a class (much simpler and more flexible than handling callbacks)
class AbstractMeepModel(meep.Callback): #{{{
    def __init__(self):
        meep.Callback.__init__(self)
        self.double_vec = None   # (callback function to be redirected to the desired function)
        self.return_value = True  

    def register_local(self, param, val):
        """ Adds a parameter as an attribute of the model (either number or float), adds it also to the simulation name"""
        setattr(self, param, val)
        ## prepare the parameter to be added into name (if not conversible to float, add it as a string)
        try: 
            self.simulation_name += ("_%s=%.2e") % (param, float(val))
            self.parameterstring += "#param %s,%.4e\n" % (param, val)
        except: 
            self.simulation_name += ("_%s=%s") % (param, val)
            self.parameterstring += "#param %s,%s\n" % (param, val)

    def register_locals(self, params):
        """ Scans through the parameters and calls register_local() for each """
        self.parameterstring = ""
        ## First look up for the "VIP" parameters that should come first in the name:
        preferred_params = ['resolution', 'comment', 'frequency', 'simtime']
        for param in preferred_params:
            if params.get(param): 
                val = params.get(param)
                self.register_local(param, val)
        ## Then add all remaining parameters of the model
        for (param, val) in params.iteritems():
            if param != 'self' and param not in preferred_params:
                self.register_local(param, val)

    def eps(self, r):
        """ Scans through materials and adds the high-frequency part of permittivity for each of them. 
        This is why materials should never overlap. """
        for mat in self.materials:
            if mat.where(r): return mat.eps
        else: return 1.

    def build_polarizabilities(self, structure):
        """ 
        This is a helper to define the susceptibilities using the callback.
        It goes through all polarizabilities for all materials. 
        Applicable for time-domain simulation only, because dispersive model is not implemented for 
        frequency-domain simulation yet.
        """
        avail_cbs = [meep.DBL5, meep.DBL4, meep.DBL3, meep.DBL2, meep.DBL1,]
        avail_cb_setters = [meep.set_DBL5_Callback, meep.set_DBL4_Callback, meep.set_DBL3_Callback, 
                meep.set_DBL2_Callback, meep.set_DBL1_Callback,]
        for material in self.materials:
            meep.master_printf("\tAdding material: %s with epsilon: %s at frequency %.4g Hz\n" % 
                    (material.name, analytic_eps(material, self.srcFreq).__str__(), self.srcFreq))
            for polariz in material.pol:
                if avail_cbs == []: meep.master_printf("Error: too many polarizabilities defined (5) Reduce their number.")
                next_cb, next_cb_setter = avail_cbs.pop(), avail_cb_setters.pop()
                self.return_value = polariz['sigma']
                self.double_vec = material.where  ## redirect the double_vec() function callback
                next_cb_setter(self.__disown__())    
                if "lorentzian_susceptibility" in dir(meep):
                    ## for meep 1.2 or newer
                    structure.add_susceptibility(next_cb, meep.E_stuff, 
                            meep.lorentzian_susceptibility(polariz['omega']/c, polariz['gamma']/c))
                    #else:todo: fix in python-meep
                        #print dir(meep)
                        #structure.add_susceptibility(next_cb, meep.E_stuff, 
                                #meep.drude_susceptibility(polariz['omega']/c, polariz['gamma']/c)) 

                else:
                    ## for meep 1.1.1 or older
                    structure.add_polarizability(next_cb, polariz['omega']/c, polariz['gamma']/c)


    def TestMaterials(self):
        """ Call the where() function for each material, in order to make sure there are no errors
        (SWIG callback does not report where the error occured, it just crashes) """
        for material in self.materials: 
            for x in np.linspace(-self.size_x/2, self.size_x/2, 10):
                for y in np.linspace(-self.size_y/2, self.size_y/2, 10):
                    for z in np.linspace(-self.size_z/2, self.size_z/2, 10):
                        if material.where(meep.vec(x, y, z)): print "teotueot"
#}}}
def process_param(args):#{{{
    """ Parse command-line parameters

    Some of them control the simulation, but all remaining will be passed to the model
    """
    sim_param = {   'frequency_domain':False,
                    'frequency':       None,
                    'MaxIter':         5000,
                    'MaxTol':          1e-2,
                    'BiCGStab':        16 }
    model_param = {}
    for namevalue in args:
        name, value = namevalue.split("=")
        if name == "frequency": 
            sim_param['frequency']          = float(value)
            sim_param['frequency_domain' ]  = True
        elif name == "maxtol": MaxTol = float(value)
        elif name == "maxiter": MaxIter = int(value)
        else:           ## all other parameters will be passed to the model:
            try: model_param[name] = float(value)
            except ValueError: model_param[name] = value
    return sim_param, model_param
#}}}

## Geometrical primitives to help defining the geometry
def in_xslab(r,cx,d):#{{{
    return (abs(r.x()-cx) < d/2)
def in_yslab(r,cy,d):
    return (abs(r.y()-cy) < d/2)
def in_zslab(r,cz,d):
    return (abs(r.z()-cz) < d/2)
def in_xcyl(r,cy,cz,rad):
    return ((r.y()-cy)**2+(r.z()-cz)**2) < rad**2
def in_ycyl(r,cx,cz,rad):
    return ((r.x()-cx)**2+(r.z()-cz)**2) < rad**2
def in_zcyl(r,cx,cy,rad):
    return ((r.x()-cx)**2+(r.y()-cy)**2) < rad**2
def in_sphere(r,cx,cy,cz,rad):
    return ((cx-r.x())**2 + (cy-r.y())**2 + (cz-r.z())**2)**.5 < rad
def in_ellipsoid(r,cx,cy,cz,rad):
    return ((cx-r.x())**2 + (cy-r.y())**2 + (cz-r.z())**2)**.5 < rad
#def in_xcone(r,cy,cz,cx,d):
    #return (abs(r.x()-cx) < d/2)
#}}}

## Use the same dispersive materials for time- and frequency-domain simulation
def permittivity2conductivity(complex_eps, freq):#{{{
    """
    Complex permittivity can express also the conductivity of the sample (in the same
    manner as dielectric losses) with the corresponding relation:
        complex_eps = real_eps - 1j conductivity / (frequency * 2*pi * epsilon_0)
    Therefore it should be inverted for classic D-conductivity:
        conductivity = -
    In order to simulate any lossy medium with the freq-domain solver, we invert this relation
    to obtain a (nondispersive) conductivity for one frequency. But it does not give the same results
    as time-domain simulation.

        What we know: 
        function of c, f, 2pi, eps0, eps.im/eps.r
        should give dimension 1  to feed unitless meep
        should be proportional to eps.im/eps.r
        should give ca. 50000 for omega = 2pi * 800 GHz and eps.im/eps.r=0.02
           => should give 2.5e6 for (eps.im/eps.r=1)
        should be proportional to frequency omega
            => should give 5e-7 for omega = 1 and  (eps.im/eps.r=1)
        already was pre-divided for meep by c = 3e8  (which counts here)
            => in real life it shall be  3e-6 * 3e8 = 148
        should be proportional to epsilon0 [C/Vsm], which is similar in units to conductivity
            => if epsilon0 was 1, it would be 1.7e13 -> roughly c**2
    """
    # return complex_eps.imag * freq * 2*np.pi * epsilon_0 * complex_eps.real  ## orig. idea
    # return complex_eps.imag * freq * 2*np.pi * epsilon_0 * complex_eps.real ## also wrong
    # return complex_eps.imag / complex_eps.real * 2*np.pi * c
    #return complex_eps.imag / complex_eps.real * 6.28*freq * 8.85e-12 * c
    magic_constant = 1.65e13       ## A. K. A. bulgarian constant...
    return complex_eps.imag / complex_eps.real * 6.28 * freq / c * 8.85e-12 * magic_constant
#}}}
def analytic_eps(mat, freq):#{{{
    complex_eps = mat.eps
    for polariz in mat.pol:
        complex_eps += polariz['sigma'] * polariz['omega']**2 / (polariz['omega']**2 - freq**2 - 1j*freq*polariz['gamma']) 
    return complex_eps # + sum(0)
#}}}
class MyHiFreqPermittivity(meep.Callback):#{{{
    def __init__(self, model, frequency):
        meep.Callback.__init__(self)
        self.model = model
        self.frequency = frequency
    def double_vec(self, r):
        for material in self.model.materials:
            if material.where(r):
                return analytic_eps(material, self.frequency).real
        else: return 1
#}}}
class MyConductivity(meep.Callback):#{{{
    def __init__(self, model, frequency):
        meep.Callback.__init__(self)
        self.model = model
        self.frequency = frequency
    def double_vec(self, r):
        for material in self.model.materials:
            if material.where(r):
                return permittivity2conductivity(analytic_eps(material, self.frequency), self.frequency)
        else: return 0
#}}}
def plot_eps(to_plot, filename="epsilon.png", plot_conductivity=False, freq_range=(1e9, 1e17), mark_freq=[]):#{{{
    """ Plots complex permittivity of the materials to a PNG file

    Accepts list of materials
    """
    from scipy.constants import epsilon_0
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    frequency = 10**np.arange(np.log10(freq_range[0]), np.log10(freq_range[1]), .001)
    colors = ['#554400', '#004400', '#003366', '#000088', '#440077', 
              '#661100', '#AA8800', '#00AA00', '#0099DD', '#2200FF', 
              '#8800DD', '#BB3300']

    subplotnumber = 2 if plot_conductivity else 1
    for material in list(to_plot):
        plt.subplot(subplotnumber,1,1)
        if colors: color = colors.pop()
        else: color = 'black'
        label = material.name
        eps = analytic_eps(material, frequency)
        #R = abs((1-eps**.5)/(1+eps**.5))**2     ## Intensity reflectivity
        plt.plot(frequency, np.real(eps), color=color, label=label, ls='-')
        plt.plot(frequency, np.imag(eps), color=color, label="", ls='--')
        plt.ylabel(u"solid: Re($\\varepsilon_r$), dashed: Im($\\varepsilon_r$) ")
        #plt.ylabel(u"Intensity reflectivity")
        for mfreq in mark_freq: plt.plot([mfreq,mfreq], [-1,1]) 

        plt.yscale('symlog')
        plt.xscale('log')
        plt.legend(); 
        #plt.ylim(ymin=1e-2); 
        plt.grid(True)

        if plot_conductivity:
            plt.subplot(subplotnumber,1,2)
            label = ""
            cond = eps * frequency * epsilon_0 / 1j
            plt.plot(frequency, np.real(cond), color=color, label=label, ls='-')
            plt.plot(frequency, np.imag(cond), color=color, label="", ls='--')
            plt.ylabel(u"$Re\\sigma$ (solid), Im$\\sigma$ (dashed)")
            plt.yscale('symlog'); plt.xscale('log'); plt.legend(); plt.grid(True)

    ## Finish the graph 
    plt.xlabel(u"Frequency [Hz]") 
    #plt.xlim((, ))
    plt.savefig(filename, bbox_inches='tight')
#}}}
def lorentzian_unstable_check_new(model, dt, quit_on_warning=True): #{{{
    for mat in model.materials:
        eps_ts = analytic_eps(mat, 1/dt/np.pi)  
        if np.real(eps_ts)<0:
            meep.master_printf("Warning: for material '%s', the permittivity is negative at timestepping frequency eps(1/pi/dt)=eps(%g)=%s.\n" % \
                        (mat.name, 1/dt/np.pi, eps_ts.__str__()));
            if quit_on_warning: quit()
        for pol in mat.pol:
            omega_0, gamma = pol['omega'], pol['gamma']
            if (omega_0 > gamma/2):
                z2 = np.sqrt(gamma*gamma + 4*omega_0*omega_0)/2
            else:
                z2 = gamma/2 + np.sqrt(gamma*gamma - 4*omega_0*omega_0)/2
            if z2 > 1/dt/np.pi:
                meep.master_printf("Warning: for material '%s', the oscillator pole having magnitude |z|=%g will be probably unstable when 1/pi/dt=%g.\n" % \
                        (mat.name, z2, 1/dt/np.pi));
                if quit_on_warning: quit()
#}}}

## Useful to make 3D snapshots of the fields
def run_bash(cmd, anyprocess=False): #{{{
    if meep.my_rank() == 0 or anyprocess:
        #meep.master_printf("CMD: "  + cmd+ "\n")
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out = p.stdout.read().strip()
        return out
#}}}
class SnapshotMaker(): #{{{
    """
    Saves the field vectors to a HDF5/VTK file at a given time.

    This function unfortunately requires command-line tools from Linux system.
    """
    def __init__(self, field=None, snapshot_times=[], outputdir=None, volume=None):
        """ Remember the time when to take a snapshot. A list of two or more numbers may be provided. """
        self.field = field
        self.snapshot_times = snapshot_times
        self.outputdir = outputdir
        self.volume = volume
        if not os.path.exists(self.outputdir): run_bash("mkdir -p %s" % self.outputdir)
        meep.master_printf("Registered %d time points(s) to take snapshot\n" % len(self.snapshot_times))

    def poll(self, now=None):
        """ Check if the time has come to make a snapshot.  """
        if len(self.snapshot_times)>0 and now > self.snapshot_times[0]:
            self.take_snapshot(now)
            self.snapshot_times[0:1] = [] ## Remove the current snapshot now from the list of remaining times
    def take_snapshot(self, now):
        """
        This function creates a directory, exports the dielectric structure, 
        then exports the electric and magnetic field vectors to separate HDF5 files. 
        Finally it converts the fields to VTK format (to be readable e. g. by Mayavi2).
        Note that it is not only called by self.poll(), but I may be also called at will.
        """
        meep.master_printf(" * Saving field snapshot at time=%.4e... \n" % now)

        ## Export the dielectric structure so that we can visually verify it
        structure_filename = os.path.join(self.outputdir, 'structure') ## name of the file to save eps (without extension)
        #if meep.my_rank() == 0 and not os.path.exists(self.outputdir): os.mkdir(self.outputdir)
        if not os.path.exists(self.outputdir): run_bash("mkdir -p %s" % self.outputdir, anyprocess=True) ## NO DATA
        #if not os.path.exists(structure_filename):

        #if meep.my_rank() ==0 :
        #if not "outputfile" in locals():
        outputfile = meep.prepareHDF5File(structure_filename+'.h5') 
        self.field.output_hdf5(meep.Dielectric, self.volume, outputfile) 
        del(outputfile)

        ## Export the fields snapshot
        snapshotfiles = []
        for (component, compname) in [(meep.Ex, "Ex"), (meep.Ey, "Ey"), (meep.Ez, "Ez"), 
                (meep.Hx, "Hx"), (meep.Hy, "Hy"), (meep.Hz, "Hz")]:
            snapshotfiles.append("%s/snapshot_%s_t%e.h5" % (self.outputdir, compname, now))
            snapshotfile = meep.prepareHDF5File(snapshotfiles[-1])
            self.field.output_hdf5(component, self.volume, snapshotfile, 1) 
            del(snapshotfile) 
        ## Convert the files for Mayavi2; join E-fields and H-fields components into vector fields
        run_bash("h5tovtk %s.h5 -o %s.vtk" % (structure_filename, structure_filename))
        run_bash("h5tovtk %s -t 0 -o %s/t%0.3e_Evec.vtk" % (" ".join(snapshotfiles[:3]), self.outputdir, now)) ## todo: use path.join
        run_bash("h5tovtk %s -t 0 -o %s/t%0.3e_Hvec.vtk" % (" ".join(snapshotfiles[3:]), self.outputdir, now))
        # (Note: the -t parameter selects real part only)
#}}}
class SliceMaker(): #{{{
    """
    Saves the field vectors from a slice to a HDF5 file during the simulation. The slice may be specified either
    by two-dimensional meep.volume() object, or by a normal ("x", "y" or "z") and the slice position on this normal.

    After simulation ends, the data are exported as GIF animation and 3D VTK file.
    Optionally, you may disable these output formats and/or enable output to many many PNG files or to a 3D HDF file.

    This function unfortunately requires command-line tools from Linux system.
    """
    def __init__(self, field=None, component=meep.Ex, timebounds=(0, np.inf), timestep=0, 
            volume=None, normal=None, position=None, model=None, pad=0,
            outputdir="", name=None, outputPNGs=False, outputGIF=True, outputHDF=False, outputVTK=False):
        """  """
        self.field = field
        self.outputdir = outputdir
        self.component = component
        self.timebounds = timebounds
        self.timestep = timestep

        self.outputPNGs = outputPNGs
        self.outputGIF  = outputGIF
        self.outputHDF  = outputHDF
        self.outputVTK =  outputVTK

        if volume:
            self.volume = volume
            meep.master_printf("Will record slices at times %.3g, %.3g ... %.3g s \n" % (timebounds[0], timebounds[0]+timestep, timebounds[1]))
        else:
            #if not position: 
                #raise RuntimeError("Specify the position of the cut plane (on the axis perpendicular to it)")
            if normal=="x":
                self.volume = meep.volume(
                        meep.vec(position, -model.size_y/2+pad, -model.size_z/2+pad), 
                        meep.vec(position, model.size_y/2-pad, model.size_z/2-pad)) 
            elif normal=="y":
                self.volume = meep.volume(
                        meep.vec(-model.size_x/2+pad, position, -model.size_z/2+pad), 
                        meep.vec(model.size_x/2-pad, position, model.size_z/2-pad)) 
            elif normal=="z":
                self.volume = meep.volume(
                        meep.vec(-model.size_x/2+pad, -model.size_y/2+pad, position), 
                        meep.vec(model.size_x/2-pad, model.size_y/2-pad, position)) 
            #else: 
                #print normal
                #raise RuntimeError("Specify the normal parameter as 'x', 'y' or 'z'")
            meep.master_printf("Will record slices at %s=%.4f m, at times %g, %g ... %g s \n" \
                    % (normal, position, timebounds[0], timebounds[0]+timestep, timebounds[1]))

        self.outputdir = outputdir
        if not name: self.name = "Slice_%s%.4f" % (normal, position)
        else: self.name = name
        self.images_number = 0
        self.last_slice_time = 0.
        #slices=[]
        #slices.append({"name":"Ex_xz_slice",    "component":meep.Ex,                        "geom":

        if not os.path.exists(outputdir): run_bash("mkdir -p %s" % outputdir, anyprocess=True)
        self.openfile = meep.prepareHDF5File("%s.h5" % (os.path.join(self.outputdir, self.name)))

    def poll(self, now):
        """ Check if the now has come to add a new slice  """
        if (now-self.last_slice_time > self.timestep) and (now > self.timebounds[0]) and (now < self.timebounds[1]):
            self.images_number += 1 
            self.field.output_hdf5(self.component, self.volume, self.openfile, 1) 
            self.last_slice_time = now

    def finalize(self):
        #run_bash("cd %s; rm *png" % outputdir)
        meep.master_printf("\n\nImages to gif\n")
        del(self.openfile)
        run_bash("cd %s; h5topng -t 0:%d -R -Zc dkbluered -a yarg %s.h5 -S 1" % (self.outputdir, self.images_number-1, self.name))
        if self.outputPNGs or self.outputGIF: run_bash("cd %s; convert -compress None -delay 10 *png %s.gif" % (self.outputdir, self.name))
        if not self.outputPNGs: 
            run_bash("cd %s; rm %s*.png" % (self.outputdir, self.name))
        if self.outputVTK: run_bash("cd %s; h5tovtk %s.h5 -o %s.vtk" % (self.outputdir, self.name))
        if not self.outputHDF: run_bash("cd %s; rm %s.h5" % (self.outputdir, self.name))
#}}}

## Print the progress and estimated time
class Timer():#{{{
    def __init__(self, simtime):
        self.starttime = time.time()
        self.simtime = simtime
        meep.master_printf("\tSimulation time: %e [s] = %e time units\n" % (simtime, simtime*c))
        self.last_reported_time = 0
    def get_time(self):
        return time.time()-self.starttime
    def print_progress(self, now):
        if now > 0 and (now-self.last_reported_time > self.simtime*.1): 
            meep.master_printf("Progress %.2f of expected total %d s\n" % (now / self.simtime, (self.simtime / now * self.get_time())))
            self.last_reported_time = now
#}}}
def notify(title, run_time=None):#{{{
    """
    Shows a bubble with notification that your results are about to be ready!
    Requires python-notify installed, otherwise just quits

    Note: you may also use similar call with libnotify-bin from your bash scripts:
        run_bash('notify-send -t 3000 -i "face-glasses" "MEEP simulation finished %s" "%s"' % (timestring, title))
    """
    if meep.my_rank() != 0: return
    try: 
        if run_time: timestring = "in %d s" % int(run_time)
        else: timestring = ""
        import pynotify
        pynotify.init("image")
        n = pynotify.Notification("MEEP simulation finished %s" % (timestring), title, "face-glasses")
        n.show()
    except:
        pass
#}}}

## Obtain and process the s-parameters of the structure (works for both time- and freq-domain)
def get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
        frequency_domain=False, frequency=None, pad_zeros=0.0, maxf=np.inf, minf=0, Kx=0, Ky=0):#{{{
    """ Returns the frequency, s11 (reflection) and s12 (transmission) in frequency domain """
    ## TODO allow omitting second monitor (-> returns s12=None)


    t, Ex1 = monitor1_Ex.get_waveforms()
    t, Hy1 = monitor1_Hy.get_waveforms()
    t, Ex2 = monitor2_Ex.get_waveforms()
    t, Hy2 = monitor2_Hy.get_waveforms()

    try:
        import matplotlib
        from matplotlib import pyplot as plt
        plt.figure(figsize=(15,10))
        plt.plot(t, abs(Ex1), label="Ex1")
        plt.plot(t, abs(Hy1), label="Hy1")
        plt.plot(t, abs(Ex2), label="Ex2")
        plt.plot(t, abs(Hy2), label="Hy2")
        #plt.ylim(1, 1e6)
        plt.legend()
        plt.yscale("log")
        plt.savefig("timedomain_debug.png")
    except:
        print "Timedomain plot failed"

    ## Obtain the electric and magnetic fields spectra
    if frequency_domain:            ## No need for FFT in frequency domain, just copy the value
        freq = np.array([frequency])
        (Ex1f, Hy1f, Ex2f, Hy2f) = (Ex1, Hy1, Ex2, Hy2)
    else:
        ## Optionally extend the data range by zeros (for stable eff param retrieval)
        def pad(arr, length):
            appended = np.zeros(int(length*len(arr))) #  * np.exp(- np.arange(0, len(,
            return np.append(arr, appended)
        if pad_zeros: 
            Ex1, Hy1, Ex2, Hy2  =  map(lambda x: pad(x, pad_zeros), (Ex1, Hy1, Ex2, Hy2))

        ## Calculate the Fourier transform of the recorded time-domain waves
        numpoints = len(Ex1)
        #fftfreq(signal.size, Sample spacing.[d])	Return the Discrete Fourier Transform sample frequencies.
        freq = np.arange(0., int(numpoints/2)) / numpoints / (t[1]-t[0]) # take positive frequency range only
        ## TODO Positive frequency range should be  separated just by truncation below, as 'minf => 0'

        #fftshift(x[, axes])	Shift the zero-frequency component to the center of the spectrum.
        (Ex1f, Hy1f, Ex2f, Hy2f) = map(lambda x: np.fft.fft(np.real(x))[0:int(numpoints/2)], (Ex1, Hy1, Ex2, Hy2))
        


        ## Truncate the data ranges to allowed radiating angles, and possibly to minf<freq<maxf
        truncated = np.logical_and((Ky**2+Kx**2)<(2*np.pi*freq/c)**2, freq>minf, freq<maxf)
        (Ex1f, Hy1f, Ex2f, Hy2f, freq) = map(lambda x: x[truncated], (Ex1f, Hy1f, Ex2f, Hy2f, freq))

    ## Prepare the angles at which the wave propagates (dependent on frequency, Kx and Ky)
    #rho = np.arcsin(Kx / (2*np.pi*freq/c)) TODO
    beta0 = np.arcsin((Kx**2+Ky**2)**.5 / (2*np.pi*freq/c))
    print 'beta0', beta0

    ## Separate the forward and backward wave in frequency domain 
    ##    (Efield+Hfield)/2 ->    forward wave amplitude, 
    ##    (Efield-Hfield)/2 ->    backward wave amplitude
    #in1, out1 =  (Ex1f+Hy1f)/2, (Ex1f-Hy1f)/2 ## old: works only for perp. incidence beta0=0
    #in2, out2 =  (Ex2f-Hy2f)/2, (Ex2f+Hy2f)/2
    in1, out1 =  (Ex1f+Hy1f/np.cos(beta0))/2, (Ex1f-Hy1f/np.cos(beta0))/2
    in2, out2 =  (Ex2f-Hy2f/np.cos(beta0))/2, (Ex2f+Hy2f/np.cos(beta0))/2
    ## Todo optimize cos(arcsin x ) = sqrt(1-x**2)

    ## TEMPORARY Plot spectral profile
    import matplotlib
    from matplotlib import pyplot as plt
    plt.figure(figsize=(15,10))
    plt.plot(freq, abs(in1), label="in1")
    plt.plot(freq, abs(out1), label="out1")
    plt.plot(freq, abs(in2), label="in2")
    plt.plot(freq, abs(out2), label="out2")
    plt.xlim(0, 1e12)
    plt.legend()
    plt.yscale("log")
    plt.savefig("ampli_debug_band.png")


    ## Check if PML works well (optional)
    #meep.master_printf("PML reflection max: %.4e" % max(in2 / in1))

    ## Get the s-parameters 
    s11 = out1 / in1
    s12 = out2 / in1

    ## Return the S-parameters (i. e. complex reflection and transmission)
    return freq, s11, s12
#}}}
def get_phase(complex_data):#{{{
    """ Unwraps and shifts the phase from Fourier transformation """
    if len(complex_data) <= 1: return np.angle(complex_data)
    phase = np.unwrap(np.angle(complex_data))
    center_phase = phase[min(5, len(phase)-1)] ## 5 is chosen to avoid zero freq.
    return phase-(round(center_phase/2/np.pi)*2*np.pi)
#}}}
class AmplitudeMonitorPlane():#{{{
    """ Calculates an average of electric field and perpendicular magnetic field.

    I asked for a similar field-averaging function built in MEEP, but it seems not to be implemented yet.
    http://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg04447.html

    Note this implementation requires the planes are in vacuum (where impedance = 1.0)
    """
    def __init__(self, comp=None, size_x=None, size_y=None, z_position=None, Kx=0, Ky=0):
        self.comp=comp
        self.size_x = size_x
        self.size_y = size_y
        self.z_position = z_position
        self.Kx = Kx
        self.Ky = Ky

        self.t = []
        self.waveform = []

    def average_field(self, field):
        """
        Average field component in some plane, return amplitudes 
        This function is ineffective - it should be implemented in C++ in meep itself

        5x5 grid is usually optimal (no visible difference between 10x10 grid and 5x5 grid)

        TODO:  This class implements a workaround for unavailable amplitude averaging in python-meep.
        In this implementation, the geometrical averaging is ineffective and inflexible here. 
        """
        xcount, ycount = (1,1)
        field_sum = 0 
        # The mode function has the form of an oblique plane wave
        #for x in [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]:
            #for y in [y0*self.size_y/ycount+(self.size_y/2/ycount)-self.size_y/2 for y0 in range(ycount)]:
                #field_sum += (field.get_field(self.comp, meep.vec(x, y, self.z_position)) *  
                                #np.exp(1j*(self.Kx*x + self.Ky*y)) )
        #return field_sum/(xcount*ycount)


        return field.get_field(self.comp, meep.vec(0,0, self.z_position))
        ## New way (removes explicit cycle, few percent faster)
        xr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        yr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        xm, ym = np.meshgrid(xr,yr)
        points = zip(xm.flatten(), ym.flatten())
        sum_ = sum(map(lambda pos: field.get_field(self.comp, meep.vec(pos[0], pos[1], self.z_position)), points))
        return sum_/(xcount*ycount)

        ## Yet newer way
        #xr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        #yr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        #xm, ym = np.meshgrid(xr,yr)
        #v = meep.vec(0,0, self.z_position)
        #points = zip(xm.flatten(), ym.flatten())
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, self.comp, v), points))
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, self.comp, _meep_mpi.new_vec(0,0,0)), points))
        #cp =self.comp
        #zp = self.z_position
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, cp, _meep_mpi.new_vec(pos[0], pos[1], zp)), points))
        #return sum_/(xcount*ycount)
        


    #def NEW_average_field_xy_plane(field, component, zpos, model): ## TODO use the internal integration of fields by MEEP
        # TODO rewrite:
        # (fields &f, linear_integrand_data &d, const volume &v, component cgrid)
        # f.integrate(0, 0, linear_integrand, (void *) &d, v)
        #integrate(meep::fields *,int,meep::component const *,meep::field_function,void *,meep::volume const &,double *)
        #v = meep.volume(
                #meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness), 
                #meep.vec(model.size_x/2, model.size_y/2, -model.size_z/2+model.pml_thickness))
        #return field.integrate(1, [component], meep.one_vec, [], v)
      #Possible C/C++ prototypes are:
        #meep::fields::integrate(int,meep::component const *,meep::field_function,void *,meep::volume const &,double *)
        #meep::fields::integrate(int,meep::component const *,meep::field_function,void *,meep::volume const &)
        #meep::fields::integrate(int,meep::component const *,meep::field_rfunction,void *,meep::volume const &,double *)
        #meep::fields::integrate(int,meep::component const *,meep::field_rfunction,void *,meep::volume const &)
        #fields::integrate(int num_fvals, const component *components,
				  #field_function integrand,
				  #void *integrand_data_,
				  #const volume &where,
				  #double *maxabs
    
    def record(self, field=None):
        """ 
        Useful for time-domain simulation only
        """
        self.t.append(field.time()/c)
        self.waveform.append(self.average_field(field))

    def get_waveforms(self):
        """ Return the recorded waveform (in time domain) """
        if len(self.t) <= 1:
            t, result_wform = np.array(self.t), np.array(self.waveform)
        else:
            t = np.array(self.t[:-1])
            ## The FDTD calculation introduces half-step time shift between Ex and Hy. Compensated by averaging the Hy field
            ## with its value in a next timestep. The error is reduced from O1 to O2.
            ## See http://ab-initio.mit.edu/wiki/index.php/Synchronizing_the_magnetic_and_electric_fields
            if meep.is_magnetic(self.comp) or meep.is_B(self.comp):
                result_wform = np.array(self.waveform[:-1])/2. + np.array(self.waveform[1:])/2.
            else: 
                result_wform = np.array(self.waveform[:-1])
            
        return t, result_wform 
         ## time, 
        ## TODO this will have to be modified in order to account for oblique incidence
        ## TODO take into account the medium impedance (... divide the Hfield)
#}}}
class AmplitudeMonitorPoint(AmplitudeMonitorPlane):#{{{
    """ Calculates an average of electric field and perpendicular magnetic field.

    I asked for a similar field-averaging function built in MEEP, but it seems not to be implemented yet.
    http://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg04447.html

    Note this implementation requires the planes are in vacuum (where impedance = 1.0)
    """
    def __init__(self, Ecomp=None, Hcomp=None, pos=None):
        self.Ecomp=Ecomp
        self.Hcomp=Hcomp
        self.pos = pos          ## where to record the field
        self.t = []
        self.Efield = []
        self.Hfield = []

    def get_amplitude(self, field, component):
        """ Record field in some point. No averaging here, but (inherited) time recording is pretty useful. """
        #count = 5
        #field_sum = 0 
        #for x in [x0*self.size_x/count+(self.size_x/2/count)-self.size_x/2 for x0 in range(count)]:
            #for y in [y0*self.size_y/count+(self.size_y/2/count)-self.size_y/2 for y0 in range(count)]:
                #field_sum += field.get_field(component, meep.vec(x, y, self.z_position))
        return field.get_field(component, self.pos)
#}}}

def savetxt(name="output.dat", freq=None, s11=None, s12=None, model=None, polar_notation=True): #{{{
    """ Saves the s-parameters to a file including comments """
    with open(model.simulation_name+".dat", "w") as outfile: 
        outfile.write("#Parameters Parameters\n")
        outfile.write("#param layer_thickness[m],%.6e\n" % (model.monitor_z2 - model.monitor_z1))
        if "interesting_frequencies" in dir(model): interest_freq = model.interesting_frequencies 
        else: interest_freq = (0, model.srcFreq+model.srcWidth)
        outfile.write("#param plot_freq_min[Hz],%.3e\n" % interest_freq[0])
        outfile.write("#param plot_freq_max[Hz],%.3e\n" % interest_freq[1])
        outfile.write("#param simulation_orig_name,%s\n" % model.simulation_name)
        outfile.write(model.parameterstring)
        if polar_notation:
            ## Convert to polar notation
            s11amp, s12amp, s11phase, s12phase = abs(s11), abs(s12), get_phase(s11), get_phase(s12)
            ## Save polar
            outfile.write("#x-column Frequency [Hz]\n#Column Reflection amplitude\n#Column Reflection phase\n" + \
                        "#Column Transmission amplitude\n#Column Transmission phase\n")
            np.savetxt(outfile, zip(freq, s11amp, s11phase, s12amp, s12phase), fmt="%.8e")
        else:
            ## Save cartesian
            # TODO should save in such format that PKGraph understands complex values
            outfile.write("#x-column Frequency [Hz]\n#Column Reflection Re\n#Column Reflection Im\n" + \
                        "#Column Transmission Re\n#Column Transmission Im\n")
            np.savetxt(outfile, zip(freq, s11.real, s11.imag, s12.real, s12.imag), fmt="%.8e")

#}}}
def get_simulation_name(argindex=1): #{{{
    """Get the name of the last simulation run.

    Priority: 1) parameter, 2) last_simulation_name.txt, 3) working directory"""
    cwd = os.getcwd()
    if len(sys.argv)>argindex and sys.argv[argindex] != "-"  and __name__ == "__main__": 
        print "Parameter passed:", sys.argv[argindex]
        last_simulation_name = sys.argv[argindex]
    elif os.path.exists(os.path.join(cwd, 'last_simulation_name.txt')):
        print "Loading from", os.path.join(cwd, 'last_simulation_name.txt')
        last_simulation_name = os.path.join(cwd, open(os.path.join(cwd, 'last_simulation_name.txt'),'r').read().strip())
    else:
        print "Error: No input file provided and 'last_simulation_name.txt' not found!"
        last_simulation_name = cwd
    if (last_simulation_name[-4:] == ".dat"): last_simulation_name = last_simulation_name[:-4] # strip the .dat extension
    return  last_simulation_name
#}}}
def load_rt(filename, layer_thickness=None, plot_freq_min=None, plot_freq_max=None, truncate=True): #{{{
    """ Loads the reflection and transmission spectra and simulation settings 

    Returns:
    * frequency axis
    * reflection s11 and transmission s12 as complex np arrays

    Compatible with the PKGraph text data file with polar data: 
    * parameters in header like: #param name,value
    * column identification like: #column Ydata
    * data columns in ascii separated by space
    Expects polar data with columns: frequency, s11 ampli, s11 phase, s12 ampli, s12 phase
    """
    with open(filename+'.dat') as datafile:
        for line in datafile:
            if line[0:1] in "0123456789": break         # end of file header
            value = line.replace(",", " ").split()[-1]  # the value of the parameter will be separated by space or comma
            if ("layer_thickness" in line) and (layer_thickness == None): d = float(value)
            if ("plot_freq_min" in line) and (plot_freq_min == None): plot_freq_min = float(value)
            if ("plot_freq_max" in line) and (plot_freq_max == None): plot_freq_max = float(value)
    xlim = (plot_freq_min, plot_freq_max)
    (freq, s11amp, s11phase, s12amp, s12phase) = \
            map(lambda a: np.array(a, ndmin=1), np.loadtxt(filename+".dat", unpack=True)) 

    ## Limit the frequency range to what will be plotted (recommended)
    if truncate:
        (d0,d1) = np.interp((plot_freq_min, plot_freq_max), freq, range(len(freq)))
        (freq, s11amp, s11phase, s12amp, s12phase) = \
                map(lambda a: a[int(d0):int(d1)], (freq, s11amp, s11phase, s12amp, s12phase))
    return freq, s11amp, s11phase, s12amp, s12phase, xlim, (d, plot_freq_min, plot_freq_min)
#}}}

