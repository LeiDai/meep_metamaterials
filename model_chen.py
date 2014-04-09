#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab
"""
"""

class ChenSRR_model(meep_utils.AbstractMeepModel): #{{{
    """ Connected metallic double split-rings on tunable dielectric substrate
    FD 2014-02-20
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=3e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            spacing=100e-6, monzd=50e-6, epsilon=100):

        self.DSOthick = 1000e-6
        monzd = self.DSOthick + 60e-6
        self.monzd = monzd

        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "ChenSRR"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Initialization of materials used
        if 'TiO2' in comment:           self.materials = [meep_materials.material_TiO2_THz(where = self.where_diel)]
        elif 'DielLossless' in comment: self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=epsilon, loss=0.0)]
        elif 'Diel' in comment:         self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=epsilon, loss=0.5)]
        elif 'Metalonly' in comment:    self.materials = []
        else:         
            self.materials =  [meep_materials.material_dielectric(where = self.where_diel, eps=epsilon, loss=0.5)]
            self.materials += [meep_materials.material_dielectric(where = self.where_DSO, eps=20,      loss=0.05)]
        self.materials += [meep_materials.material_Metal_THz(where = self.where_metal)]

        ## Dimension constants for the simulation
        self.pml_thickness = 10e-6
        self.size_x = spacing*2
        self.size_y = spacing
        self.size_z = cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 2000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 1000e9) 


        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14), plot_conductivity=True)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_metal(self, r):



        for cellz in self.cell_centers():

            #if (in_yslab(r, cy=self.resolution/2, d=self.resolution*2)): 
                #h = 50e-6
                #g = -50e-6
                #if (in_ycyl(r, cx=h+g, cz=0, rad=30e-6) and not in_ycyl(r, cx=h+g, cz=0, rad=25e-6)): ## outer ring upper
                    #if ((r.x()>0) or not in_zslab(r, cz=0, d=25e-6)):       ## the notch in outer
                        #return self.return_value
                #if (in_ycyl(r, cx=h+g, cz=0, rad=20e-6) and not in_ycyl(r, cx=h+g, cz=0, rad=15e-6)): ## outer ring upper
                    #if ((r.x()<0) or not in_zslab(r, cz=0, d=15e-6)):       ## the notch in outer
                        #return self.return_value

            if (in_zslab(r, cz=self.resolution*2 + self.DSOthick/2, d=self.resolution*2) and not in_zslab(r, cz=0, d=self.DSOthick)): 
                h = 50e-6
                g = -10e-6
#
                if (in_zcyl(r, cx=h+g, cy=0, rad=30e-6) and not in_zcyl(r, cx=h+g, cy=0, rad=25e-6)): ## outer ring upper
                    if ((abs(r.x()-g)>h) or not in_yslab(r, cy=0, d=125e-6)):       ## the notch in outer
                        return self.return_value
                if (in_zcyl(r, cx=h+g, cy=0, rad=20e-6) and not in_zcyl(r, cx=h+g, cy=0, rad=15e-6)): ## inner ring
                    if ((abs(r.x()-g)<h) or not in_yslab(r, cy=0, d=150e-6)):       ## the notch in inner
                        return self.return_value
 
                if (in_zcyl(r, cx=-h+g, cy=0, rad=30e-6) and not in_zcyl(r, cx=-h+g, cy=0, rad=25e-6)): ## outer ring lower
                    if ((abs(r.x()-g)>h) or not in_yslab(r, cy=0, d=125e-6)):       ## the notch in outer
                        return self.return_value
                if (in_zcyl(r, cx=-h+g, cy=0, rad=20e-6) and not in_zcyl(r, cx=-h+g, cy=0, rad=15e-6)): ## inner ring
                    if ((abs(r.x()-g)<h) or not in_yslab(r, cy=0, d=150e-6)):       ## the notch in inner
                        return self.return_value

                if ((abs(r.x()-g)>(h+25e-6)) and in_yslab(r, cy=0, d=5e-6)):       ## the connection to outer SRRs
                        return self.return_value
                if ((abs(r.x()-g)<(h-15e-6)) and in_yslab(r, cy=0, d=5e-6)):       ## the connection to inner SRRs
                        return self.return_value
#
                if (in_xslab(r, cx=g, d=10e-6)):       ## the horiz connection
                        return self.return_value
                if (in_xslab(r, cx=2*h+g, d=10e-6) or in_xslab(r, cx=-2*h+g, d=10e-6)):       ## the horiz connection
                        return self.return_value
        return 0
    def where_DSO(self, r):
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=self.DSOthick/2, d=self.resolution*4)): 
                return 0
            if (in_zslab(r, cz=0, d=self.DSOthick)): 
                return self.return_value
        return 0
    def where_diel(self, r):


        for cellz in self.cell_centers():
            if (in_zslab(r, cz=self.DSOthick/2, d=self.resolution*4)): 
                return self.return_value
        return 0
#}}}

class ChenI_model(meep_utils.AbstractMeepModel): #{{{
    """ Connected metallic "I" dipole resonators on tunable dielectric substrate
    FD 2014-02-20
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=20e-6,
            monzd=60e-6, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        if 'TiO2' in comment:           self.materials = [meep_materials.material_TiO2_THz(where = self.where_diel)]
        elif 'DielLossless' in comment: self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=epsilon, loss=0.0)]
        elif 'Diel' in comment:         self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=epsilon, loss=0.05)]
        else:                           self.materials = []
        self.materials += [meep_materials.material_Metal_THz(where = self.where_metal)]

        ## Dimension constants for the simulation
        self.pml_thickness = 10e-6
        self.size_x = 40e-6
        self.size_y = 30e-6
        self.size_z = 80e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 2000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 1000e9) 

        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14), plot_conductivity=True)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_metal(self, r):
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=self.resolution/2, d=self.resolution*4) and not in_zslab(r, cz=20e-6, d=40e-6)): 
                h = 50e-6
                g = 00e-6
                if (in_xslab(r, cx=g, d=36e-6) and in_yslab(r, cy=0, d=4e-6)):       ## the vertical conductor of I-resonator
                        return self.return_value
                if (in_xslab(r, cx=g, d=4e-6)):                                     ## the horiz inconnection of resonatorsn
                        return self.return_value
                if (in_yslab(r, cy=0, d=10e-6)): 
                    if (in_xslab(r, cx=g+(36e-6-4e-6)/2, d=4e-6) or in_xslab(r, cx=g-(36e-6-4e-6)/2, d=4e-6)):       ## the horiz capacitor pads
                        return self.return_value
        return 0
    def where_diel(self, r):
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=20e-6, d=40e-6)): 
                return self.return_value
            #if (in_zslab(r, cz=cellz, d=self.radius)): # and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                #return self.return_value
            #if abs(cellz)>self.yspacing*1.5:
                #if (in_zslab(r, cz=cellz, d=self.radius)): # and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                    #return self.return_value
            #else:
                #if (in_zslab(r, cz=cellz, d=self.radius*2)): # and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                    #return self.return_value
        return 0
#}}}


            ## SINGLE SRR
            #if (in_yslab(r, cy=self.resolution/2, d=self.resolution*2)): 
                #h = 50e-6
                #g = -50e-6
                #if (in_ycyl(r, cx=h+g, cz=0, rad=30e-6) and not in_ycyl(r, cx=h+g, cz=0, rad=20e-6)): ## outer ring upper
                    #if ((r.x()>0) or not in_zslab(r, cz=0, d=25e-6)):       ## the notch in outer
                        #return self.return_value


    #def where_metal(self, r):
        #for cellz in self.cell_centers():
            #if (in_zslab(r, cz=self.resolution/2, d=self.resolution*2) and not in_zslab(r, cz=20e-6, d=40e-6)): 
                #if (in_zcyl(r, cx=0, cy=0, rad=30e-6) and not in_zcyl(r, cx=0, cy=0, rad=25e-6)): ## outer ring
                    #if ((r.x()<0) or not in_yslab(r, cy=0, d=25e-6)):       ## the notch in outer
                        #return self.return_value
                #if (in_zcyl(r, cx=0, cy=0, rad=20e-6) and not in_zcyl(r, cx=0, cy=0, rad=15e-6)): ## inner ring
                    #if ((r.x()>0) or not in_yslab(r, cy=0, d=15e-6)):       ## the notch in inner
                        #return self.return_value
                #if ((r.x()<-25e-6) and in_yslab(r, cy=0, d=5e-6)):       ## the connection to outer
                        #return self.return_value
                #if ((r.x()>25e-6) and in_yslab(r, cy=0, d=5e-6)):       ## the connection to inner
                        #return self.return_value
        #return 0
