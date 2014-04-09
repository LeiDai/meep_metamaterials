#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab
"""
"""

class XRectWireFN_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=50e-15, resolution=3e-9, cells=1, monzc=0e-9, Kx=0, Ky=0, padding=25e-9,
            thick=4e-9, width=4e-9, xspacing=50e-9, monzd=100e-9):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        if 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=20, loss=0.001)]
        else:
            self.materials = [meep_materials.material_Au(where = self.where_wire)]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 10e-9
        self.size_x = xspacing
        self.size_y = resolution/1.8
        self.size_z = 100e-9+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e12     
        self.srcFreq = 4e12 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 1000e12) 

        #meep_utils.plot_eps(self.materials)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_wire(self, r):
        for cellz in self.cell_centers():
            #if (in_xcyl(r, cy=self.resolution/4, cz=cellz, rad=self.radius)): 
                #return self.return_value
            if (in_zslab(r, cz=cellz, d=self.thick) and in_xslab(r, cx=self.resolution/4, d=self.width)):  # XXX
                return self.return_value
        return 0
#}}}

class PKCutSheet_model_test(meep_utils.AbstractMeepModel): #{{{
    """  Thin high-permittivity dielectric sheet along E,k vectors
    FD 2013-06-04

         |PMLPMLPML| <------  Absorbing boundaries on X, Y faces
    ^    |         |
    |    +mmmmmmmmm+ <-- monitor plane 2
    z    |         | <------  Bloch-periodic boundaries on X, Y faces
         |         |    
         |XXXXXXXXX| <---sheet of micromachined STO
         |         |
         |         |
    |    +mmmmmmmmm+ <-- monitor plane 1 
    x    |         |
    -->  +sssssssss+ <-- source
         |PMLPMLPML|                 
    The purpose of computing is insight, not numbers. -- Richard W. Hamming

    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=1000e-12, resolution=5e-6, cells=1, monzc=0e-6, Kx=0, Ky=0,
            xs=50e-6, ys=50e-6, zs=50e-6, xyspacing=100e-6, monzd=100e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "PKCTest"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2), 
                ]

        ## Dimension constants for the simulation
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 100e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc, (monzd*cells)/2+monzc)  
        self.simtime = simtime      # [s]
        self.srcWidth = 1000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    #def where_metal(self, r):
        #for cellz in self.cell_centers():
            #if (in_xslab(r, cx=0e-6, d=self.wtth) or in_yslab(r, cy=0e-6, d=self.wtth)) and \
                    #in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                #return self.return_value
        #return 0

    def where_TiO2(self, r):
        for cellz in self.cell_centers():
            zz = in_zslab(r, cz=cellz,    d=self.zs)  
            yy = in_yslab(r, cy=0,          d=self.ys)
            xx = in_xslab(r, cx=0,          d=self.xs)
            if (zz and yy and xx): 
                return self.return_value
        return 0

#}}}
