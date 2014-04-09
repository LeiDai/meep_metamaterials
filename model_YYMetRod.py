#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab
"""
"""
class YYMR(meep_utils.AbstractMeepModel): #{{{
    """ 
            MMxxxxxxMM
            MMxxxxxxMM
     ^ z    MMxxxxxxMM
     |      MMxxxxxxMM
     |      MMxxxxxxMM
     +-----> y

    FD 2014-02
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.01),
                meep_materials.material_Metal_THz(where = self.where_metal)
                ]

        #if 'TiO2' in comment:
            #self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        #elif 'STO' in comment:
            #self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        #elif 'DielLossless' in comment:
            #self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.0)]
        #elif 'Diel' in comment:
            #self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.01)]
        #else:
            #self.materials = [meep_materials.material_Metal_THz(where = self.where_wire)]
#

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = resolution/1.8
        self.size_y = yspacing
        self.size_z = 200e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 5000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 2000e9) 

        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14), plot_conductivity=True)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_diel(self, r):
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=cellz, d=self.radius)) and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                return self.return_value
        return 0
    def where_metal(self, r):
        if where_diel(self, r):
            return 0
        for cellz in self.cell_centers():
            #if (in_zslab(r, cz=cellz, d=self.radius)) and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                #return self.return_value
        return 0
#}}}

