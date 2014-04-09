#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xslab, in_yslab, in_zslab
"""
"""


class PKCutSheet_model(meep_utils.AbstractMeepModel): 
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
         |PMLPMLPML|                 #}}}

     ,gPPRg,     
    dXP^^\XYb    
    8XL ,/X(8   
    YbAAAMMMP    
     "8ggg8"     

    The purpose of computing is insight, not numbers. -- Richard W. Hamming

    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=200e-12, resolution=6e-6, hole=0e-6, cells=1, monzc=0e-6, Kx=0, Ky=0,
            cutw=20e-6, cutl=80e-6, thick=50e-6, holecy=0e-6, period=100e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "PKCutSheet"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2), 
                ]

        ## Dimension constants for the simulation
        monzd = thick+100e-6            ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = period, period, 200e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 50e-6
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
            xb = in_xslab(r, cx=0, d=self.size_x-self.cutw) ## define the junctions
            yb = in_yslab(r, cy=0, d=self.size_y-self.cutw) 
            xs = in_xslab(r, cx=0, d=self.size_x-self.cutl) ## define the central square
            ys = in_yslab(r, cy=0, d=self.size_y-self.cutl)
            xh = in_xslab(r, cx=0, d=self.hole)         ## define square hole in the middle
            yh = in_yslab(r, cy=self.holecy, d=self.hole)
            if (in_zslab(r, cz=cellz, d=self.thick) and (xs or ys or (xb and yb))) and not (xh and yh): 
                return self.return_value
        return 0

#}}}


