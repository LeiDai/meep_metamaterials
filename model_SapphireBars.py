#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xslab, in_yslab, in_zslab
import  meep_utils
"""
"""

class SapphireBars(meep_utils.AbstractMeepModel): #{{{
    """  Array of sapphire bars along the E-field
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=200e-12, resolution=5e-6, cells=1, monzc=0e-6,
            BarWidth=50e-6, BarThick=50e-6, BarPeriod=100e-6, YCellShift=0, XCut=0,
            Kx=0, Ky=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SapphireBars"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Initialization of materials used
        #self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)]
        self.materials = [
                meep_materials.material_Sapphire_THz(where = self.where_Sapphire), 
                meep_materials.material_Sapphire_THz(where = self.where_Sapphire, ordinary=0)]
                #meep_materials.material_TiO2_THz(where = None),]
                #meep_materials.material_dielectric(where = self.where_substr, eps=4.0),
                #meep_materials.material_Metal_THz(where = self.where_metal)]

        meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14))

        ## Dimension constants for the simulation
        monzd = BarPeriod          ## monitor z-distance
        self.monzd = monzd
        #self.size_x, self.size_y, self.size_z = resolution, BarPeriod, 500e-6 + cells*self.monzd
        self.size_x, self.size_y, self.size_z = BarPeriod, BarPeriod, 500e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 30e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc, (monzd*cells)/2+monzc)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 500e9, 1000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0., 1800e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_metal(self, r):
        for cellz in self.cell_centers():
            #zz = in_zslab(r, cz=cellz + (self.BarPeriod-self.BarThick)/2, d=self.BarThick)  
            zz = in_zslab(r, cz=cellz, d=self.BarThick+10e-6) and not in_zslab(r, cz=cellz, d=self.BarThick) 
            yy = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift, d=self.BarWidth)
            yy2 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift + self.size_y, d=self.BarWidth)
            yy3 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift - self.size_y, d=self.BarWidth)
            xx = in_xslab(r, cx=self.resolution/4, d=self.size_x - 10e-6)
            if (zz and (yy or yy2 or yy3) and xx): 
                return self.return_value
        return 0
    def where_Sapphire(self, r):
        for cellz in self.cell_centers():
            #zz = in_zslab(r, cz=cellz + (self.BarPeriod-self.BarThick)/2, d=self.BarThick)  
            zz = in_zslab(r, cz=cellz, d=self.BarThick)  
            yy = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift, d=self.BarWidth)
            yy2 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift + self.size_y, d=self.BarWidth)
            yy3 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift - self.size_y, d=self.BarWidth)
            #xx = True
            if (zz and (yy or yy2 or yy3)): 
                return self.return_value
        return 0
    def where_substr(self, r):
        return 0
        for cellz in self.cell_centers():
            if in_zslab(r, cz=cellz-self.BarThick/2, d=self.monzd-self.BarThick): 
                return self.return_value
        return 0
#}}}
