#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xslab, in_yslab, in_zslab
"""
"""


class CKEBars_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of dielectric bars along the E-field, proposed by Christelle (June 2013)  #{{{
    FD 2013-06-03

         |PMLPMLPML|
    ^    +---------+ <-- monitor plane 2
    |    |         |
    z    |         | <-- Bloch-periodic boundaries on X, Y faces
         |  +---+  |    
         |  |XXX|  |<------- bar of TiO2
         |--+---+--|<-- mylar substrate 
         |---------|
    |    |         |
    x    +---------+ <-- monitor plane 1 
    -->  +=========+ <-- source
         |PMLPMLPML|                 #}}}
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=1000e-12, resolution=5e-6, cells=1, monzc=0e-6,
            BarWidth=100e-6, BarThick=100e-6, SubstrThick=50e-6, BarPeriod=200e-6, YCellShift=0,
            Kx=0, Ky=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "CKEBars"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Initialization of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2), 
                meep_materials.material_dielectric(where = self.where_substr, eps=1.0)] ## XXX

        ## Dimension constants for the simulation
        monzd = BarThick + SubstrThick            ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = resolution/1.8, BarPeriod, 500e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 100e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc, (monzd*cells)/2+monzc)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 500e9, 1000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    #def where_metal(self, r):
        #for cellz in self.cell_centers():
            #if (in_xslab(r, cx=0e-6, d=self.wtth) or in_yslab(r, cy=0e-6, d=self.wtth)) and \
                    #in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                #return self.return_value
        #return 0

    def where_TiO2(self, r):
        # callback used for each polarizability of each material (materials should never overlap)
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=cellz+self.SubstrThick/2,d=self.BarThick) and
                    in_yslab(r, cy=0,d=self.BarWidth)): 
                return self.return_value
        return 0

    def where_substr(self, r):
        for cellz in self.cell_centers():
            if in_zslab(r, cz=cellz-self.BarThick/2,d=self.SubstrThick): 
                return self.return_value
        return 0
#}}}

class CKEBars_model_test(meep_utils.AbstractMeepModel): #{{{
    """  Array of dielectric bars along the E-field, proposed by Christelle (June 2013)  #{{{
    Testing version with additional features
    FD 2013-06-03

         |PMLPMLPML|
    ^    +---------+ <-- monitor plane 2
    |    |         |
    z    |         | <-- Bloch-periodic boundaries on X, Y faces
         |  +---+  |    
         |  |XXX|  |<------- bar of TiO2
         |  +---+  |        <-- (no substrate in this case)
         |         |
         |         |
    y    +---------+ <-- monitor plane 1 
    -->  +=========+ <-- source
         |PMLPMLPML|                 #}}}
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
        self.simulation_name = "CKEBars"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Initialization of materials used
        #self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)]
        self.materials = [meep_materials.material_Sapphire_THz(where = self.where_TiO2), 
                meep_materials.material_dielectric(where = self.where_substr, eps=4.0)]

        ## Dimension constants for the simulation
        monzd = BarPeriod          ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = resolution, BarPeriod, 500e-6 + cells*self.monzd
        #self.size_x, self.size_y, self.size_z = BarPeriod, BarPeriod, 500e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 100e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc, (monzd*cells)/2+monzc)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 500e9, 1000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_TiO2(self, r):
        for cellz in self.cell_centers():
            #zz = in_zslab(r, cz=cellz + (self.BarPeriod-self.BarThick)/2, d=self.BarThick)  
            zz = in_zslab(r, cz=cellz, d=self.BarThick)  
            yy = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift, d=self.BarWidth)
            yy2 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift + self.size_y, d=self.BarWidth)
            yy3 = in_yslab(r, cy=(cellz/self.monzd)*self.YCellShift - self.size_y, d=self.BarWidth)
            #xx = True
            #xx = in_xslab(r, cx=self.resolution/4, d=self.size_x - self.XCut)
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
