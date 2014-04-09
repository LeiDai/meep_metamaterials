#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab
"""
"""

class dielbar_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of dielectric bars or blocks #{{{
    FD 2013-07-13

         |PMLPMLPML|
    ^    +---------+ <-- monitor plane 2
    |    |         |
    z    |         | <-- Bloch-periodic boundaries on X, Y faces
         |  +---+  |    
         |  |XXX|  |<------- bar of TiO2
         |  +---+  |        
         |         |
         |         |
    y    +---------+ <-- monitor plane 1 
    -->  +=========+ <-- source
         |PMLPMLPML|                 #}}}
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, cells=1, monzc=0e-6, Kx=0, Ky=0,
            xs=55e-6, ys=95e-6, zs=25e-6, xyspacing=150e-6, monzd=100e-6, padding=0e-6, metalspacing=0e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "DielectricBar"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        #self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)]
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=4.)]
        self.materials = [
                meep_materials.material_dielectric(where = self.where_TiO2, eps=12.), 
                #meep_materials.material_STO(where = self.where_TiO2), 
                #meep_materials.material_STO_THz(where = self.where_TiO2), 
                #meep_materials.material_STO_hiloss(where = self.where_TiO2), 
                #meep_materials.material_Metal_THz(where = self.where_metal)
                ]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.size_x = resolution/1.8 if (xs==np.inf) else xyspacing
        self.size_y = resolution/.9 if (ys==np.inf) else xyspacing
        self.size_z = 300e-6+cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 30e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc, (monzd*cells)/2+monzc)  
        self.simtime = simtime      # [s]
        self.srcWidth = 2000e9     
        self.srcFreq = 1000e9     
        #self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 3000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback


        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14))

    def where_TiO2(self, r):
        for cellz in self.cell_centers():
            zz = in_zslab(r, cz=cellz,    d=self.zs)  
            yy = in_yslab(r, cy=0,          d=self.ys)
            xx = in_xslab(r, cx=0,          d=self.xs)
            if (zz and yy and xx): 
                return self.return_value
        return 0
    def where_metal(self, r): ## XXX
        for cellz in self.cell_centers():
            zz = in_zslab(r, cz=cellz+self.zs/2+10e-6/2 + self.metalspacing, d=10e-6)  
            xx = in_xslab(r, cx=0,    d=self.xyspacing-10e-6)  
            yy = in_yslab(r, cy=0,          d=self.xyspacing-10e-6)
            if (zz and yy and xx): 
                return self.return_value
        return 0
#}}}

class XCylWire_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            radius=10e-6, yspacing=100e-6, monzd=100e-6, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters



        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'DielLossless' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.0)]
        elif 'DielLoLoss' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.005)]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.05)]
        else:
            self.materials = [meep_materials.material_Metal_THz(where = self.where_wire)]


        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = resolution/1.8
        self.size_x = yspacing #XXX
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

    def where_wire(self, r):
        for cellz in self.cell_centers():
            if (in_xcyl(r, cy=self.resolution/4, cz=cellz, rad=self.radius)): 
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

class YCylWire_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            radius=20e-6, xspacing=100e-6, monzd=100e-6, epsilon=600):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "YCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.001)]
        else:
            self.materials = [meep_materials.material_Metal_THz(where = self.where_wire)]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = xspacing
        self.size_y = resolution/1.8
        self.size_z = 200e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        self.monitor_z1 =-(monzd*cells)/2+self.monzc - self.padding 
        self.monitor_z2 = (monzd*cells)/2+self.monzc + self.padding  
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        self.srcFreq = 1e12 #4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 3000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_wire(self, r):
        for cellz in self.cell_centers():
            if (in_ycyl(r, cx=self.resolution/4, cz=cellz, rad=self.radius)): 
                return self.return_value
            #if (in_zslab(r, cz=cellz, d=self.radius) and in_yslab(r, cy=.7e-6, d=60e-6-self.radius)):  # XXX
                #return self.return_value
        return 0
#}}}
class XRectWire_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="STO", simtime=100e-12, resolution=2e-6, cells=1, Kx=0, Ky=0,
            width=64e-6, thick=26e-6, yspacing=96e-6, monzd=96e-6, epsilon=600, padding=50e-6):

        ## XXX
        padding=100e-6
        ## XXX

        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XRectWire"    ## 
        self.monzc=0e-6 
        self.register_locals(locals())          ## Remember the parameters

        self.padding=50e-6
        self.monzc=0e-6 
        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'DielLossless' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0)]
        elif 'DielLoloss' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.01)]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=.1)]
        else:
            self.materials = [meep_materials.material_Metal_THz(where = self.where_wire)]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = resolution/1.8
        self.size_y = yspacing
        self.size_z = 200e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        self.monitor_z1 =-(monzd*cells)/2+self.monzc - self.padding 
        self.monitor_z2 = (monzd*cells)/2+self.monzc + self.padding  
                
                 
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        #self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.srcFreq = 1e12
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_wire(self, r):
        for cellz in self.cell_centers():
            #if (in_xcyl(r, cy=self.resolution/4, cz=cellz+self.resolution/4, rad=self.radius)): 
                #return self.return_value
            if (in_zslab(r, cz=cellz, d=self.thick) and in_yslab(r, cy=.7e-6, d=self.width)):  # XXX
                return self.return_value
        return 0
#}}}
class XRectWireMet_model(meep_utils.AbstractMeepModel): #{{{
    """  
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cells=1, Kx=0, Ky=0,
            width=64e-6, thick=26e-6, yspacing=96e-6, monzd=96e-6, epsilon=100, padding=50e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XRectWire"    ## 
        self.monzc=0e-6 
        self.register_locals(locals())          ## Remember the parameters

        self.padding=50e-6
        self.monzc=0e-6 
        ## Initialization of materials used
        self.materials = [ meep_materials.material_dielectric(where = self.where_slab, eps=epsilon, loss=0.01),
                meep_materials.material_Metal_THz(where = self.where_met)]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = resolution/1.8
        self.size_y = yspacing
        self.size_z = 120e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        self.monitor_z1 =-(monzd*cells)/2+self.monzc - self.padding 
        self.monitor_z2 = (monzd*cells)/2+self.monzc + self.padding  
                
                 
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        #self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.srcFreq = 1e12
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_slab(self, r):
        for cellz in self.cell_centers():
            if (in_zslab(r, cz=cellz, d=self.thick) and in_yslab(r, cy=0e-6, d=self.width)):  # XXX
                return self.return_value
        return 0
    def where_met(self, r):
        for cellz in self.cell_centers():
            #if (in_xcyl(r, cy=self.resolution/4, cz=cellz+self.resolution/4, rad=self.radius)): 
                #return self.return_value

            ## Metallisation along Z
            #if ((in_zslab(r, cz=cellz, d=self.thick+4*self.resolution) and not in_zslab(r, cz=cellz, d=self.thick)) \
                    #and in_yslab(r, cy=.7e-6, d=2*self.width)):  # XXX
                #return self.return_value
            if (in_zslab(r, cz=cellz, d=self.thick)) \
                    and (in_yslab(r, cy=0e-6, d=self.width+max(3e-6,2*self.resolution)) and not in_yslab(r, cy=0e-6, d=self.width)):  # XXX
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
        self.monitor_z1 =-(monzd*cells)/2+self.monzc - self.padding 
        self.monitor_z2 = (monzd*cells)/2+self.monzc + self.padding  
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


class Fishnet_model(meep_utils.AbstractMeepModel): #{{{
    """  
    FD 2013-09-20
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=7e-6, cells=1, monzc=0e-6, Kx=0, Ky=0,
            xs=50e-6, ys=50e-6, zs=12e-6, xyspacing=100e-6, monzd=100e-6, padding=50e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "PKCTest"    ## 
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_Metal_THz(where = self.where_metal), 
                ]

        ## Dimension constants for the simulation
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 200e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 30e-6
        self.monitor_z1 =-(monzd*cells)/2+self.monzc - self.padding 
        self.monitor_z2 = (monzd*cells)/2+self.monzc + self.padding  
        self.simtime = simtime      # [s]
        self.srcWidth = 4000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 5000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback


    #def where_metal(self, r):
        #for cellz in self.cell_centers():
            #if (in_xslab(r, cx=0e-6, d=self.wtth) or in_yslab(r, cy=0e-6, d=self.wtth)) and \
                    #in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                #return self.return_value
        #return 0

    def where_metal(self, r):
        for cellz in self.cell_centers():
            zz = in_zslab(r, cz=cellz,    d=self.zs)  
            yy = not in_yslab(r, cy=self.resolution/2,   d=self.ys)
            xx = not in_xslab(r, cx=self.resolution/2,   d=self.xs)
            if (zz and (yy or xx)): 
                return self.return_value
        return 0

#}}}

class XCylWire_model_test(meep_utils.AbstractMeepModel): #{{{
    """  Array of metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=100e-6,
            radius=4e-6, spacing=100e-6, monzd=100e-6, epsilon=600, xlen=50e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "XCylWire"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        self.materials = [
                meep_materials.material_Metal_THz(where = self.where_particle),
                #meep_materials.material_TiO2_THz(where = self.where_particle),
                meep_materials.material_dielectric(where = self.where_filling, eps=epsilon)]

        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 100e-6
        self.size_x = spacing
        self.size_y = spacing
        self.size_z = 400e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding , (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 1000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_particle(self, r):
        for cellz in self.cell_centers():
            if (in_xcyl(r, cy=self.resolution/4, cz=cellz+self.resolution/4, rad=self.radius) and 
                in_xslab(r, cx=self.resolution/4, d=self.xlen)) : 
                return self.return_value
        return 0
    def where_filling(self, r):
        for cellz in self.cell_centers():
            if (in_xcyl(r, cy=self.resolution/4, cz=cellz+self.resolution/4, rad=self.radius) and 
                not in_xslab(r, cx=self.resolution/4, d=self.xlen)) : 
                return self.return_value
        return 0
#}}}

class Wedge_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="TiO2", simtime=10e-12, resolution=3e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            radius=10e-6, yspacing=100e-6, zspacing=100e-6, monzd=200e-6, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "Wedge"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'DielLossless' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.0)]
        elif 'DielLoLoss' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.005)]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.05)]
        else:
            self.materials = [meep_materials.material_Metal_THz(where = self.where_wire)]


        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 20e-6
        self.size_x = resolution/1.8
        self.size_y = 2200e-6
        self.size_z = 2000e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 2000e9) 

        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14), plot_conductivity=True)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_wire(self, r):
        y,z = r.y(), r.z()
        if z < self.size_z*(-.3) or z+y/2>0: 
            return 0
        yy,zz = np.arccos(np.cos(y/self.yspacing*np.pi*2))*self.yspacing/(np.pi*2), np.arccos(np.cos(z/self.zspacing*np.pi*2))*self.zspacing/(np.pi*2)
        if (yy**2+zz**2)**.5 < self.radius:
            return self.return_value
        return 0
#}}}
