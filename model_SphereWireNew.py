#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
"""
"""

class SphereWire_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of spheres with a metallic mesh around to provide possibly negative epsilon.   #{{{
    FD 2013-01-10

         |PMLPMLPML|
    ^    +---------+ <-- monitor plane 2
    |    |         |
    z    |_       _| <-- Bloch-periodic boundaries on X, Y faces      
         | \     /.|    
   wire---->| o |  |
   mesh  |_/     \_|<-- dielectric TiO2 sphere
         |         |
    x    +---------+ <-- monitor plane 1 
    -->  +=========+ <-- source
         |PMLPMLPML|

    An X-Y view:
         +--------+--------+
         |__/     I     \__| 	<--- sphere
         |        I        |
         |        I        |
         |        I        |
         |========+========| 	<---- wire cross, possibly with crossing-split and near-sphere-split
         |        I        |
         |        I        |
         |__      I      __|
         |  \     I     /  |
         +--------+--------+
                  <-> wireconnw = WireConnectorWidth (increased capacitive coupling between wires)

    By default, using the geometry from Yakiyama2012
    Note: The plasma frequency of the wire array may be calculated as:
        a, r = 1e-4, 10e-6          ## wire spacing and diameter [Pendry1996]
        f = (6.28*3e8**2/a**2/numpy.log(a/r))**.5  / 6.28 / 1e12
        >>> 789e9 ## Hz
    """
#}}}
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, cells=1, monzc=0e-6, padding=50e-6,
            radius=25e-6, wzofs=0e-6, spacing=75e-6, wlth=10e-6, wtth=10e-6, Kx=0, Ky=0, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SphereWire"    ## 
        monzd=spacing
        
        
        self.register_locals(locals())          ## Remember the parameters
        #wlth=1e-6
        #wtth=90e-6 
        #wlth=6e-6
        #wtth=1000e-6 
        

        ## Definition of materials used
        #self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=4.)]
        #self.materials = [meep_materials.material_TiO2_THz_HIEPS(where = self.where_TiO2), 
                #meep_materials.material_Metal_THz(where = self.where_metal)]
        ## Initialization of materials used
        if 'TiO2NC' in comment:
            self.materials = [meep_materials.material_TiO2_NC_THz(where = self.where_TiO2), ]
                    #meep_materials.material_Metal_THz(where = self.where_metal) ]
        elif 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2), 
                    meep_materials.material_Metal_THz(where = self.where_metal) ]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=epsilon),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]
        else:
            self.materials = [meep_materials.material_STO_THz(where = self.where_TiO2),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]


        ## constants for the simulation
        self.pml_thickness = 20e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 2000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0e9, 1000e9) 

        ## Dimension constants for the simulation
        self.size_x = spacing 
        self.size_y = spacing 
        self.size_z = 60e-6 + cells*monzd + 2*self.pml_thickness + 2*self.padding
        #self.size_x, self.size_y, self.size_z = resolution/1.8, spacing, size_z ## two dimensional case XXX

        self.TestMaterials()
        print "CellCenters:", self.cell_centers()

    def where_diel(self, r):
        for cellz in self.cell_centers():
            if in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                return self.return_value
        return 0

    def where_metal(self, r):
        return 0              ## nometal XXX

        for cellz in self.cell_centers():
            if in_xcyl(r, cy=0e-6, cz=self.wzofs+cellz, rad=self.wlth) or in_ycyl(r, cx=0e-6, cz=self.wzofs+cellz, rad=self.wlth):
                return self.return_value
        return 0              ## nometal
        for cellz in self.cell_centers():

            #if (in_xslab(r, cx=0e-6, d=self.wtth) or in_yslab(r, cy=0e-6, d=self.wtth)) and \
            if (in_yslab(r, cy=0e-6, d=self.wtth)) and \
                    in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                return self.return_value
        return 0

    def where_TiO2(self, r):
        # callback used for each polarizability of each material (materials should never overlap)
        for cellz in self.cell_centers():
            ## XXX ## avoid the metal volume
            #if (in_xslab(r, cx=0e-6, d=self.wtth) or in_yslab(r, cy=0e-6, d=self.wtth)) and \
                    #in_zslab(r, cz=self.wzofs+cellz, d=self.wlth):
                #return 0                

            # Y-bars
            #if   in_ycyl(r, cx=-self.size_x/2, cz=cellz, rad=self.radius) \
              #or in_ycyl(r, cx= self.size_x/2, cz=cellz, rad=self.radius):
                #return self.return_value            # (do not change this return value)
            # X-bars
            #if   in_xcyl(r, cy=-self.size_y/2, cz=cellz, rad=self.radius) \
              #or in_xcyl(r, cy= self.size_y/2, cz=cellz, rad=self.radius):
                #return self.return_value            # (do not change this return value)
            ## /XXX 

            if   in_sphere(r, cx=0, cy=0, cz=cellz, rad=self.radius):
               return self.return_value            # (do not change this return value) XXX

            #if   in_sphere(r, cx=-self.size_x/2, cy=-self.size_y/2, cz=cellz, rad=self.radius) \
              #or in_sphere(r, cx= self.size_x/2, cy=-self.size_y/2, cz=cellz, rad=self.radius) \
              #or in_sphere(r, cx=-self.size_x/2, cy= self.size_y/2, cz=cellz, rad=self.radius) \
              #or in_sphere(r, cx= self.size_x/2, cy= self.size_y/2, cz=cellz, rad=self.radius):
               #return self.return_value            # (do not change this return value)
        return 0
#}}}


class SphereFishnet_model(meep_utils.AbstractMeepModel): #{{{
    """  
    Array of diel spheres embedded in metallized mylar foil with holes
    Enables electric field tuning of resonances?
    FD 2013-10-20
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=100e-12, resolution=5e-6, cells=1, monzc=0e-6, padding=50e-6,
            radius=25e-6, wzofs=0e-6, spacing=75e-6, thick=20e-6, Kx=0, Ky=0, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SphereWire"    ## 
        monzd=spacing
        
        
        self.register_locals(locals())          ## Remember the parameters
        #wlth=1e-6
        #wtth=90e-6 
        #wlth=6e-6
        #wtth=1000e-6 
        

        ## Definition of materials used
        #self.materials = [meep_materials.material_dielectric(where = self.where_diel, eps=4.)]
        #self.materials = [meep_materials.material_TiO2_THz_HIEPS(where = self.where_TiO2), 
                #meep_materials.material_Metal_THz(where = self.where_metal)]
        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2), 
                    meep_materials.material_dielectric(where = self.where_diel, eps=2.),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]
        elif 'Diel' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=epsilon),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]
        else:
            self.materials = [meep_materials.material_STO_THz(where = self.where_TiO2),
                    meep_materials.material_dielectric(where = self.where_diel, eps=2.),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]


        ## constants for the simulation
        self.pml_thickness = 20e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 2000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0e9, 1000e9) 

        ## Dimension constants for the simulation
        self.size_x = spacing 
        self.size_y = spacing 
        self.size_z = 60e-6 + cells*monzd + 2*self.pml_thickness + 2*self.padding
        #self.size_x, self.size_y, self.size_z = resolution/1.8, spacing, size_z ## two dimensional case XXX


        self.r2 = radius + 6e-6

        self.TestMaterials()
        print "CellCenters:", self.cell_centers()

    def where_diel(self, r):
        for cellz in self.cell_centers():
            if not (in_zcyl(r, cx=-self.size_x/2, cy=-self.size_y/2, rad=self.r2) \
              or in_zcyl(r, cx= self.size_x/2, cy=-self.size_y/2,    rad=self.r2) \
              or in_zcyl(r, cx=-self.size_x/2, cy= self.size_y/2,    rad=self.r2) \
              or in_zcyl(r, cx= self.size_x/2, cy= self.size_y/2,    rad=self.r2)) \
              and in_zslab(r, cz=cellz, d=self.thick):
                    return self.return_value
        return 0

    def where_metal(self, r):
        for cellz in self.cell_centers():
            if not (in_zcyl(r, cx=-self.size_x/2, cy=-self.size_y/2, rad=self.r2) \
              or in_zcyl(r, cx= self.size_x/2, cy=-self.size_y/2,    rad=self.r2) \
              or in_zcyl(r, cx=-self.size_x/2, cy= self.size_y/2,    rad=self.r2) \
              or in_zcyl(r, cx= self.size_x/2, cy= self.size_y/2,    rad=self.r2)) \
              and (in_zslab(r, cz=cellz, d=self.thick+self.resolution*4) and not in_zslab(r, cz=cellz, d=self.thick)):
                return self.return_value
        return 0           

    def where_TiO2(self, r):
        for cellz in self.cell_centers():
            if   in_sphere(r, cx=-self.size_x/2, cy=-self.size_y/2, cz=cellz+self.wzofs, rad=self.radius) \
              or in_sphere(r, cx= self.size_x/2, cy=-self.size_y/2, cz=cellz+self.wzofs, rad=self.radius) \
              or in_sphere(r, cx=-self.size_x/2, cy= self.size_y/2, cz=cellz+self.wzofs, rad=self.radius) \
              or in_sphere(r, cx= self.size_x/2, cy= self.size_y/2, cz=cellz+self.wzofs, rad=self.radius):
               return self.return_value
        return 0
#}}}

class SimpleEllipsoid_model(meep_utils.AbstractMeepModel): #{{{
    """  
    FD 2014-01-28
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cells=1, monzc=0e-6, spacing=75e-6,
            xradius=15e-6, yradius=15e-6, zradius=15e-6, padding=0, Kx=0, Ky=0, metalpos=1000):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SE"    
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)]

        ## Dimension constants for the simulation
        monzd = spacing            ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = spacing, spacing, 60e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 10e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 7000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0., 5000e9) 

        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_TiO2(self, r):
        # callback used for each polarizability of each material (materials should never overlap)
        for cellz in self.cell_centers():
            #xd, yd, zd = (cx-r.x()), (cy-r.y()), (cz-r.z())
            xd, yd, zd = r.x(), r.y(), r.z()
            if ((xd)**2/self.xradius**2 + (yd)**2/self.yradius**2 + (zd)**2/self.zradius**2)**.5 < 1 :
                return self.return_value
        return 0
#}}}

class SphereElliptic_model(meep_utils.AbstractMeepModel): #{{{
    """  
    FD 2014-01-28
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=10e-12, resolution=2e-6, cells=1, monzc=0e-6, spacing=120e-6,
            radius1=15e-6, radius2=12e-6, padding=0, Kx=0, Ky=0, metalpos=1000):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SE"    
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2),
                    meep_materials.material_Metal_THz(where = self.where_metal) ]

        ## Dimension constants for the simulation
        monzd = spacing            ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = spacing, spacing, 60e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 10e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 7000e9     # [Hz], note: "Last source time" = 10/srcWidth
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
            if r.z() > self.metalpos:
                return self.return_value
        return 0

    def where_TiO2(self, r):
        # callback used for each polarizability of each material (materials should never overlap)
        for cellz in self.cell_centers():
            #xd, yd, zd = (cx-r.x()), (cy-r.y()), (cz-r.z())
            xd, yd, zd = r.y()*.866 - r.x()*.5, r.y()*.5 + r.x()*.866, r.z()
            if ((xd+yd)**2/2./self.radius1**2 + (xd-yd)**2/2./self.radius2**2 + zd**2/(self.radius1*self.radius2))**.5 < 1 :
                return self.return_value
        return 0

    #def where_substr(self, r):
        #for cellz in self.cell_centers():
            #if in_zslab(r, cz=cellz-self.BarThick/2,d=self.SubstrThick): 
                #return self.return_value
        #return 0
#}}}


class SphereMitrof_model(meep_utils.AbstractMeepModel): #{{{
    """  
    FD 2014-01-28
    """
    def cell_centers(self):#{{{
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: adds surrounding two cells!!
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)#}}}
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cells=1, monzc=0e-6, spacing=30e-6,
            radius=15e-6, ex=2, padding=0, Kx=0, Ky=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SE"    
        self.register_locals(locals())          ## Remember the parameters


        ## Definition of materials used
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)] 

        ## Dimension constants for the simulation
        monzd = spacing            ## monitor z-distance
        self.monzd = monzd
        self.size_x, self.size_y, self.size_z = spacing, spacing, 60e-6 + cells*self.monzd

        ## constants for the simulation
        self.pml_thickness = 10e-6
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 4000e9     # [Hz], note: "Last source time" = 10/srcWidth
        self.interesting_frequencies = (0., 5000e9) 

        print self.ex
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
            if in_ellipsoid(r, cx=0, cy=0, cz=0, rad=self.radius, ex=self.ex):
                return self.return_value
        return 0

    #def where_substr(self, r):
        #for cellz in self.cell_centers():
            #if in_zslab(r, cz=cellz-self.BarThick/2,d=self.SubstrThick): 
                #return self.return_value
        #return 0
#}}}
