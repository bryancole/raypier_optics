from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Bool

from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor

from raypier.core.cmaterials import OpaqueMaterial

from raypier.bases import Traceable
from raypier.prisms import Extrusion
from raypier.splines import Extruded_interpolant
from raypier.core.cfaces import PolygonFace, ExtrudedPlanarFace

from raypier.core.ctracer import FaceList
import numpy as np 

class LDLF(Extrusion):
    name = "Linear Dialectric Light Funnle"
    slat_width= Float #width of slats
    ap_width = Float #width of exit apperture
    slant = Float #angle of slanted sides
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       Item('n_inside'),
                       Item('length'),
                       Item('slat_width'),
                       Item('ap_width'),
                       Item('slant')
                       )
                       )

    @on_trait_change("slat_width, ap_width, slant")
    def config_profile(self):
        theta = self.slant*np.pi/180.
        l = self.slat_width
        h = l*np.sin(theta)
        dis = l*np.cos(theta)
        w = self.ap_width/2
        self.profile = [(-w,0),
                        (w,0),
                        (w+dis,h),
                        (-w-dis,h)]

class CylindricalLDLF(Extruded_interpolant):
    #subclassing the Bezier extrusion to build a multipart optical element.
    name = "Cylindrical Linear Dialectric Light Funnle"
    entry_width = Float
    radius = Float
    slant = Float
    depth = Float
    aperture_material = OpaqueMaterial()


    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl

    def make_faces(self):
        radius = self.radius
        opening = .5*self.entry_width
        d=self.depth
        slant = self.slant
        z1 = self.z_height_1
        z2 = self.z_height_2
        m=self.material

        if radius**2 < opening**2:
            print("error: radius is too small (<.5 width)")
            return

        if d/np.tan(np.radians(slant)) > opening:
            print("The Cyclindrical LDLF is too deep, slanted sides intersect.")
            return

        sign = 1
        if radius < 0:
            sign = -1
        
        #lens_profile
        x = np.linspace(-opening,opening)
        displace = np.sqrt(radius**2-opening**2)
        ys = sign*(np.sqrt(radius**2-x**2)-displace)
        self.profile = [x,ys]
        self.trace_ends=True
        self.trace_top = False
        sides = Extruded_interpolant.make_faces(self)    #list of faces

        #one flat, slanty side
        x1 = opening
        y1 = 0
        x2 = opening-(d/np.tan(np.radians(slant)))
        y2 = -d

        sides.append(ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=x1, y1=y1, 
                    x2=x2, y2=y2, material=m))

        #the other, symmetrical one
        sides.append(ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=-x1, y1=y1, 
                    x2=-x2, y2=y2, material=m,invert_normal=True))

        #the target
        sides.append(ExtrudedPlanarFace(owner=self, z1=z1, z2=z2, x1=x2, y1=y2, 
                    x2=-x2, y2=y2, material=self.aperture_material))
        
        
        if self.trace_ends:
            #warning, just polygon approximations.  not water tight.
            profile = self.get_real_profile()
            profile[0] = np.append(profile[0], [x1,x2,-x2,-x1],1)
            profile[1] = np.append(profile[1], [y1,y2,y2,y1],1)
            base = PolygonFace(owner=self, z_plane=z1,
                        xy_points=profile, material=m)
            top = PolygonFace(owner=self, z_plane=z2, material=m,
                        xy_points=profile, invert_normal=True)
            sides.extend([base, top])

        return sides

    def _pipeline_default(self):
        #pretty much just copied this from extruded_interpolant, and added three profile points
        source = self.data_source
        def execute():

            opening = .5*self.entry_width
            d=self.depth
            slant = self.slant

            x1 = opening
            y1 = 0
            x2 = opening-(d/np.tan(np.radians(slant)))
            y2 = -d

            profile = self.get_real_profile()
            profile[0] = np.append(profile[0], [x2,-x2],1)
            profile[1] = np.append(profile[1], [y2,y2],1)
            xy = np.column_stack((profile[:][0],profile[:][1]))
            z = np.ones(xy.shape[0]) * self.z_height_1
            points = np.column_stack((xy,z))
            
            cells = [list(range(len(z))),]
            
            output = source.poly_data_output
            output.points = points
            output.polys = cells
        source.set_execute_method(execute)
        
        self.extrude.scale_factor = self.z_height_2 - self.z_height_1  #mm, put here because it wasn;t being initialized
        if self.trace_ends:
            print("drew ends")
            self.extrude.capping = True 
        extrude = self.extrude
        extrude.input = source.output
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=extrude.output_port, 
                                      transform=t)
        return transf


