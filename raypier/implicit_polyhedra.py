
from traits.api import Instance, List, Float, on_trait_change
from traitsui.api import View, Item, VGroup, InstanceEditor
from tvtk.api import tvtk

from raypier.core.ctracer import FaceList
from raypier.core.cfaces import ImplicitBoundedPlanarFace
from raypier.core.cimplicit_surfs import Intersection, Union, Difference
from raypier.implicit_surfaces import Plane
from raypier.bases import Optic, Traceable
from raypier.materials import BaseOpticalMaterial, OpticalMaterial
from raypier.core.cmaterials import FullDielectricDispersiveMaterial

from itertools import combinations
import numpy as np
from numpy.linalg import solve, det
from collections import defaultdict
from scipy.spatial import Delaunay
from scipy.spatial.qhull import QhullError


def get_intersection_of_planes(*planes : "tuple(p1,p2,p3)") -> "None | tuple(x,y,z)":
    A = np.array([list(p.normal) for p in planes])
    b = np.array([np.dot(p.origin,p.normal) for p in planes])
    if abs(det(A)) < 1e-5:
        return None 
    return solve(A,b)


class Polyhedron(Optic):
    """
    Defines the 3D volumn enclosed by N planes.
    The plane normals are the *outward* normals.
    """
    name = "Polyhedron"
    planes = List(Plane)
    _tolerance = Float(1e-12)
    
    data_source = Instance(tvtk.ProgrammableSource, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('n_inside'),
                       )
                       )
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = self.make_faces()
        return fl
    
    def make_faces(self):
        planes = self.planes
        faces = []
        for plane in planes:
            bounding_c_planes = [p.c_surface for p in planes if p is not plane]
            face = ImplicitBoundedPlanarFace(origin=plane.origin,
                                             normal=plane.normal,
                                             boundary=Intersection(*bounding_c_planes)
                                             )
            faces.append(face)
            
        for f in faces:
            f.material=self.material
        return faces
    
    def get_vertices(self):
        """Find the vertices from the list of planes"""
        planes = self.planes
        tol = self._tolerance
        triples = combinations(planes, 3)
        imp_vol = Intersection(*[p.c_surface for p in planes])
        for triple in triples:
            pt : "tuple[float,float,float]" = get_intersection_of_planes(*triple)
            if pt is None:
                continue
            dist = imp_vol.evaluate(*pt)
            if dist <= tol: #negative is inside
                yield pt, triple
                
    def get_points_and_cells(self):
        dd = defaultdict(list)
        verts =[]
        cells=[]
        for pt_id, (pt, triple) in enumerate(self.get_vertices()):
            #print(pt_id, pt, [p.origin for p in triple])
            verts.append(pt)
            for p in triple:
                dd[p].append(pt_id)
        for p, p_verts in dd.items():
            n_verts = len(p_verts)
            if n_verts < 3: #Can't make a face with under 3 vertices
                continue
            if n_verts==3:
                # Don't both worry about the ordering. We'll use vtkPolyDataNormals to sort it out later...
                cells.append(p_verts)
            else:
                p_pts = np.array([verts[i] for i in p_verts])
                ### We use Priciple Component Analysis to transform into 2D U,V coords.
                centre = p_pts.mean(axis=0)
                cov = np.cov(p_pts-centre, rowvar=False)
                eig_vals, eig_vecs = np.linalg.eig(cov)
                order = np.argsort(eig_vals)
                ax1 = eig_vecs[:,order[-1]]
                ax2 = eig_vecs[:,order[-2]]
                u = np.dot(p_pts,ax1)
                v = np.dot(p_pts,ax2)
                pts_2d = np.column_stack([u,v])
                ### Then triangulate
                try:
                    tri = Delaunay(pts_2d)
                    for t in tri.simplices:
                        cells.append([p_verts[idx] for idx in t])
                except QhullError:
                    print(f"Failed to triangulate plane {p}\n"
                          f"with points: {p_pts}\n"
                          f"UV pts:\n{u}\n{v}\n"
                          f"Ax:\n{ax1}\n{ax2}\n"
                          f"eig_vals: {eig_vals}\n"
                          f"eig_vecs: {eig_vecs}")
                    
        return verts, cells
                
                
    def _pipeline_default(self):
        #self.config_profile()
        #self.config_pipeline()
        source = self.data_source
        def execute():
            print("RECAL VTK DATA")
            verts, cells = self.get_points_and_cells()
            print(verts)
            print(cells)
            output = source.poly_data_output
            output.points = verts #points, actually
            output.polys = cells
            
        source.set_execute_method(execute)
        print("Made PIPELINE")
        
        tr = tvtk.PolyDataNormals(input_connection=source.output_port)
        
        t = self.transform
        transf = tvtk.TransformFilter(input_connection=tr.output_port, 
                                      transform=t)
        return transf
    
    
class DispersivePolyhedron(Polyhedron):
    dispersion = Instance(BaseOpticalMaterial)
    
    @on_trait_change("n_inside, n_outside")
    def n_changed(self):
        pass
    
    def _dispersion_default(self):
        return OpticalMaterial(glass_name="N-BK7")
    
    def _material_default(self):
        dispersion_curve = self.dispersion.dispersion_curve
        m = FullDielectricDispersiveMaterial(dispersion_inside=dispersion_curve)
        return m
    

class RightAnglePrism(DispersivePolyhedron):
    """
    A right-angle prism which can include angle- and pyramid-errors on the side-faces.
    """
    ### The short side of the prism
    side_length = Float(10.0)
    ### height along the "extrusion axis"
    height = Float(10.0)
    
    angle_error = Float(0.0) #in arcmin
    pyramid_error = Float(0.0) #in arcmin
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('side_length'),
                       Item('height'),
                       Item('angle_error'),
                       Item('pyramid_error'),
                       Item('dispersion', style="custom", editor=InstanceEditor())
                       )
                       )
    
    
    
    def make_step_shape(self):
        from raypier.step_export import make_general_extrusion
        s_len = self.side_length
        h_len = np.sqrt(2*(s_len**2))/2
        y = self.height / 2.
        
        points = [(-h_len,-y,0.),(0.,-y,-h_len),(h_len,-y,0.)]
        length = self.height
        vector = (0.,1.,0.)
        
        shape = make_general_extrusion(self.centre, 
                                      self.direction, 
                                      self.x_axis, 
                                      points, length, vector)
        return shape, 'purple3'
    
    def _planes_default(self):
        return self.make_planes()
    
    @on_trait_change("side_length, height, angle_error, pyramid_error")
    def on_planes_update(self):
        self.data_source.modified()
        self.planes = self.make_planes()
        self.faces.faces = self.make_faces()
        self.update=True
    
    def make_planes(self):
        y = self.height / 2.
        side = self.side_length
        sa = np.sqrt((side**2)/2)
        root2 = np.sqrt(2)
        y_err_rads = np.pi*(self.pyramid_error/(60*180))
        y_err = root2*np.tan(y_err_rads)
        side_angle = np.pi*((45.0 + self.angle_error/60)/180)
        x = np.cos(side_angle)
        z = -np.sin(side_angle)
        p1 = Plane(origin=(0.,y,0.), normal=(0.,1.,0.))
        p2 = Plane(origin=(0.,-y,0.), normal=(0.,-1.,0.))
        p3 = Plane(origin=(0.,0.,0.), normal=(0.,0.,1.))
        p4 = Plane(origin=(0.,0.,-sa), normal=(x,y_err,z))
        p5 = Plane(origin=(0.,0.,-sa), normal=(-1.,0.,-1.))
        return [p1,p2,p3,p4,p5]
    
    
class Rhomboid(DispersivePolyhedron):
    name = "Rhomboid Polyhedra"
    abstract = False
    height = Float(15.0) #distance between parallel faces
    width = Float(15.0) #width of parallel faces
    slant = Float(45.0) #angle of the oblique faces
    z_height_1 = Float(0.0)
    z_height_2 = Float(20.0)
    
    angle_error = Float(0.0) #in arcmin
    pyramid_error = Float(0.0) #in arcmin
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('n_inside'),
                       Item('z_height_1'),
                       Item('z_height_2'),
                       Item('height'),
                       Item('width'),
                       Item('slant'),
                       Item('angle_error'),
                       Item('pyramid_error')
                       )
                       )
    
    @on_trait_change("height, width, slant, z_height_1, z_height_2, angle_error, pyramid_error")
    def on_planes_update(self):
        self.data_source.modified()
        self.planes = self.make_planes()
        self.faces.faces = self.make_faces()
        self.update=True
        
    def _planes_default(self):
        return self.make_planes()
    
    def make_planes(self):
        z_min, z_max = self.z_height_1,self.z_height_2
        if z_min > z_max:
            z_max, z_min = z_min, z_max
        p1 = Plane(origin=(0.,0.,z_min), normal=(0.,0.,-1.))
        p2 = Plane(origin=(0.,0.,z_max), normal=(0.,0.,1.))
        h = self.height/2
        w = self.width/2
        p3 = Plane(origin=(0.,-h,0.), normal=(0.,-1.,0.))
        p4 = Plane(origin=(0.,h,0.), normal=(0.,1.,0.))
        angle = self.slant*np.pi/180
        angle_err = (np.pi*self.angle_error/(2*60.*180))
        a = np.cos(angle + angle_err)
        b = np.sin(angle + angle_err)
        z_err = np.sin(np.pi*self.pyramid_error/(2*60.*180))
        p5 = Plane(origin=(-w,0.,0.), normal=(-a,b,z_err))
        a = np.cos(angle - angle_err)
        b = np.sin(angle - angle_err)
        p6 = Plane(origin=(w,0.,0.), normal=(a,-b,z_err))
        return [p1,p2,p3,p4,p5,p6]
