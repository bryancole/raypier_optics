class BeamStop(VTKOptic):
    name = "Beam Stop"
    height = Float #distance between parallel faces
    width = Float #width of parallel faces
    
    _polydata = Property(depends_on=['height', 'width'])
    
    ray_count = Property(Int, depends_on=['intersections_items'])
    ave_ray_length = Float(0.0)
    ave_power = Float(0.0)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('ray_count', style="readonly",height=-25),
                       Item('ave_ray_length', label="path length",
                            style="readonly",height=-25),
                       Item('ave_power', label="incident power",
                            style="readonly",height=-25)
                       )
                       )
    
    @cached_property
    def _get__polydata(self):
        h = self.height/2
        w = self.width/2
        points = [(-w,-h,0),
                  (w,-h,0),
                  (w,h,0),
                  (-w,h,0)]
        cells = [(0,1,2),
                 (2,3,0)
                 ]
        pd = tvtk.PolyData(points=numpy.array(points), 
                           polys=numpy.array(cells))
        return pd
    
    def _get_ray_count(self):
        return len(self.intersections)
    
    def eval_children(self, seg, point, cell_id):
        return []
    
    def update_complete(self):
        try:
            ave = sum(seg.cum_length for seg in self.intersections)/len(self.intersections)
        except ZeroDivisionError:
            ave = 0.0
        self.ave_ray_length = ave
        
#        for seg in self.intersections:
#            print seg.E1_amp, seg.E2_amp
        
        try:
            pwr = sum(self.eval_pwr(seg) for seg in self.intersections)/len(self.intersections)
        except ZeroDivisionError:
            pwr = 0.0
        self.ave_power = pwr
        
    @staticmethod
    def eval_pwr(seg):
        E1 = seg.E1_amp
        E2 = seg.E2_amp
        P1 = (abs(E1)**2)
        P2 = (abs(E2)**2)
        return P1 + P2
    