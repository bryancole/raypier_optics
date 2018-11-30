'''
Created on 23 Nov 2018

@author: bryan
'''
import pyximport
pyximport.install()

from raytrace.cmaterials import BaseDispersionCurve

import sqlite3
import pkg_resources
import numpy


GLASS_DATABASE_PATH = "material_data/glass_dispersion_database.db"
MATERIAL_DATABASE = pkg_resources.resource_filename("raytrace", GLASS_DATABASE_PATH)


class NondispersiveCurve(BaseDispersionCurve):
    def __init__(self, refractive_index=1.37, absorption=0.0):
        formula_id=0
        coefs = numpy.array([refractive_index,])
        wavelen_min=0.0
        wavelen_max=1000000.0
        super(NondispersiveCurve,self).__init__(formula_id,
                                                  coefs,
                                                  absorption,
                                                  wavelen_min,
                                                  wavelen_max
                                                  )
        
class FusedSilica(BaseDispersionCurve):
    def __init__(self, absorption=0.0):
        formula_id=1
        coefs = numpy.array([0.0,
                             0.6961663,
                             0.0684043,
                             0.4079426,
                             0.1162414,
                             0.8974794,
                             9.896161])
        wavelen_min=0.21
        wavelen_max=6.7
        super(FusedSilica,self).__init__(formula_id,
                                                  coefs,
                                                  absorption,
                                                  wavelen_min,
                                                  wavelen_max
                                                  )


class NamedDispersionCurve(BaseDispersionCurve):
    def __init__(self, name=None, book=None, filename=None, absorption=0.0):
        filters = {"name": name, "book": book, "filename": filename}
        names, vals = zip(*[('(%s=?)'%(k,),v) for k,v in filters.items()\
                             if v is not None])
        if not names:
            raise ValueError("No material data identifier given. "
                             "Must give at least one name, book or filename")
        where = " AND ".join(names)
        sql = "select * from dispersion where %s"%(where,)
        
        conn = sqlite3.connect(MATERIAL_DATABASE)
        print MATERIAL_DATABASE
        print "SQL", sql, vals
        rows = list(conn.execute(sql, vals))
        if len(rows) > 1:
            raise ValueError("Names given return multiple Material entries.")
        
        row=rows[0]
        fname, name, book, formula_id, wavelen_min, wavelen_max, n_coefs = row[:7]
        coefs = numpy.array(row[7:7+n_coefs], dtype=numpy.double)
        
        super(NamedDispersionCurve,self).__init__(formula_id,
                                                  coefs,
                                                  absorption,
                                                  wavelen_min,
                                                  wavelen_max
                                                  )
        
        
if __name__=="__main__":
    
    BK7 = NamedDispersionCurve("N-LAK22")
    print BK7.formula_id, list(BK7.coefs)
    print BK7.evaluate_n([1.0, 1.5])
        