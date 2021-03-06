'''
Created on 23 Nov 2018

:author: bryan
'''
import pyximport
pyximport.install()

from raypier.core.cmaterials import BaseDispersionCurve

import sqlite3
import pkg_resources
import numpy


GLASS_DATABASE_PATH = "material_data/glass_dispersion_database.db"
MATERIAL_DATABASE = pkg_resources.resource_filename("raypier", GLASS_DATABASE_PATH)


class NondispersiveCurve(BaseDispersionCurve):
    """
    A Dispersion curve for a non-dispersive material with a given refractive index and absorption.
    """
    def __init__(self, refractive_index=1.37, absorption=0.0):
        self._refractive_index=refractive_index
        self._absorption=absorption
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
        
    def __repr__(self):
        return f"<Nondispersion Curve: ri={self._refractive_index}, absorption={self._absorption}>"
        
class FusedSilica(BaseDispersionCurve):
    """
    A Dispersion curve for fused silica.
    """
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
    """
    A Dispersion curve obtained from the materials database (http://refractiveindex.info).
    """
    def __init__(self, name=None, book=None, filename=None, absorption=0.0):
        filters = {"name": name, "book": book, "filename": filename}
        names, vals = list(zip(*[('(%s=?)'%(k,),v) for k,v in list(filters.items())\
                             if v is not None]))
        if not names:
            raise ValueError("No material data identifier given. "
                             "Must give at least one name, book or filename")
        where = " AND ".join(names)
        sql = "select * from dispersion where %s"%(where,)
        
        conn = sqlite3.connect(MATERIAL_DATABASE)
        #print(MATERIAL_DATABASE)
        #print("SQL", sql, vals)
        rows = list(conn.execute(sql, vals))
        if len(rows) > 1:
            raise ValueError("Names given return multiple Material entries.")
        
        row=rows[0]
        fname, name, book, formula_id, wavelen_min, wavelen_max, n_coefs = row[:7]
        
        self._data = row
        
        coefs = numpy.array(row[7:7+n_coefs], dtype=numpy.double)
        
        super(NamedDispersionCurve,self).__init__(formula_id,
                                                  coefs,
                                                  absorption,
                                                  wavelen_min,
                                                  wavelen_max
                                                  )
        
    def __repr__(self):
        data = self._data
        return f"<Named Dispersion: name={data[1]}, formula={data[3]}, coefs={data[6]}>"
        
    @classmethod
    def get_glass_names(cls):
        sql  = "select name from dispersion"
        conn = sqlite3.connect(MATERIAL_DATABASE)
        rows = list(conn.execute(sql))
        return [r[0] for r in rows]
        
        
if __name__=="__main__":
    
    BK7 = NamedDispersionCurve("N-LAK22")
    print(BK7.formula_id, list(BK7.coefs))
    print(BK7.evaluate_n([1.0, 1.5]))
    
    print([a for a in NamedDispersionCurve.get_glass_names() if a.startswith("D")])
        