#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages

setup(
    name="raytrace",
    version="0.1dev",
    packages=find_packages(),
    scripts = [], #no stand-alone application yet
    
    zip_safe = True, #why not!
    
    install_requires = ["numpy >= 1.1",
        #"VTK",
        #"wxPython", #maybe can avoid an explicit dependancy on wx
        "Traits >= 3.0",
        "Mayavi >= 3.0", #TVTK is distributed as part of Mayavi
        "TraitsGUI >= 3.0"
        ],
        
    package_data = {}, #none, yet
    
    author = "Bryan Cole",
    author_email = "bryan.cole@teraview.com",
    description = "A optical ray-tracing package, for design, optimisation and \
visualisation of mirror/lens systems. Good support for off-axis aspheric \
surfaces is a priority",
    license = "GPL3",
    keywords = "science engineering optics ray-tracing physcics",
    url = "python-raytrace@googlegroups.com",
    #download_url="http://bitbucket.org/bryancole/raytrace/get/"
    )
    