#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
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

import warnings

try:
    from traits.etsconfig.api import ETSConfig
    target = 'qt4' #was 'wx'
    ETSConfig.toolkit = target
except ImportError: 
    ### If ETS is not present, we don't want to fail immediately, in case the user only want's to use raypier.core.
    warnings.warn(f"Couldn't set ETS toolkit backend to {target}.")

