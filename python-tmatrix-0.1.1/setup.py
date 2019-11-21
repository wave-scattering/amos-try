#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2010  Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of python-tmatrix
#
#    This program is free software: you can redistribute it and/or modify
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

__version__ = '0.1.1'
__title__ = 'T-matrix codes by M.I. Mishchenko, L.D. Travis, and D.W. Mackowski'
__mod__ = 'python-tmatrix'
__author__ = 'Ovidio Peña Rodríguez'
__email__ = 'ovidio@bytesfall.com'
__url__ = 'http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html'

from numpy.distutils.core import setup, Extension
import numpy as np

setup(name = __mod__,
      version = __version__,
      description = __title__,
      long_description="""A Python implementation of T-matrix codes by Michael I. Mishchenko, Larry D. Travis, \
and Daniel W. Mackowski (the used versions of the code employ DOUBLE PRECISION \
variables) for the following cases: \
 * The calculation of the amplitude and phase matrices for a particle with an axially \
   symmetric shape in fixed or randomly oriented positions. The calculations are \
   applicable to spheroids, finite circular cylinders, Chebyshev particles, and \
   generalized Chebyshev particles (distorted water drops). \
 * The calculation of light scattering by randomly oriented two-sphere clusters with \
   touching or separated components. \
 * The calculation of the efficiency factors and scattering matrix for an ensemble of \
   spheres in fixed or randomly oriented positions.""",
      author = __author__,
      author_email = __email__,
      maintainer = __author__,
      maintainer_email = __email__,
      keywords = ['Light scattering', 'Axisymmetric particles', 'Efficiency factors', 'Cross-sections'],
      url = __url__,
      license = 'GPL',
      platforms = 'any',
      packages = ['tmatrix'],
      ext_package = 'tmatrix',
      ext_modules = [Extension('_tmfixed', ['src/tmfixed.pyf', 'src/lpd.f', 'src/ampld.lp.f'],
                         f2py_options = ['--include_paths', np.get_include()]),
                     Extension('_tmrandom', ['src/tmrandom.pyf', 'src/lpd.f', 'src/tmd.lp.f'],
                         f2py_options = ['--include_paths', np.get_include()]),
                     Extension('_tmbisphere', ['src/tmbisphere.pyf', 'src/bisphere.f'],
                         f2py_options = ['--include_paths', np.get_include()]),
                     Extension('_tmnsfixed', ['src/tmnsfixed.pyf', 'src/scsmfo1b.f'],
                         f2py_options = ['--include_paths', np.get_include()]),
                     Extension('_tmnsrandom', ['src/tmnsrandom.pyf', 'src/scsmtm1.f'],
                         f2py_options = ['--include_paths', np.get_include()])]
)

