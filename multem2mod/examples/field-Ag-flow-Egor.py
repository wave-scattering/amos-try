#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This file is part of python-scattnlay
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
#    The only additional remark is that we expect that all publications
#    describing work using this software, or all commercial products
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This test case calculates the electric field in the
# E-k plane, for an spherical Ag nanoparticle.

from fieldplot import fieldplot
#from polarplot import polarplot
from scattnlay import fieldnlay
from scattnlay import scattnlay
import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import scattnlay
# # a)
#WL=400 #nm
#core_r = WL/20.0
#epsilon_Ag = -2.0 + 10.0j

# # b)
#WL=400 #nm
#core_r = WL/20.0
#epsilon_Ag = -2.0 + 1.0j
crossplane='XZ'
#crossplane='YZ'
#crossplane='XY'

# c)
#WL=548 #nm
#WL=1000 #nm
#WL=750 #nm
#WL=688 #nm

factor=1.5
# core_r = 120
# index_Ag = 4.0
# factor=12.33
# index_Ag = 1.001
nm = 1.0
WL_units='nm'

# WL=456.33
# core_r = 90
# index_Ag = 4.615+0.131j

# WL=1305.
# core_r = 900.
# index_Ag = 2.5

# WL=434.6
# core_r = 300.
# index_Ag = 2.5


# core_r = 400.
# kr = 7.72745 # hybrid state
# # kr = 7.09 #left border of wide band
# # kr = 8.23 #right border of wide band
# # kr = 13.5 #
# WL = 2*np.pi*core_r/kr
# index_Ag = 1.41125

core_r = 400.
# kr = 4.48175 # hybrid state
kr = 1.476
WL = 2*np.pi*core_r/kr
# index_Ag = 1.7222
# index_Ag = 4.532
index_Ag = 4.0

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)

isFieldStrealines = True

# isFieldStrealines = False
#flow_total = 9
# flow_total = 21
# flow_total = 41
flow_total = 0

# Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
field_to_plot='Eabs'
#field_to_plot='Pabs'
#field_to_plot='angleEx'

# ["Eabs","Habs"]
#   [480, 570, 410, 370, 444, 531,375.5, 436.1, 327]
#   [-1,1,2,3,4,5]:
#   [0,1]
# ["XZ","YZ","XY", "XYZ"]

#plot_params = [
#     #E2
#     [ 0,"Eabs", 900, -1, "XZ", True]
#     ,[ 0,"Eabs", 900, 1, "XZ", True]
#     ,[ 0,"Eabs", 735, -1, "XZ", True]
   #     ,[ 0,"Eabs", 735, 1, "XZ", True]

#     #M2
#     ,[ 1,"Habs", 950, -1, "YZ", True]
#     ,[ 1,"Habs", 950, 1, "YZ", True]
#     ,[ 1,"Habs", 479, -1, "YZ", True]
#     ,[ 1,"Habs", 479, 1, "YZ", True]
#     #E4
#     ,[ 0,"Eabs", 600, -1, "XZ", True]
#     ,[ 0,"Eabs", 600, 2, "XZ", True]
#     ,[ 0,"Eabs", 545, -1, "XZ", True]
#     ,[ 0,"Eabs", 545, 2, "XZ", True]
#     #M4
#     ,[ 1,"Habs", 500, -1, "YZ", True]
#     ,[ 1,"Habs", 500, 2, "YZ", True]
#     ,[ 1,"Habs", 397.2, -1, "YZ", True]
#     ,[ 1,"Habs", 397.2, 2, "YZ", True]
    # #E8
    # ,[ 0,"Eabs", 480, -1, "XZ", True]
    # ,[ 0,"Eabs", 480, 3, "XZ", True]
    # ,[ 0,"Eabs", 444, -1, "XZ", True]
    # ,[ 0,"Eabs", 444, 3, "XZ", True]
    # #M8
    # ,[ 1,"Habs", 570, -1, "YZ", True]
    # ,[ 1,"Habs", 570, 3, "YZ", True]
    # ,[ 1,"Habs", 531, -1, "YZ", True]
    # ,[ 1,"Habs", 531, 3, "YZ", True]
    # #E16
    # ,[ 0,"Eabs", 410, -1, "XZ", True]
    # ,[ 0,"Eabs", 410, 4, "XZ", True]
    # ,[ 0,"Eabs", 375.5, -1, "XZ", True]
    # ,[ 0,"Eabs", 375.5, 4, "XZ", True]
    # #M16
    # ,[ 1,"Habs", 480, -1, "YZ", True]
    # ,[ 1,"Habs", 480, 4, "YZ", True]
    # ,[ 1,"Habs", 436.1, -1, "YZ", True]
    # ,[ 1,"Habs", 436.1, 4, "YZ", True]
    # #E32
    # ,[ 0,"Eabs", 370, -1, "XZ", True]
    # ,[ 0,"Eabs", 370, 5, "XZ", True]
    # ,[ 0,"Eabs", 327, -1, "XZ", True]
    # ,[ 0,"Eabs", 327, 5, "XZ", True]
    # #647 E2
    # ,[ 0,"Eabs", 647, -1, "XZ", True]
    # ,[ 0,"Eabs", 647, 1, "XZ", True]
    #516.1 M2
    #[ 1,"Habs", 516.1, -1, "YZ", True]
    # [ 1,"Habs", 516.1, 1, "YZ", True]
    # #511.1 E4
    # ,[ 0,"Eabs", 511.1, -1, "XZ", True]
    # ,[ 0,"Eabs", 511.1, 2, "XZ", True]
    # #427.4 M4
    # ,[ 1,"Habs", 427.4, -1, "YZ", True]
    # ,[ 1,"Habs", 427.4, 2, "YZ", True]
    # #424.2 E8
    # ,[ 0,"Eabs", 424.2, -1, "XZ", True]
    # ,[ 0,"Eabs", 424.2, 3, "XZ", True]

    #kr = 2 pi N_host r /lambda
   #  [ 1,"Habs", 2*np.pi*core_r/1.117, 1, "YZ", True]
   # ,[ 0,"Eabs", 2*np.pi*core_r/1.117, 1, "XZ", True]
   # ,[ 1,"Habs", 2*np.pi*core_r/1.945, 1, "YZ", True]
   # ,[ 0,"Eabs", 2*np.pi*core_r/1.945, 1, "XZ", True]
   #  #kr = 2 pi N_host r /lambda
   #  [ 1,"Habs", 671.19017, -1, "YZ", True]
   # ,[ 0,"Eabs", 671.19017, -1, "XZ", True]
#]
# plot_params = [
#     #E8
#     # [ 0,"Eabs", 508, -1, "XZ", True]
#     # ,[ 0,"Eabs", 508, 2, "XZ", True]
#     [ 0,"Eabs", 521, -1, "XZ", True]
#     ,[ 0,"Eabs", 521, 2, "XZ", True]
#     # ,[ 0,"Eabs", 548, -1, "XZ", True]
#     # ,[ 0,"Eabs", 548, 2, "XZ", True]
#      ,[ 0,"Eabs", 700, -1, "XZ", True]
#      ,[ 0,"Eabs", 700, 2, "XZ", True]
#      ,[ 0,"Eabs", 647, -1, "XZ", True]
#      ,[ 0,"Eabs", 647, 1, "XZ", True]
#      ,[ 0,"Eabs", 1025, -1, "XZ", True]
#      ,[ 0,"Eabs", 1025, 1, "XZ", True]
#     # ,[ 1,"Habs", 658, -1, "YZ", True]
#     # ,[ 1,"Habs", 658, 1, "YZ", True]
#     # ,[ 1,"Habs", 523, -1, "YZ", True]
#     # ,[ 1,"Habs", 523, 1, "YZ", True]
#     # ,[ 0,"Eabs", 550, 1, "XZ", True]
#     # ,[ 0,"Eabs", 550, 2, "XZ", True]
#     # ,[ 0,"Eabs", 550, 3, "XZ", True]
#     ]
#plot_params = [ [ 0,"Eabs", 405+i*10, 5, "XZ", True] for i in range(5)]
#plot_params = [ [ 1,"Habs", 500+i*1, 2, "YZ", True] for i in range(30)]
#plot_params = [
    # E2
    #[ 0,"Eabs", 735, -1, "XZ", True]
    # ,[ 0,"Eabs", 735, 1, "XZ", True]
    # ,[ 0,"Eabs", 647, -1, "XZ", True]
    # ,[ 0,"Eabs", 647, 1, "XZ", True]
    # ,[ 0,"Eabs", 1000, -1, "XZ", True]
    # ,[ 0,"Eabs", 1000, 1, "XZ", True]
    # # M2
    # ,[ 1,"Habs", 997, -1, "YZ", True]
    # ,[ 1,"Habs", 997, 1, "YZ", True]
    # ,[ 1,"Habs", 479, -1, "YZ", True]
    # ,[ 1,"Habs", 479, 1, "YZ", True]
    # ,[ 1,"Habs", 516.1, -1, "YZ", True]
    # ,[ 1,"Habs", 516.1, 1, "YZ", True]
    # ,[ 1,"Habs", 1000, -1, "YZ", True]
    # ,[ 1,"Habs", 1000, 1, "YZ", True]
    # # E4
    # ,[ 0,"Eabs", 548, -1, "XZ", True]
    # ,[ 0,"Eabs", 548, 2, "XZ", True]
    # ,[ 0,"Eabs", 350, -1, "XZ", True]
    # ,[ 0,"Eabs", 350, 2, "XZ", True]
    # ,[ 0,"Eabs", 511.1, -1, "XZ", True]
    # ,[ 0,"Eabs", 511.1, 2, "XZ", True]
    # ,[ 0,"Eabs", 680, -1, "XZ", True]
    # ,[ 0,"Eabs", 680, 2, "XZ", True]
    # # M4
    # ,[ 1,"Habs", 688, -1, "YZ", True]
    # ,[ 1,"Habs", 688, 2, "YZ", True]
    # ,[ 1,"Habs", 397.2, -1, "YZ", True]
    # ,[ 1,"Habs", 397.2, 2, "YZ", True]
    # ,[ 1,"Habs", 427.4, -1, "YZ", True]
    # ,[ 1,"Habs", 427.4, 2, "YZ", True]
    # ,[ 1,"Habs", 680, -1, "YZ", True]
    # ,[ 1,"Habs", 680, 2, "YZ", True]
    # # E8
    # ,[ 0,"Eabs", 443.2, -1, "XZ", True]
    # ,[ 0,"Eabs", 443.2, 3, "XZ", True]
    # ,[ 0,"Eabs", 300.5, -1, "XZ", True]
    # ,[ 0,"Eabs", 300.5, 3, "XZ", True]
    # ,[ 0,"Eabs", 424.2, -1, "XZ", True]
    # ,[ 0,"Eabs", 424.2, 3, "XZ", True]
    # ,[ 0,"Eabs", 600, -1, "XZ", True]
    # ,[ 0,"Eabs", 600, 3, "XZ", True]
    # # M8
    # ,[ 1,"Habs", 531, -1, "YZ", True]
    # ,[ 1,"Habs", 531, 3, "YZ", True]
    # ,[ 1,"Habs", 511, 3, "YZ", True]
    # ,[ 1,"Habs", 491, 3, "YZ", True]
    # ,[ 1,"Habs", 471, 3, "YZ", True]
    # ,[ 1,"Habs", 451, 3, "YZ", True]
    # ,[ 1,"Habs", 431, 3, "YZ", True]
    # ,[ 1,"Habs", 336.6, -1, "YZ", True]
    # ,[ 1,"Habs", 336.6, 3, "YZ", True]
    # ,[ 1,"Habs", 600, -1, "YZ", True]
    # ,[ 1,"Habs", 600, 3, "YZ", True]
    # # E16
    # ,[ 0,"Eabs", 375.5, -1, "XZ", True]
    # ,[ 0,"Eabs", 375.5, 4, "XZ", True]
    # ,[ 0,"Eabs", 264, -1, "XZ", True]
    #  ,[ 0,"Eabs", 264, 4, "XZ", True]
    # ,[ 0,"Eabs", 500, -1, "XZ", True]
    # ,[ 0,"Eabs", 500, 4, "XZ", True]
    # M16
    # ,[ 1,"Habs", 436.1, -1, "YZ", True]
    # ,[ 1,"Habs", 436.1, 4, "YZ", True]
    #,[ 1,"Habs", 436.1, -1, "YZ", True]
#    ,[ 1,"Habs", 436.1, 4, "YZ", True]
    #,[ 1,"Habs", 435.1, -1, "YZ", True]
#    ,[ 1,"Habs", 435.1, 4, "YZ", True]
    # ,[ 1,"Habs", 380, -1, "YZ", True]
    # ,[ 1,"Habs", 380, 4, "YZ", True]
   # ,[ 1,"Habs", 292.9, -1, "YZ", True]
    # ,[ 1,"Habs", 292.9, 4, "YZ", True]
    # ,[ 1,"Habs", 500, -1, "YZ", True]
    # ,[ 1,"Habs", 500, 4, "YZ", True]
    # # E32
    # ,[ 0,"Eabs", 327, -1, "XZ", True]
    # ,[ 0,"Eabs", 327, 5, "XZ", True]
    # ,[ 0,"Eabs", 236.6, -1, "XZ", True]
    # ,[ 0,"Eabs", 236.6, 5, "XZ", True]
    # [ 0,"Eabs", 2036.6, -1, "XZ", True]
    # ,[ 0,"Eabs", 2036.6, 1, "XZ", True]
    # [ 0,"Eabs", 644.7, -1, "XZ", True]
    # ,[ 0,"Eabs", 644.7, 1, "XZ", True]
    # [ 0,"Eabs", 688.9, -1, "XZ", True]
    # ,[ 0,"Eabs", 688.9, 1, "XZ", True]
    #  ,[ 1,"Eabs", 688.9, 1, "XZ", True]
    #  ,[ -1,"Eabs", 688.9, 1, "XZ", True]
    #  ,[ 0,"Eabs", 688.9, 2, "XZ", True]
    #  ,[ 1,"Eabs", 688.9, 2, "XZ", True]
    #  ,[ -1,"Eabs", 688.9, 2, "XZ", True]

    # ,[ 0,"Eabs", 671, -1, "XZ", True]
    # ,[ 0,"Eabs", 671, 1, "XZ", True]
    #  ,[ 1,"Eabs", 671, 1, "XZ", True]
    #  ,[ -1,"Eabs", 671, 1, "XZ", True]
    #  ,[ 0,"Eabs", 671, 2, "XZ", True]
    #  ,[ 1,"Eabs", 671, 2, "XZ", True]
    #  ,[ -1,"Eabs", 671, 2, "XZ", True]

    # ,[ 0,"Eabs", 648, -1, "XZ", True]
    # ,[ 0,"Eabs", 648, 1, "XZ", True]
    #  ,[ 1,"Eabs", 648, 1, "XZ", True]
    #  ,[ -1,"Eabs", 648, 1, "XZ", True]
    #  ,[ 0,"Eabs", 648, 2, "XZ", True]
    #  ,[ 1,"Eabs", 648, 2, "XZ", True]
    #  ,[ -1,"Eabs", 648, 2, "XZ", True]


    # ,[ 0,"Eabs", 610.2, -1, "XZ", True]
    # ,[ 0,"Eabs", 610.2, 1, "XZ", True]
    #  ,[ 1,"Eabs", 610.2, 1, "XZ", True]
    #  ,[ -1,"Eabs", 610.2, 1, "XZ", True]
    #  ,[ 0,"Eabs", 610.2, 2, "XZ", True]
    #  ,[ 1,"Eabs", 610.2, 2, "XZ", True]
    #  ,[ -1,"Eabs", 610.2, 2, "XZ", True]



    # [ 0,"Eabs", 660.5, -1, "XZ", True]
    # ,[ 0,"Eabs", 660.5, 1, "XZ", True]
    # ,[ 0,"Eabs", 800, -1, "XZ", True]
    # ,[ 0,"Eabs", 800, 1, "XZ", True]
    # ,[ 0,"Eabs", 384.5, -1, "XZ", True]
    # ,[ 0,"Eabs", 384.5, 1, "XZ", True]
    # ,[ 0,"Eabs", 406.3, -1, "XZ", True]
    # ,[ 0,"Eabs", 406.3, 1, "XZ", True]
    # [ 1,"Habs", 426.1, -1, "YZ", True]
    # ,[ 1,"Habs", 426.1, 4, "YZ", True]
    # ,[ 1,"Habs", 2036.6, 2, "YZ", True]
    # ,[ 1,"Habs", 2036.6, 3, "YZ", True]
    # ,[ 0,"Eabs", 500, -1, "XZ", True]
    # ,[ 0,"Eabs", 500, 5, "XZ", True]
#]
# plot_params = [
# #     #E2

#     [ 0,"Eabs", WL, -1, "XZ", True]
#    ,[ 0,"Eabs", WL, 1, "XZ", True]
# #     [ 0,"Eabs", 950, 1, "XZ", True]
# #    # ,[ 0,"Eabs", 371.4, 1, "XZ", True]
# #    # ,[ 0,"Eabs", 390, 1, "XZ", True]
# #    # ,[ 0,"Eabs", 648.9, 1, "XZ", True]
# #    # ,[ 0,"Eabs", 669.7, 1, "XZ", True]
# # #     #M2
# #      ,[ 1,"Habs", 950, 1, "YZ", True]
# #      ,[ 1,"Habs", 509.4475, 1, "YZ", True]
#      ,[ 1,"Habs", WL, -1, "YZ", True]
#      ,[ 1,"Habs", WL, 1, "YZ", True]
# # #     #E4
#      ,[ 0,"Eabs", WL, 2, "XZ", True]
# #    ,[ 0,"Eabs", 600, 2, "XZ", True]
# #    ,[ 0,"Eabs", 491.8345, 2, "XZ", True]
# # #     #M4
#      ,[ 1,"Habs", WL, 2, "YZ", True]
# #    ,[ 1,"Habs", 600, 2, "YZ", True]
# #    ,[ 1,"Habs", 414.0485, 2, "YZ", True]
# #     # #E8
#      ,[ 0,"Eabs", WL, 3, "XZ", True]
# #     ,[ 0,"Eabs", 600, 3, "XZ", True]
# #     ,[ 0,"Eabs", 389.8564, 3, "XZ", True]
# #     # #M8
#      ,[ 1,"Habs", WL, 3, "YZ", True]
# #     ,[ 1,"Habs", 600, 3, "YZ", True]
# #     ,[ 1,"Habs", 374.9290, 3, "YZ", True]
# #     # #E16
#      ,[ 0,"Eabs", WL, 4, "XZ", True]
# #    ,[ 0,"Eabs", 500, 4, "XZ", True]
# #    ,[ 0,"Eabs", 350.0382, 4, "XZ", True]
# #     # #M16
#      ,[ 1,"Habs", WL, 4, "YZ", True]
# #     ,[ 1,"Habs", 500, 4, "YZ", True]
# #     ,[ 1,"Habs", 319.8906, 4, "YZ", True]
# #     # #E32
#      ,[ 0,"Eabs", WL, 5, "XZ", True]
# #     ,[ 0,"Eabs", 500, 5, "XZ", True]
# #     ,[ 0,"Eabs", 313.7671, 5, "XZ", True]
# ]

plot_params = [


    [ 0,"Eabs", WL, -1, "XZ", isFieldStrealines]
   ,[ 1,"Habs", WL, -1, "YZ", isFieldStrealines]
   ,[ 0,"Eabs", WL, 2, "XZ", isFieldStrealines]
   ,[ 1,"Habs", WL, 2, "YZ", isFieldStrealines]
   # ,[ -1,"Pabs", WL, -1, "XYZ", isFieldStrealines]


#     #     #E2
#      ,[ 0,"Eabs", WL, 1, "XZ", isFieldStrealines]
# # # #     #M2
#      ,[ 1,"Habs", WL, 1, "YZ", isFieldStrealines]
# # #     #E4
#      ,[ 0,"Eabs", WL, 2, "XZ", isFieldStrealines]
# # #     #M4
#      ,[ 1,"Habs", WL, 2, "YZ", isFieldStrealines]
# #     # #E8
#      ,[ 0,"Eabs", WL, 3, "XZ", isFieldStrealines]
# #     # #M8
#      ,[ 1,"Habs", WL, 3, "YZ", isFieldStrealines]
# #     # #E16
#      ,[ 0,"Eabs", WL, 4, "XZ", isFieldStrealines]
# #     # #M16
#      ,[ 1,"Habs", WL, 4, "YZ", isFieldStrealines]
# #     # #E32
#      ,[ 0,"Eabs", WL, 5, "XZ", isFieldStrealines]
]

# dipole - 1 -- 2
# quad   - 2 -- 4
# octo   - 3 -- 8
# hex    - 4 -- 16
# 32     - 5 -- 32

# npts = 150
# ext = ".png"
npts = 351
ext = ".pdf"
i=0
font = {'family' : 'monospace',
            # 'weight' : 'bold',
            'size'   : '15'}
mp.rc('font', **font)
mp.rcParams['axes.linewidth'] = 2
#mp.rcParams['grid.linewidth'] = 1.5
mp.rcParams['legend.fontsize'] = 15
# mp.rcParams['axes.labelweight'] = 'bold'
#mp.rcParams['lines.linewidth'] = 1.5

for mode_type, field_to_plot, WL, mode_n,  crossplane, isStream in plot_params :
    i = i+1
    x[0] = 2.0*np.pi*core_r/WL#/4.0*3.0
    m[0] = index_Ag/nm
    # comment='bulk-NP-WL'+str(WL)+WL_units
    comment='bulk-NP-kr'+str(x[0])+"-index"+str(m[0])
    fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
    fig.tight_layout()
    minlength=0.25
    if i>50: minlength=0.3
    fieldplot(fig, axs, x,m, WL, comment, WL_units, crossplane,
              field_to_plot, npts, factor, flow_total,
              subplot_label=' ',is_flow_extend=False,
              mode_n= mode_n, mode_type=mode_type,
              isStream = isStream,
              minlength=minlength,
              inner_only = False
    )
    fig.subplots_adjust(hspace=0.3, wspace=-0.1)
    if mode_n == -1:
        mode = "All"
        mt=""
    else:
        mode = str(2**mode_n)
        mt = "t"+str(mode_type)
    st = ""
    if isStream == True: st = "_stream"
    print(mode_type, field_to_plot, WL, mode_n,  crossplane, isStream)
    print(mode, mt)
    plt.savefig("Egor2/"+"%02d"%(i)+
                comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-"+crossplane+"-"
                +field_to_plot+"-mode"+mode+mt+st+ext,pad_inches=0.02, bbox_inches='tight')
    plt.clf()
    plt.close()
    # polarplot(x,m,"Egor2/"+"%02d"%(i)+"-TSCS-"+
    #             comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-mode"+mode+ext,
    #           mode_n=mode_n)
