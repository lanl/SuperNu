# © 2023. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
# Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
# Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
# National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
# The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
# irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
# copies to the public, perform publicly and display publicly, and to permit others to do so.
#!/usr/bin/env python

import numpy as np
import argparse
import matplotlib.pyplot as plt

#-- parse field
parser = argparse.ArgumentParser(description=
                                 'Unpack and plot grid data from SuperNu output.')
parser.add_argument('f',type=str,help='Variable data file.')
parser.add_argument('g',type=str,help='Grid/mesh file.')
parser.add_argument('t',type=str,help='Time step file.')
parser.add_argument('--time',type=float,default=0.0, help='Time to plot [day].')
parser.add_argument('--title',type=str, default='Data', help='Title for plot.')
parser.add_argument('--variable',type=str, default='Variable',
                    help='Label for colorbar.')

#-- parse arguments
args = parser.parse_args()

#-- open grid file
fo = open(args.g,'r')

#-- read grid headers
hd1 = fo.readline().split()
hd2 = fo.readline().split()
hd3 = fo.readline().split()

#-- read x,y,z grids
s_xgrid = fo.readline().split()
s_ygrid = fo.readline().split()
s_zgrid = fo.readline().split()

#-- read (i,j,k) -> cell index map
ijk_cell_map = []
for line in fo:
    ijk_cell_map.append([int(s) for s in line.split()])

#-- close file
fo.close()

#-- extract real and compression dimensions from headers
geom = int(hd1[1])
real_dims = [int(hd2[1]), int(hd2[2]), int(hd2[3])]
comp_dims = [int(hd3[2]), int(hd3[3])]

#-- convert x,y,z grids to numpy array of floats
xgrid = np.array([float(sx) for sx in s_xgrid])
ygrid = np.array([float(sy) for sy in s_ygrid])
zgrid = np.array([float(sz) for sz in s_zgrid])

#-- reshape index map
ijk_cell_map = np.array(ijk_cell_map)
ijk_cell_map = ijk_cell_map.reshape(tuple(reversed(real_dims)))

#-- obtain number of cells
ncell = np.amax(np.unique(ijk_cell_map))

#-- determine if there was more than one zero-mass cell
has_void = False
if(ncell < real_dims[0] * real_dims[1] * real_dims[2]): has_void = True

#-- read grid-based data
data = np.loadtxt(args.f)

#-- get number of time steps
nt = np.size(data, 0) // comp_dims[0]

#-- reshape the data and truncate cell padding
data = data.reshape((nt, comp_dims[0]*comp_dims[1]))
data = data[:,:ncell]

#-- obtain the time step from the nearest time
s2day = 86400
tgrid = np.loadtxt(args.t,skiprows=1)
t = 0.5 * (tgrid[1:] + tgrid[:nt])
it = np.argmin(abs(t - 86400 * args.time))

#-- using the index map, unpack the data
data_full = np.zeros(tuple(reversed(real_dims)))
for k in range(real_dims[2]):
    for j in range(real_dims[1]):
        for i in range(real_dims[0]):
            icell = ijk_cell_map[k,j,i] - 1
            if(icell == ncell - 1 and has_void):
                data_full[k,j,i] = 0.0
            else:
                data_full[k,j,i] = data[it,icell]


#-- scale x, y velocity by speed of light
c = 2.998e10
x = 0.5 * (xgrid[1:] + xgrid[:real_dims[0]]) / c
y = 0.5 * (ygrid[1:] + ygrid[:real_dims[1]]) / c

#-- reflect the data about the y-axis
x_refl = np.zeros(2 * real_dims[0])
x_refl[:real_dims[0]] = -x[::-1]
x_refl[real_dims[0]:] = x
data_refl = np.zeros((real_dims[2], real_dims[1], 2*real_dims[0]))
data_refl[:,:,:real_dims[0]] = data_full[:,:,::-1]
data_refl[:,:,real_dims[0]:] = data_full

#-- plot
X, Y = np.meshgrid(x_refl, y)
plt.pcolor(X, Y, data_refl[0,:,:])
plt.colorbar(label=args.variable)
plt.xlabel('R Velocity [c]')
plt.ylabel('Z Velocity [c]')
plt.title(args.title+' , '+str(t[it]/86400)+' day')
plt.show()
