#!/usr/bin/env python

import numpy as np
import itertools
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import calc_tau

#-- parse field
parser = argparse.ArgumentParser(description=
                                 'Unpack and plot grid data from SuperNu output.')
parser.add_argument('f',type=str,help='Variable data file.')
parser.add_argument('g',type=str,help='Grid/mesh file.')
parser.add_argument('t',type=str,help='Time step file.')
parser.add_argument('wl',type=str,help='Flux wavelength file.')
parser.add_argument('--wavelen',type=float,default=0.0, help='Wavelength to select group to colormap [Ang].')
parser.add_argument('--time',type=float,default=0.0, help='Time to plot [day].')
parser.add_argument('--title',type=str, default='Data', help='Title for plot.')
parser.add_argument('--variable',type=str, default='Variable',
                    help='Label for colorbar.')
parser.add_argument('--plot_optdepth_path',action="store_true",help='Plot optical depth over path.')
parser.add_argument('--coor1',type=float,nargs=2,default=[0.0,0.0], help='Starting point of path [c]')
parser.add_argument('--coor2',type=float,nargs=2,default=[0.1,0.1], help='Ending point of path [c]')

#-- parse arguments
args = parser.parse_args()

#------------------------------------------------------------------------------#
# Load grid de-compression data
#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#
# Load time steps and find index of time nearest to user time.
#------------------------------------------------------------------------------#

#-- obtain the time step from the nearest time
s2day = 86400
tgrid = np.loadtxt(args.t,skiprows=1)
nt = np.size(tgrid) - 1
t = 0.5 * (tgrid[1:] + tgrid[:nt])
it = np.argmin(abs(t - 86400 * args.time))

#------------------------------------------------------------------------------#
# Load flux wavelength grid
#------------------------------------------------------------------------------#

fo = open(args.wl,'r')
fo.readline()
flx_wl = np.array(list(map(float,fo.readline().split())))
fo.close()
ng = np.size(flx_wl) - 1
wlc = 0.5*(flx_wl[1:]+flx_wl[:ng])
ig = np.argmin(abs(wlc - 1e-8*args.wavelen))

#------------------------------------------------------------------------------#
# Decompress spatial gridding of multigroup data
#------------------------------------------------------------------------------#

#-- read grid-based data
data = []
with open(args.f, 'r') as infile:
    l1 = int(it*ng*comp_dims[0])
    l2 = int((it+1)*ng*comp_dims[0])
    print()
    print('l1,l2 = ', l1,l2)
    print()
    for line in itertools.islice(infile, l1, l2):
        data.append(list(map(float,line.split())))

#-- reshape the data and truncate cell padding
data = np.array(data)
data = data.reshape((ng, comp_dims[0]*comp_dims[1]))

#-- using the index map, unpack the data
real_dims.insert(0, ng)
data_full = np.zeros(tuple(reversed(real_dims)))
for k in range(real_dims[3]):
    for j in range(real_dims[2]):
        for i in range(real_dims[1]):
            icell = ijk_cell_map[k,j,i] - 1
            if(icell == ncell - 1 and has_void):
                data_full[k,j,i,:] = 0.0
            else:
                data_full[k,j,i,:] = data[:,icell]


#------------------------------------------------------------------------------#
# Plot colormap of data (todo: extend for 1D and 3D)
#------------------------------------------------------------------------------#

#-- scale x, y velocity by speed of light
c = 2.998e10
x = 0.5 * (xgrid[1:] + xgrid[:real_dims[1]]) / c
y = 0.5 * (ygrid[1:] + ygrid[:real_dims[2]]) / c

#-- reflect the data about the y-axis
x_refl = np.zeros(2 * real_dims[1])
x_refl[:real_dims[1]] = -x[::-1]
x_refl[real_dims[1]:] = x
data_refl = np.zeros((real_dims[3], real_dims[2], 2*real_dims[1], real_dims[0]))
data_refl[:,:,:real_dims[1],:] = data_full[:,:,::-1,:]
data_refl[:,:,real_dims[1]:,:] = data_full

#-- plot
ig = 0
f1 = plt.figure(1)
X, Y = np.meshgrid(x_refl, y)
plt.pcolor(X, Y, data_refl[0,:,:,ig])
plt.colorbar(label=args.variable+' , {:.2f}'.format(1e8*wlc[ig])+' Ang')
plt.xlabel('R Velocity [c]')
plt.ylabel('Z Velocity [c]')
plt.title(args.title+' , {:.2f}'.format(t[it]/86400)+' day')

#-- extract a line-out
x0, y0 = args.coor1[0], args.coor1[1]
x1, y1 = args.coor2[0], args.coor2[1]
lgth = np.hypot(x1-x0, y1-y0)
xi0, yi0 = np.argmin(abs(x_refl-x0)), np.argmin(abs(y-y0))
xi1, yi1 = np.argmin(abs(x_refl-x1)), np.argmin(abs(y-y1))
lgthi = int(np.hypot(xi1-xi0, yi1-yi0))
xi, yi = np.linspace(xi0, xi1, lgthi), np.linspace(yi0, yi1, lgthi)
zi = data_refl[0,yi.astype(np.int), xi.astype(np.int), :]

plt.plot([x0, x1], [y0, y1], 'ro-',label='Path')
plt.legend()

f2 = plt.figure(2)
xl = np.linspace(0.0,lgth,np.size(zi[:,ig]))
plt.plot(xl, zi[:,ig])
plt.xlabel('Path Velocity [c]')
plt.ylabel(args.variable)
plt.title(args.title+' , {:.2f}'.format(t[it]/86400)+' day')

if('output.grd_cap'==args.f and args.plot_optdepth_path):
    #-- calculate optical depth
    nx = np.size(xl)
    taus, dift = calc_tau.calc_tau(t[it],nx,c*xl,ng,flx_wl,zi.T)
    #-- plot optical depth
    xl = np.insert(xl,0,0.0)
    wlc_obs = flx_wl[1:]*np.exp(-xl[nx])*1e8
    xlc = 0.5 * (xl[1:] + xl[:nx])
    VV, WW = np.meshgrid(xlc, wlc_obs)
    #f3 = plt.figure(3)
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    pcm=ax.pcolor(VV,WW,taus[:,:],cmap='RdBu_r',norm=LogNorm(1e-3,1e3))
    plt.xlabel('Path Velocity Differential [c]')
    plt.ylabel('Wavelength [$\AA$]')
    cbar=plt.colorbar(pcm)
    clabel = 'Optical Depth'
    cbar.set_label(clabel,rotation=270,labelpad=15)


#-- show all figures
plt.show()
