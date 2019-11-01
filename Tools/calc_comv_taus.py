#!/usr/bin/env python

#-- import libraries
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Plot comoving optical depth.')
parser.add_argument('sim_dir',type=str,help='Simulation directory.')
#-- parsing arguments
args = parser.parse_args()

#-- constants
pc_c = 2.998e10 # [cm/s]

#--------------------------------------------------------------------#
# Load input structure and grid
#--------------------------------------------------------------------#

instr = np.loadtxt(args.sim_dir+'/input.str', skiprows=3)
grd_vel = instr[:,0]
nx = np.size(grd_vel)
fo = open(args.sim_dir+'/output.flx_grid','r')
fo.readline()
flx_wl = np.array(list(map(float,fo.readline().split())))
fo.close()
ng = np.size(flx_wl) - 1
t = np.loadtxt(args.sim_dir+'/output.tsp_time',skiprows=1)
nt = np.size(t) - 1
time = 0.5 * (t[:nt] + t[1:])

#--------------------------------------------------------------------#
# Load opacity vs time step, spatial cell, wavelength group
#--------------------------------------------------------------------#

#-- multigroup absorption opacity [1/cm]
opac = np.loadtxt(args.sim_dir+'/output.grd_cap')
opac = np.reshape(opac,(nt,ng,nx))

#-- optical depth array
taus = np.zeros((nt,ng,nx))
dift = np.zeros((nt,ng,nx))

for it in range(nt):
    for jg in range(ng):
        #-- backtrack ray from rightmost wavelength
        tau = 0.0
        dif = 0.0
        wl0 = flx_wl[ng]
        v0 = grd_vel[nx-1]
        #-- initialize ray counters
        ix = nx-2
        ig = jg
        #-- calculate optical depth from surface to center
        while (ix >= 0 and ig >= 0):
            dvx = v0-grd_vel[ix]
            dvg = pc_c*(wl0/flx_wl[ig] - 1.0)
            #-- accrue optical depth
            tau += time[it] * min(dvx,dvg)*opac[it,ig,ix]
            dif += opac[it,ig,ix] * (time[it] * min(dvx,dvg))**2 / pc_c
            #-- update ray parameter
            if (dvx < dvg):
                # set new ray state
                v0 = grd_vel[ix]
                wl0 /= (1.0 + dvx/pc_c)
                # set the tau for this cell
                taus[it,jg,ix] = tau
                dift[it,jg,ix] = dif
                # decrement cell index
                ix -= 1
            else:
                # set new ray state
                wl0 = flx_wl[ig]
                v0 -= dvg
                # decrement group index
                ig -= 1


#-- save optical depth and supporting data in simulation directory
np.savetxt(args.sim_dir+'/output.grd_tau', np.reshape(taus,(nt*ng,nx)), fmt='%.4e')
np.savetxt(args.sim_dir+'/output.grd_dift', np.reshape(dift,(nt*ng,nx)), fmt='%.4e')
