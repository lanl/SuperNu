#!/usr/bin/env python

#-- import libraries
import numpy as np
import argparse
import calc_tau

parser = argparse.ArgumentParser(description='Plot comoving optical depth.')
parser.add_argument('sim_dir',type=str,help='Simulation directory.')
#-- parsing arguments
args = parser.parse_args()

#--------------------------------------------------------------------#
# Load input structure and grid
#--------------------------------------------------------------------#

#-- assumes 1D, with only right-cell-edge velocities
instr = np.loadtxt(args.sim_dir+'/input.str', skiprows=3)
grd_vel = instr[:,0]
nx = np.size(grd_vel)
grd_vel = np.insert(grd_vel, 0, 0.0)
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
    taus[it,:,:], dift[it,:,:] = calc_tau.calc_tau(time[it], nx, grd_vel, ng,
                                                   flx_wl, opac[it,:,:])


#-- save optical depth and supporting data in simulation directory
np.savetxt(args.sim_dir+'/output.grd_tau', np.reshape(taus,(nt*ng,nx)), fmt='%.4e')
np.savetxt(args.sim_dir+'/output.grd_dift', np.reshape(dift,(nt*ng,nx)), fmt='%.4e')
