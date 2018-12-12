#!/usr/bin/python

#-- libraries
import argparse
import numpy as np
#import matplotlib.pyplot as plt
#-- custom
import elem_data

#-- functions
def masshelper(q):
    "antiderivative: int (1-q^2)^3 q^2 dq"
    return q**3/3.-3.*q**5/5.+3.*q**7/7.-q**9/9.

beta = 3
def plaw_masshelper(q):
    "power law mass distribution, with radius power = -3"
    ret_val = q**(-beta)
    return ret_val

#-- constants
rgas = 8.3145e7 # ideal gas constant [erg/mol/K]
c = 2.997e10 # light speed [cm/s]
msol = 1.989e33 # solar mass [g]
arad = 7.5657e-15 # radiation constant [erg/cc/K^4]

#-- parse field
parser = argparse.ArgumentParser(description='1 element, 1D KN input.str')
parser.add_argument('--t',type=float,help='Light curve start time [s].')
parser.add_argument('--nr',type=int,help='Number of radial cells.')
parser.add_argument('--el',type=str,help='Symbol of selected element.')
parser.add_argument('--mej',type=float,help='Ejecta mass.')
parser.add_argument('--vmax',type=float,help='Max ejecta velocity [c].')
parser.add_argument('--outputname',type=str,help='Structure name')
parser.add_argument('--save',action='store_true',help='Save structure')

#-- parse arguments
args = parser.parse_args()

#-- defaults
if(args.t==None): t = 1e4
else: t=args.t # free expansion time [s]
if(args.nr==None): nr = 64
else: nr=args.nr # number of radial cells
if(args.el==None): el = 'Pd' # Paladium
else: el = args.el
if(args.mej==None): mej = 1.4e-2*msol
else: mej = args.mej*msol
if(args.vmax==None): vmax = 0.25*c
else: vmax = args.vmax*c

#-- reference constants
t0 = 10.0
temp0 = 5.7e7
ye0 = 0.39

#-- helper variables
r0=vmax*t0 # initial radius [cm]
rho0=315*mej/(64*np.pi*r0**3) # 10 s max density [g/cc]

#-- uniform velocity grid
vv = np.linspace(0.,vmax,nr+1)

#-- cell center
vm = .5*(vv[1:]+vv[:nr])
#-- spatial grid
xx = t*vv
xm = .5*(xx[1:]+xx[:nr])
#-- cell volumes
vol = 4*np.pi*(xx[1:]**3-xx[:nr]**3)/3
#-- mass integral per cell
mass=4*np.pi*rho0*(vmax*t0)**3*(masshelper(vv[1:]/vmax)-masshelper(vv[:nr]/vmax))
#-- check summ of mass over ccells
print('sum(mass) [Msol]: ',np.sum(mass)/msol)

#-- normalize density, mass
hlp = mej/np.sum(mass)
print('mass discrepancy ratio: ',hlp)

#-- temperature
temp = temp0*(t0/t)*(1.-(vm/vmax)**2)

#-- print mass-averaged velocity
print('mass averaged velocity [c] = ', np.sum(mass*vm)/(c*np.sum(mass)))

#-- locate selected element in elem_data
iz = 0
elemd = []
for elem in elem_data.elem_info:
    iz+=1
    #-- assuming atomic number sorted and no isotopes in elem_info
    if(el in elem):
        elemd.append(elem[1])
        elemd.append(iz)
        elemd.append(elem[0])

#-- if element strings are files, load them
if ('.dat' in el):
    elsd = np.loadtxt(el)
else:
    elsd = np.array([[elem_data.eldict[el],1.0]])
#-- re-normalize
elsd[:,1] /= np.sum(elsd[:,1])

#-- save
if(args.save):
    print('saving intput.str ...')
    nelem=np.size(elsd,0)
    #-- raw structure array
    cye=np.tile(ye0,(nr))
    massfr=np.tile(elsd[:,1],(nr,1))
    ncol=nelem + 4
    raw=np.array(nr*[ncol*[0.0]])
    raw[:,0]=vv[1:]
    raw[:,1]=mass
    raw[:,2]=temp
    raw[:,3]=cye
    raw[:,4:]=massfr
    #-- format input structure
    #-- headers
    hd1='spherical\n'
    hd2=str(nr)+' 1 1 '+str(ncol)+' '+str(nelem)+'\n'
    #-- column labels
    lbls=[' rightvel','mass','temp','ye']
    hd3=lbls[0].rjust(1)
    for lbl in lbls[1:]: hd3+=lbl.rjust(12)
    #-- add element label to header
    for eld in elsd[:,0]: hd3+=elem_data.eldict_inv[eld].lower().rjust(12)
    #-- total header
    hd=hd1+hd2+hd3
    #-- incorporate el symbols in file name
    np.savetxt('input.str_x'+str(nr),raw,fmt='% .4e',
               delimiter=' ',header=hd)
