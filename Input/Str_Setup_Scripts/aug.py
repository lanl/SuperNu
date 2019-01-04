#!/bin/python

#-----------------------------------------------------
# Format Rosswog input structures and add initial gas
# temperature.
#-----------------------------------------------------

#-- modules
import numpy as np
import argparse

#-- parameters
parser = argparse.ArgumentParser(description='Format Rosswog'\
                                 ' input structures and add'\
                                 ' initial gas temperature.')
parser.add_argument('f',type=str,help='Select raw data.')
parser.add_argument('--igeom',type=int,help='Select geometry.')
parser.add_argument('--t',type=float,help='Select initial'\
                    ' expansion time [s].')
parser.add_argument('--el',type=str,help='Select element'\
                    ' (element symbol).')
parser.add_argument('--model',type=str,help='Select model'\
                    ' (replaces element symbol in file name).')
parser.add_argument('--opac',type=float,help='Grey dyn. ej. opacity [cm^2/g]')
parser.add_argument('--wind_el',type=str,help='Add wind element')
parser.add_argument('--wind_mass',type=float,help='Wind mass [solar mass]')
parser.add_argument('--wind_vmax',type=float,help='Max wind velocity [c]')
parser.add_argument('--wind_opac',type=float,help='Grey wind opacity [cm^2/g]')
parser.add_argument('--add_dyn_frac',action='store_true',help='Add column'\
                    ' of dynamical ejecta fraction to structure.')
#-- parameters for rescaling dynamical ejecta speed, mass
parser.add_argument('--dyn_mass_rescale',type=float,help='Select dynamical ejecta'\
                    ' mass rescaling factor.')
parser.add_argument('--dyn_vmax_rescale',type=float,help='Select dynamical ejecta'\
                    ' velocity rescaling factor.')
#-- parse parameters
args = parser.parse_args()

#-- constants
c = 2.997e10 # light speed [cm/s]
msol = 1.989e33 # solar mass [g]

#-- defaults
if(args.t==None): t=10000.
else: t=args.t
if(args.el==None): el='Sm'
else: el=args.el
if(args.igeom==None): igeom=11
else: igeom=args.igeom
#-- spherical analytic wind parameters
if(args.wind_mass==None): mw=5e-3*msol
else: mw=args.wind_mass*msol
if(args.wind_vmax==None): vmaxw=.16*c
else: vmaxw=args.wind_vmax*c
#-- there can be no wind
wind_el = args.wind_el
#-- dynamical ejecta rescaling factors
if(args.dyn_mass_rescale==None): dyn_mass_rescale = 1.0
else: dyn_mass_rescale = args.dyn_mass_rescale
if(args.dyn_vmax_rescale==None): dyn_vmax_rescale = 1.0
else: dyn_vmax_rescale = args.dyn_vmax_rescale

#-- raw input structure
raw=np.loadtxt(args.f,skiprows=3)

#-- set the mass index
if(igeom==11): imass = 2
elif(igeom==2): imass = 4

#-- rescale the dynamical ejecta velocity and mass
raw[:,:imass] *= dyn_vmax_rescale
raw[:,imass] *= dyn_mass_rescale
#-- print preliminary ejecta values
print("Total Mass [Msol]: ",np.sum(raw[:,imass])/msol)
print("Median R-Velocity [c]: ",0.5*np.amax(raw[:,0])/c)

if(igeom==11):
    #-- number of cells
    nx=np.size(raw,0)
    ny=1
    #-- total number of cells
    ncell=nx
    #-- volume [cc]
    vol=4*np.pi*(raw[:,1]**3-raw[:,0]**3)*t**3
elif(igeom==2):

    #-- extension flag
    extflag = 0

    #
    #-- x grid
    xleft = np.unique(raw[:,0])
    xright = np.unique(raw[:,1])
    cx = .5*(xleft+xright)
    #-- keep original grid
    cxold = cx.copy()
    #-- number of x-cells
    nx=np.size(cx)
    #-- check if wind is extended off x-range
    if(vmaxw>xright[nx-1]):
        #-- flag extension
        extflag = 1
        #-- extend space in x-direction
        dx = xright[nx-1]-xleft[nx-1]
        #-- number of additional cells
        nxadd = int((vmaxw-xright[nx-1])/dx)+1
        #-- reset grid values
        xleft = np.append(xleft,np.linspace(xleft[nx-1]+dx,
                                            xleft[nx-1]+nxadd*dx,nxadd))
        xright = np.append(xright,np.linspace(xright[nx-1]+dx,
                                              xright[nx-1]+nxadd*dx,nxadd))
        cx = .5*(xleft+xright)
        nx += nxadd

    #
    #-- y grid
    yleft = np.unique(raw[:,2])
    yright = np.unique(raw[:,3])
    cy = .5*(yleft+yright)
    #-- keep original grid
    cyold = cy.copy()
    #-- number of y-cells
    ny=np.size(cy)
    #-- check if wind is extended off positive y-range
    if(vmaxw>yright[ny-1]):
        #-- flag extension
        extflag = 1
        #-- extend space in y-direction
        dy = yright[ny-1]-yleft[ny-1]
        #-- number of additional cells
        nyadd = int((vmaxw-yright[ny-1])/dy)+1
        #-- reset grid values
        yleft = np.append(yleft,np.linspace(yleft[ny-1]+dy,
                                            yleft[ny-1]+nyadd*dy,nyadd))
        yright = np.append(yright,np.linspace(yright[ny-1]+dy,
                                              yright[ny-1]+nyadd*dy,nyadd))
        cy = .5*(yleft+yright)
        ny += nyadd
    #-- check if wind is extended off negative y-range
    if(vmaxw>abs(yleft[0])):
        #-- flag extension
        extflag = 1
        #-- extend space in y-direction
        dy = yright[0]-yleft[0]
        #-- number of additional cells
        nyadd = int((vmaxw-abs(yleft[0]))/dy)+1
        #-- reset grid values
        yleft = np.append(np.linspace(yleft[0]-nyadd*dy,
                                      yleft[0]-dy,nyadd),yleft)
        yright = np.append(np.linspace(yright[0]-nyadd*dy,
                                       yright[0]-dy,nyadd),yright)
        cy = .5*(yleft+yright)
        ny += nyadd

    #
    #-- if wind is extended, extend raw data
    if(extflag==1):
        #-- initialized extended structure
        ncolr = np.size(raw,1)
        if(ncolr!=5): print('ncolr != 5, ncolr = ',ncolr)
        rawext = np.array(nx*ny*[ncolr*[0.]])
        for jj in range(ny):
            #-- x-edge
            rawext[jj*nx:(jj+1)*nx,0]=xleft
            rawext[jj*nx:(jj+1)*nx,1]=xright
            #-- y-edge
            rawext[jj*nx:(jj+1)*nx,2]=yleft[jj]
            rawext[jj*nx:(jj+1)*nx,3]=yright[jj]
        #-- set masses
        nxold = np.size(cxold)
        nyold = np.size(cyold)
        for j in range(nyold):
            for i in range(nxold):
                #-- get indices in extended grids
                jj = np.argmin(abs(cy-cyold[j]))
                ii = np.argmin(abs(cx-cxold[i]))
                rawext[ii+jj*nx,4]=raw[i+j*nxold,4]
        #-- reset raw data array
        raw = rawext

    #
    #-- total number of cells
    ncell=nx*ny
    #-- volume [cc]
    vol=np.kron(yright-yleft,.5*(xright**2-xleft**2))
    vol*=2*np.pi*t**3

if(wind_el!=None):
    #-- spherical density function
    def frhow(v,vm):
        "functional from of wind density"
        rhof=(1.0-(v/vm)**2)**3
        np.putmask(rhof,v>vm,0.0)
        return rhof
    if(igeom==11): print('no 1d wind addition')
    elif(igeom==2):
        #-- sanity check
        if(vmaxw>cx.max()): print('vmaxw>cx, augmented.')
        if(vmaxw>cy.max() or vmaxw>abs(cy.min())): print('vmaxw>cy, augmented.')
        vw=np.sqrt(cx[None,:]**2+cy[:,None]**2).ravel()
        rhow=frhow(vw,vmaxw)
        #-- mass function form
        massw=rhow*vol
        #-- re-normalize to correct total mass
        hlp=mw/np.sum(massw)
        massw*=hlp

#-- density [g/cc]
massd=raw[:,imass]
mass=massd.copy()
if(wind_el!=None): mass+=massw
rho=mass/vol
#-- temperature [K]
t0=10. #[s]
temp0=5.7e7 #[K]
rho0=np.max(rho)*(t/t0)**3
temp=temp0*(rho/rho0)**(1./3)
#-- electron fraction
ye0d = 0.39
ye0w = 0.42
ye0 = ye0d*massd
if(wind_el!=None): ye0+=massw*ye0w
np.putmask(ye0,mass==0.0,ye0d)
massh=mass.copy(); massh[massh==0.0]=1.0
ye0/=massh
#-- grey opacity
opac = np.array([])
if(args.opac!=None):
    opac = args.opac*massd
    if(args.wind_opac!=None):
        opac += massw*args.wind_opac
    opac /= massh

#-- element dictionary
eldict = {'Cr':24, 'Fe':26, 'Se':34, 'Br':35, 'Zr':40,
          'Pd':46, 'Te':52, 'Ce':58, 'Nd':60,
          'Sm':62, 'U':92}
eldict_inv = {elv: elk for elk, elv in eldict.items()}

#-- if element strings are files, load them
if ('.dat' in el):
    elsd = np.loadtxt(el)
    nelemd = np.size(elsd,0)
else:
    elsd = np.array([[eldict[el],1.0]])
    nelemd = 1

#-- initialize element list
els = elsd.copy()

#-- total input structure
ncol=imass + 3 + nelemd
if(args.add_dyn_frac): ncol += 1
if(args.opac!=None): ncol += 1
if(wind_el!=None):
    if('.dat' in wind_el):
        elsw = np.loadtxt(wind_el)
        nelemw = np.size(elsw,0)
    else:
        elsw = np.array([[eldict[wind_el],1.0]])
        nelemw = 1
    #-- add to column number if wind has new elements
    ncommon = np.size(list(set(elsw[:,0]).intersection(set(elsd[:,0]))))
    naddw = nelemw - ncommon
    ncol += naddw
    #-- number of elements to add to wind element list
    naddd = nelemd - ncommon
    #-- augment and sort dynamical ejecta element list
    if(naddw > 0):
        eladd_dyn = np.array(list((set(elsw[:,0])^set(elsd[:,0])) & set(elsw[:,0])))
        eladdd = np.array(np.size(eladd_dyn)*[2*[0.0]])
        eladdd[:,0] = eladd_dyn
        elsd = np.concatenate((elsd,eladdd))
        elsd = elsd[elsd[:,0].argsort()]
    #-- augment and sort wind element list
    if(naddd > 0):
        eladd_wind = np.array(list((set(elsw[:,0])^set(elsd[:,0])) & set(elsd[:,0])))
        eladdw = np.array(np.size(eladd_wind)*[2*[0.0]])
        eladdw[:,0] = eladd_wind
        elsw = np.concatenate((elsw,eladdw))
        elsw = elsw[elsw[:,0].argsort()]

#-- generate input structure
spnstr=np.array(ncell*[ncol*[0.]])
spnstr[:,:imass]=raw[:,:imass]
spnstr[:,imass]=mass
spnstr[:,imass+1]=temp
icol = 3
if(args.add_dyn_frac): icol+=1
if(args.opac!=None):
    icol+=1
    spnstr[:,imass+icol-2] = opac
if(wind_el==None):
    #-- only dyn ejecta
    spnstr[:,imass+2]=ye0d
    spnstr[:,imass+icol:]=np.repeat(elsd[:,1:2].T,ncell,axis=0)
    if(args.add_dyn_frac): spnstr[:,imass+icol-1] = 1.0
else:
    #-- dyn ejecta + wind
    fracd = massd/massh
    fracw = massw/massh
    spnstr[:,imass+2]=ye0
    chid = np.outer(fracd,elsd[:,1])
    chiw = np.outer(fracw,elsw[:,1])
    spnstr[:,imass+icol:] = chid + chiw
    if(args.add_dyn_frac): spnstr[:,imass+icol-1] = fracd
    np.putmask(spnstr[:,imass+icol],mass==0.0,1.0)

#-- re-normalize mass fractions
for icell in range(ncell): spnstr[icell,imass+icol:] /= np.sum(spnstr[icell,imass+icol:])

if(igeom==2):
    rout = min(xright[nx-1],yright[ny-1])**2
    for icell in range(ncell):
        if(spnstr[icell,imass] > 0.0):
            routh = spnstr[icell,1]**2 + max(abs(spnstr[icell,2]),abs(spnstr[icell,3]))**2
            if(routh > rout): rout = routh
    #-- final
    rout = np.sqrt(rout)
    print("rout = ", rout, "rout / c = ", rout/c)

#-- format input structure
#-- headers
if(args.model==None):
    if(igeom==11): hd1='spherical\n'
    elif(igeom==2): hd1='cylindrical\n'
else: hd1=args.model+'\n'
nelem = ncol - imass - icol
hd2=str(nx)+' '+str(ny)+' 1 '+str(ncol)+' '+str(nelem)+'\n'
#-- column labels
if(igeom==11): lbls=['   x_left','x_right','mass','temp','ye']
elif(igeom==2): lbls=['   x_left','x_right','y_left','y_right','mass','temp','ye']
if(args.opac!=None): lbls += ['cap']
if(args.add_dyn_frac): lbls += ['dyn_fr']
hd3=lbls[0]
for lbl in lbls[1:]: hd3+=lbl.rjust(12)
#-- add element label to header
for eld in elsd[:,0]: hd3 += eldict_inv[eld].lower().rjust(12)
#-- total header
hd=hd1+hd2+hd3
#-- incorporate el symbol in file name
lbl='' #el.strip()
if(wind_el!=None and '.dat' not in wind_el): lbl+=wind_el.strip()
#-- add model name to file name (e.g. A1d, B1d, C1d, D1d)
if(args.model!=None): lbl=args.model+lbl
#-- save
flabel=lbl+'_x'+str(nx)
if(igeom==2): flabel+='y'+str(ny)
np.savetxt('input.rho_'+flabel,rho,fmt='% .4e')
np.savetxt('input.str_'+flabel,spnstr,fmt='% .4e',
           delimiter=' ',header=hd)
