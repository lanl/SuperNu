#!/usr/bin/env python

import numpy as np

#-- constants
pc_c = 2.998e10 # [cm/s]

def calc_tau(time, nx, grd_vel, ng, flx_wl, opac):
    '''Calculate comoving optical depths at sime time.'''

    #-- optical depth array
    taus = np.zeros((ng,nx))
    dift = np.zeros((ng,nx))

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
            tau += time * min(dvx,dvg)*opac[ig,ix]
            dif += opac[ig,ix] * (time * min(dvx,dvg))**2 / pc_c
            #-- update ray parameter
            if (dvx < dvg):
                # set new ray state
                v0 = grd_vel[ix]
                wl0 /= (1.0 + dvx/pc_c)
                # set the tau for this cell
                taus[jg,ix] = tau
                dift[jg,ix] = dif
                # decrement cell index
                ix -= 1
            else:
                # set new ray state
                wl0 = flx_wl[ig]
                v0 -= dvg
                # decrement group index
                ig -= 1

    #-- return taus and dift
    return taus, dift
