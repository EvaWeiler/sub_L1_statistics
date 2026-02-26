#!/usr/bin/env python
"""
Functions for index predictions.

Author: C. Moestl, R. Bailey, IWF Graz, Austria
twitter @chrisoutofspace, https://github.com/cmoestl
started April 2018, last update May 2019

Python 3.7
Packages not included in anaconda installation: sunpy, cdflib (https://github.com/MAVENSDC/cdflib)

Issues:
- ...

To-dos:
- ...

Future steps:
- ...

"""

import copy
from datetime import datetime
from dateutil import tz

import numpy as np
from matplotlib.dates import num2date, date2num
from numba import njit, jit
import astropy.time
import scipy


def calc_dst_temerin_li(time, btot, bx, by, bz, speed, speedx, density, dst1, dst2, dst3, version='2002n', linear_t_correction=False, minute_res=True):
    """Calculates Dst from solar wind input according to Temerin and Li 2002 method.
    Credits to Xinlin Li LASP Colorado and Mike Temerin.
    Calls _jit_calc_dst_temerin_li. All constants are defined in there.
    Note: vx has to be used with a positive sign throughout the calculation.

    Parameters
    ==========
    time : np.array
        Array containing time variables.
    btot : np.array
        Array containing Btot.
    bx : np.array
        Array containing Bx in coordinate system ?.
    by : np.array
        Array containing By in coordinate system ?.
    bz : np.array
        Array containing Bz in coordinate system ?.
    speed : np.array
        Array containing solar wind speed.
    speedx : np.array
        Array containing solar wind speed in x-direction.
    density : np.array
        Array containing solar wind density.
    version : str (default='2002')
        String determining which model version should be used.

    Returns
    =======
    dst_burton : np.array
        Array with calculated Dst values over timesteps time.
    """
    
    # Arrays
    dst1_=np.zeros(len(bz))
    dst2_=np.zeros(len(bz))
    dst3_=np.zeros(len(bz))
    dst_tl=np.zeros(len(bz))
    
    # Define initial values (needed for convergence, see Temerin and Li 2002 note)
    #dst1_[0:10]=-15
    #dst2_[0:10]=-13
    #dst3_[0:10]=-2

    if version == '2002':
        newparams = False
    else:
        newparams = True

    if version in ['2002', '2002n']:
        # julian_days = [sunpy.time.julian_day(num2date(x)) for x in time]
        julian_days = [astropy.time.Time(num2date(x), format='datetime', scale='utc').jd for x in time]
        return _jit_calc_dst_temerin_li_2002(time, btot, bx, by, bz, speed, speedx, density, dst1, dst2, dst3, dst_tl, julian_days, newparams=newparams)
    elif version == '2006':
        dst1_[0], dst2_[0], dst3_[0] = dst1, dst2, dst3# -10, -5, -10
        #dst1[0:2], dst2[0:2], dst3[0:2] = -10, -5, -10
        ds1995 = time - date2num(datetime(1995,1,1))
        ds2000 = time - date2num(datetime(2000,1,1))

        # YEARLY DRIFT CORRECTION TERM (NOT IN PAPER)
        if linear_t_correction:
            drift_corr = -0.014435865642103548 * ds2000 + 9.57670996872173
        else:
            drift_corr = 0.
        return _jit_calc_dst_temerin_li_2006(ds1995, ds2000, btot, bx, by, bz, speed, speedx, density, dst1_, dst2_, dst3_, minute_res=minute_res) + drift_corr


@njit(parallel=True)#"void(f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:])"
def _jit_calc_dst_temerin_li_2006(t1, t2, btot, bx, by, bz, speed, speedx, density, dst1, dst2, dst3,minute_res=True):
    """Fast(er) calculation of Dst using jit on Temerin-Li method."""

    fy = 2.*np.pi/365.24
    sun1 = np.sin(np.pi * 10.27 / 180.)
    sun2 = np.sin(np.pi * 10.27 / 180.) * np.cos(np.pi * 23.5 / 180.)
    alpha = 0.0449

    # SOLAR WIND VALUES
    # -----------------
    speedx = np.abs(speedx)
    bt = np.sqrt(by**2 + bz**2)
    bt[bt < 0.0001] = 1e-4      # Correction from dst.pro, escaping zero-division error
    bp = np.sqrt(by**2 + bx**2)
    btot = np.sqrt(bx**2. + by**2. + bz**2.)

    theta = -(np.arccos(-bz/bt) - np.pi) / 2.
    ang = np.arctan2(bx, by)

    exx = speedx * bt**0.993 * np.sin(theta)**7.29
    exx2 = density**0.493 * speedx**2.955 * bt**1.105 * np.sin(theta)**5.24 # paper
    exx3 = density**0.397 * speedx**0.576 * bt**1.413 * np.sin(theta)**8.56 # paper
    #exx2 = speedx * bt**1.105 * np.sin(theta)**5.24 # code
    #exx3 = speedx * bt**1.413 * np.sin(theta)**8.56 # code

    # TIME VALUES
    # -----------
    if minute_res == True:
        index = 400
    else: 
        index = 10
    
    dh = 0.0435 * np.cos(fy*t1 + 0.1680) - 0.0208 * np.sin(2*np.pi*t1 - 1.589)
    itest = index
    it1 = itest - np.where(t1 > (t1[itest] - 0.0486))[0][0]
    it2 = itest - np.where(t1 > (t1[itest] - 0.181))[0][0]
    it3 = itest - np.where(t1 > (t1[itest] - 0.271))[0][0]
    it4 = itest - np.where(t1 > (t1[itest] - 0.0625))[0][0]
    it5 = itest - np.where(t1 > (t1[itest] - 0.104))[0][0]
    it6 = itest - np.where(t1 > (t1[itest] - 0.0278))[0][0]
    it7 = itest - np.where(t1 > (t1[itest] - 0.139))[0][0]
    idst1t1 = itest - np.where(t1 > (t1[itest] - 0.132))[0][0]
    idst2t1 = itest - np.where(t1 > (t1[itest] - 0.0903))[0][0]
    idst1t2 = itest - np.where(t1 > (t1[itest] - 0.264))[0][0]
    #print('time indices: ', t1[it1], t1[it2], t1[it3], t1[it4], t1[it5], t1[it6], t1[it7])
    #print('time idst1t1, idst1t2: ', t1[idst1t1], t1[idst1t2])
    #print('indices: ', it1, it2, it3, it4, it5, it6, it7)

    # FUNCTION TERMS
    # --------------
    tt = t1*fy
    cosphi = sun2 * np.sin(tt + alpha) * np.sin(2.*np.pi*t1 - tt - 1.632) + \
                np.cos(tt + alpha) * (0.39 + sun1*np.cos(2*np.pi*t1 - tt - 1.632))
    cosphi5 = sun2 * np.sin(tt + alpha) * np.sin(2.*np.pi*t1 - tt + 0.27) + \
                np.cos(tt + alpha) * (0.39 + sun1*np.cos(2*np.pi*t1 - tt + 0.27))
    cosphi6 = sun2 * np.sin(tt + alpha) * np.sin(2.*np.pi*t1 - tt - 0.21) + \
                np.cos(tt + alpha) * (0.39 + sun1*np.cos(2*np.pi*t1 - tt - 0.21))
    cosphi7 = sun2 * np.sin(tt + alpha) * np.sin(2.*np.pi*t1 - tt - 0.79) + \
                np.cos(tt + alpha) * (0.39 + sun1*np.cos(2*np.pi*t1 - tt - 0.79))
    cosphi8 = sun2 * np.sin(tt + alpha) * np.sin(2.*np.pi*t1 - tt - 2.81) + \
                np.cos(tt + alpha) * (0.39 + sun1*np.cos(2*np.pi*t1 - tt - 2.81))

    sin_phi_factor = 0.95097
    tst3 = ( np.sqrt(1. - cosphi**2.) / sin_phi_factor )**-0.13
    tst4 = ( np.sqrt(1. - cosphi**2.) / sin_phi_factor )**6.54
    tst5 = ( np.sqrt(1. - cosphi5**2.) / sin_phi_factor )**5.13
    tst6 = ( np.sqrt(1. - cosphi6**2.) / sin_phi_factor )**-2.44
    tst7 = ( np.sqrt(1. - cosphi7**2.) / sin_phi_factor )**2.84
    tst8 = ( np.sqrt(1. - cosphi8**2.) / sin_phi_factor )**2.49

    fe1 = -1.703e-6 * (1. + erf(-0.09*bp * np.cos(ang - 0.015) * dh)) * \
                tst3 * ((exx - 1231.2/tst4 + np.abs(exx - 1231.2/tst4)) + \
                (exx - 3942./tst4 + np.abs(exx - 3942./tst4))) * speedx**1.307 * density**0.548
    fe2 = -5.172e-8 * exx2 * (1. + erf(0.418*bp * np.cos(ang - 0.015) * dh) )    # paper
    fe3 =  -0.0412 * exx3 * (1. + erf(1.721*bp * np.cos(ang - 0.015) * dh) )    # paper
    #fe2 = -5.172e-8 * exx2 * (1. + erf(0.418*bp * np.cos(ang - 0.015) * dh) ) * speedx**1.955 * density**0.493   
    #fe3 =  -0.0412 * exx3 * (1. + erf(1.721*bp * np.cos(ang - 0.015) * dh) ) * speedx**-0.424 * density**0.397  

    df2 = 1440. * tst7 * fe2/(-fe2 + 922.1)
    df3 = 272.9 * tst8 * fe3/(-fe3 + 60.5)

    # PRESSURE TERM
    # -------------
    pressureterm = ( 0.330*btot**2 * (1. + 0.100*density) + \
                (1.621e-4 * tst6 * speed**2 + 18.70)*density )**0.5

    # DIRECT BZ TERM
    # --------------
    directbzterm = 0.574 * tst5 * bz

    # OFFSET TERM
    # -----------
    offsetterm = 19.35 + 0.158*np.sin(fy*t2 - 0.94) + 0.01265*t2 - 2.224e-11*t2**2.
    dt = (t1[1] - t1[0])
 
    # INITIAL DST LOOP
    # ----------------
    for i in range(0,index):
        #if minute_res:
        
        #print(dt)
        # Code
        dst1[i+1] = dst1[i] + (0.005041 * (-dst1[i])**2.017 + fe1[i]) * dt
        dst2[i+1] = dst2[i] + (0.00955 * (-dst2[i])**2.269 + df2[i] * (1. + 0.01482*dst1[i] / (1. - 0.01482*dst1[i]) )) * dt
        dst3[i+1] = dst3[i] + (-5.10*dst3[i] + df3[i]) * dt

    # MAIN DST LOOP
    # -------------
    for i in range(index,len(bz)-1):

        # DST TERMS
        # ---------
        bzt1 = bz[i-it1]
        bzt2 = bz[i-it2]
        bzt3 = bz[i-it3]
        dst1t1 = dst1[i-idst1t1]
        dst2t1 = dst2[i-idst2t1]
        dst1[i+1] = dst1[i] + (5.041e-3 * (-dst1[i])**2.017 * \
                    (1. + erf(-0.010*bz[i])) + fe1[i] * \
                    (1. + erf(-0.0094*bzt1 - 0.0118*bzt2 + 0.0138*bzt3)) * \
                    np.exp(0.00313*dst1t1 + 0.01872*dst2t1)) * dt

        bzt4 = bz[i-it4]
        bzt5 = bz[i-it5]
        dst1t2 = dst1[i-idst1t2]
        dst2[i+1] = dst2[i] + (0.00955 * (-dst2[i])**2.017 * \
                    (1. + erf(-0.014*bz[i])) + df2[i] * \
                    (1 + erf(-0.0656*bzt4 + 0.0627*bzt5)) * \
                    np.exp(0.01482*dst1t2)) * dt

        bzt6 = bz[i-it6]
        bzt7 = bz[i-it7]
        dst3[i+1] = dst3[i] + (5.10 * (-dst3[i])**0.952 * \
                    (1. + erf(-0.027*bz[i])) + df3[i] * \
                    (1. + erf(-0.0471*bzt6 + 0.0184*bzt7))) * dt

    # ANNUAL VARIATIONS
    # -----------------
    dst1_ = dst1 * (1. + 0.0807*np.sin(t1*fy + 1.886))
    dst2_ = dst2 * (1. + 0.0251*np.sin(t1*fy + 3.18))
    dst3_ = dst3 * (1. + 0.0901*np.sin(t1*fy + 5.31)) * (1.-0.00007*dst1_)
    #directbzterm_ = directbzterm * (1. + 0.293 * np.sin(t1[i]*fy + 3.19)) * (1. + 0.0034*dst1_) # paper
    directbzterm_ = directbzterm * (1. + 0.293 * np.sin(t1[i]*fy + 3.19))   # code
    pressureterm_ = pressureterm * (1. + 0.0986*np.sin(t1[i]*fy - 1.383)) * (1. + 0.00184*dst1_)

    # FINAL DST
    # ---------
    dst_tl = dst1_ + dst2_ + dst3_ + pressureterm_ + directbzterm_ + offsetterm

    return dst_tl

@njit
def erf(x):
    # adjusted from https://stackoverflow.com/questions/457408/is-there-an-easily-available-implementation-of-erf-for-python
    # save the sign of x
    sign = np.sign(x)
    x = np.abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1) * t * np.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)


def calc_dst_burton(time, bz, speed, density):
    """Calculates Dst from solar wind input according to Burton et al. 1975 method.

    Parameters
    ==========
    time : np.array
        Array containing time variables.
    bz : np.array
        Array containing Bz in coordinate system ?.
    speed : np.array
        Array containing solar wind speed.
    density : np.array
        Array containing Bz in coordinate system ?.

    Returns
    =======
    dst_burton : np.array
        Array with calculated values over timesteps time.
    """

    protonmass=1.6726219*1e-27  #kg
    bzneg = copy.deepcopy(bz)
    bzneg[bz > 0] = 0
    pdyn = density*1e6*protonmass*(speed*1e3)**2*1e9  #in nanoPascal
    Ey = speed*abs(bzneg)*1e-3 #now Ey is in mV/m

    dst_burton = np.zeros(len(bz))
    Ec=0.5  
    a=3.6*1e-5
    b=0.2*100 #*100 due to different dynamic pressure einheit in Burton
    c=20  
    d=-1.5/1000 
    lrc=0
    for i in range(len(bz)-1):
        if Ey[i] > Ec:
            F = d*(Ey[i]-Ec) 
        else: F=0
        #Burton 1975 p4208: Dst=Dst0+bP^1/2-c
        # Ring current Dst
        deltat_sec = (time[i+1]-time[i])*86400  #deltat must be in seconds
        rc = lrc + (F-a*lrc)*deltat_sec
        # Dst of ring current and magnetopause currents 
        dst_burton[i+1] = rc + b*np.sqrt(pdyn[i+1]) - c
        lrc = rc

    return dst_burton


def calc_dst_obrien(time, bz, speed, density):
    """Calculates Dst from solar wind input according to OBrien and McPherron 2000 method.

    Parameters
    ==========
    time : np.array
        Array containing time variables.
    bz : np.array
        Array containing Bz in coordinate system ?.
    speed : np.array
        Array containing solar wind speed.
    density : np.array
        Array containing Bz in coordinate system ?.

    Returns
    =======
    dst_burton : np.array
        Array with calculated values over timesteps time.
    """

    protonmass=1.6726219*1e-27  #kg
    bzneg = copy.deepcopy(bz)
    bzneg[bz > 0] = 0
    pdyn = density*1e6*protonmass*(speed*1e3)**2*1e9  #in nanoPascal
    Ey = speed*abs(bzneg)*1e-3; #now Ey is in mV/m

    Ec=0.49
    b=7.26  
    c=11  #nT
    lrc=0
    dst_obrien = np.zeros(len(bz))
    for i in range(len(bz)-1):
        if Ey[i] > Ec:            #Ey in mV m
            Q = -4.4 * (Ey[i]-Ec) 
        else: Q=0
        tau = 2.4 * np.exp(9.74/(4.69 + Ey[i])) #tau in hours
        # Ring current Dst
        deltat_hours=(time[i+1]-time[i])*24 # time should be in hours
        rc = ((Q - lrc/tau))*deltat_hours + lrc
        # Dst of ring current and magnetopause currents 
        dst_obrien[i+1] = rc + b*np.sqrt(pdyn[i+1])-c; 
        lrc = rc

    return dst_obrien