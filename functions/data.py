# Standard
import os
import sys
import copy
from datetime import datetime, timedelta, timezone
from dateutil.relativedelta import relativedelta
from dateutil import tz
import gzip
import h5py
import logging
import numpy as np
import pdb
import pickle
import re
import scipy
import scipy.io
import shutil
import subprocess
from matplotlib.dates import date2num, num2date
from glob import iglob
import json
import urllib
from sunpy.time import parse_time
import pandas as pd

# External
import astropy.time

# Local
from .predict import calc_dst_temerin_li

logger = logging.getLogger(__name__)



# =======================================================================================
# ----------------------------- I. CLASS SatDat() ---------------------------------------
# =======================================================================================

class SatData():
    """Data object containing satellite data.

    Init Parameters
    ===============
    --> SatData(input_dict, source=None, header=None)
    input_dict : dict(key: dataarray)
        Dict containing the input data in the form of key: data (in array or list)
        Example: {'time': timearray, 'bx': bxarray}. The available keys for data input
        can be accessed in SatData.default_keys.
    header : dict(headerkey: value)
        Dict containing metadata on the data array provided. Useful data headers are
        provided in SatData.empty_header but this can be expanded as needed.
    source : str
        Provide quick-access name of satellite/data type for source.

    Attributes
    ==========
    .data : np.ndarray
        Array containing measurements/indices. Best accessed using SatData[key].
    .position : np.ndarray
        Array containing position data for satellite.
    .h : dict
        Dict of metadata as defined by input header.
    .state : np.array (dtype=object)
        Array of None, str if defining state of data (e.g. 'quiet', 'cme').
    .vars : list
        List of variables stored in SatData.data.
    .source : str
        Data source name.

    Methods
    =======
    .convert_GSE_to_GSM()
        Coordinate conversion.
    .convert_RTN_to_GSE()
        Coordinate conversion.
    .cut(starttime=None, endtime=None)
        Cuts data to within timerange and returns.
    .get_position(timestamp)
        Returns position of spacecraft at time.
    .get_newell_coupling()
        Calculates Newell coupling indice for data.
    .interp_nans(keys=None)
        Linearly interpolates over nans in data.
    .interp_to_time()
        Linearly interpolates over nans.
    .load_position_data(position_data_file)
        Loads position data from file.
    .make_aurora_power_prediction()
        Calculates aurora power.
    .make_dst_prediction()
        Makes prediction of Dst from data.
    .make_kp_prediction()
        Prediction of kp.
    .make_hourly_data()
        Takes minute resolution data and interpolates to hourly data points.
    .shift_time_to_L1()
        Shifts time to L1 from satellite ahead in sw rotation.

    Examples
    ========
    """

    default_keys = ['time', 'time_utc',
                    'speed', 'speedx', 'speedy', 'speedz', 'density', 'temp', 'pdyn',
                    'bx', 'by', 'bz', 'btot',
                    'br', 'bt', 'bn',
                    'dst', 'kp', 'aurora', 'ec', 'ae', 'f10.7', 'symh', 'time_shifted', 'expansion_delta',
                    'btot_err', 'br_err', 'bx_err', 'bt_err', 'by_err', 'bn_err', 'bz_err', 'speed_err', 'dst_err', 'dst_err_min', 'dst_err_max', 'symh_err', 'symh_err_min', 'symh_err_max']

    empty_header = {'DataSource': '',
                    'SourceURL' : '',
                    'SamplingRate': None,
                    'ReferenceFrame': '',
                    'FileVersion': {},
                    'Instruments': [],
                    'RemovedTimes': [],
                    'PlasmaDataIntegrity': 10
                    }

    def __init__(self, input_dict, source=None, header=None):
        """Create new instance of class."""

        # Check input data
        for k in input_dict.keys():
            if not k in SatData.default_keys: 
                raise NotImplementedError("Key {} not implemented in SatData class!".format(k))
        if 'time' not in input_dict.keys():
            raise Exception("Time variable is required for SatData object!")
        dt = [x for x in SatData.default_keys if x in input_dict.keys()]
        
        if len(input_dict['time']) == 0:
            logger.warning("SatData.__init__: Inititating empty array! Is the data missing?")
        # Create data array attribute
        data = [input_dict[x] if x in dt else np.zeros(len(input_dict['time'])) for x in SatData.default_keys]
        self.data = np.asarray(data)
        # Create array for state classifiers (currently empty)
        self.state = np.array([None]*len(self.data[0]), dtype='object')
        # Add new attributes to the created instance
        self.source = source
        if header == None:               # Inititalise empty header
            self.h = copy.deepcopy(SatData.empty_header)
        else:
            self.h = header
        self.pos = None
        self.vars = dt
        self.vars.remove('time')
        
    # -----------------------------------------------------------------------------------
    # Internal methods
    # -----------------------------------------------------------------------------------

    def __getitem__(self, var):
        if isinstance(var, str):
            if var in self.vars+['time']:
                return self.data[SatData.default_keys.index(var)]
            else:
                raise Exception("SatData object does not contain data under the key '{}'!".format(var))
        return self.data[:,var]


    def __setitem__(self, var, value):
        if isinstance(var, str):
            if var in self.vars:
                self.data[SatData.default_keys.index(var)] = value #added [0]
            elif var in SatData.default_keys and var not in self.vars:
                self.data[SatData.default_keys.index(var)] = value #added [0]
                self.vars.append(var)
            else:
                raise Exception("SatData object does not contain the key '{}'!".format(var))
        else:
            raise ValueError("Cannot interpret {} as index for __setitem__!".format(var))


    def __len__(self):
        return len(self.data[0])


    def __str__(self):
        """Return string describing object."""

        ostr = "Length of data:\t\t{}\n".format(len(self))
        ostr += "Keys in data:\t\t{}\n".format(self.vars)
        ostr += "First data point:\t{}\n".format(num2date(self['time'][0]))
        ostr += "Last data point:\t{}\n".format(num2date(self['time'][-1]))
        ostr += "\n"
        ostr += "Header information:\n"
        for j in self.h:
            if self.h[j] != None: ostr += "    {:>25}:\t{}\n".format(j, self.h[j])
        ostr += "\n"
        ostr += "Variable statistics:\n"
        ostr += "{:>12}{:>12}{:>12}\n".format('VAR', 'MEAN', 'STD')
        for k in self.vars:
            ostr += "{:>12}{:>12.2f}{:>12.2f}\n".format(k, np.nanmean(self[k]), np.nanstd(self[k]))

        return ostr
    

    # =======================================================================================
    # -------------------------------- II. DATA HANDLING ------------------------------------
    # =======================================================================================


    def cut(self, starttime=None, endtime=None):
        """Cuts array down to range defined by starttime and endtime. One limit
        can be provided or both.

        Parameters
        ==========
        starttime : datetime.datetime object
            Start time (>=) of new array.
        endtime : datetime.datetime object
            End time (<) of new array.

        Returns
        =======
        self : obj within new time range
        """

        if starttime != None and endtime == None:
            new_inds = np.where(self.data[0] >= date2num(starttime))[0]
            self.data = self.data[:,new_inds]
            if self.pos != None:
                self.pos.positions = self.pos.positions[:,new_inds]
        elif starttime == None and endtime != None:
            new_inds = np.where(self.data[0] < date2num(endtime))[0]
            self.data = self.data[:,new_inds]
            if self.pos != None:
                self.pos.positions = self.pos.positions[:,new_inds]
        elif starttime != None and endtime != None:
            new_inds = np.where((self.data[0] >= date2num(starttime)) & (self.data[0] < date2num(endtime)))[0]
            self.data = self.data[:,new_inds]
            if self.pos != None:
                self.pos.positions = self.pos.positions[:,new_inds]
        return self
    
    
    def interp_nans(self, keys=None, return_masked_array=False):
        """Linearly interpolates over nans in array.

        Parameters
        ==========
        keys : list (default=None)
            Provide list of keys (str) to be interpolated over, otherwise all.
        """

        logger.info("interp_nans: Interpolating nans in {} data".format(self.source))
        if keys == None:
            keys = self.vars
        if return_masked_array:
            orig_data = np.ma.array(copy.deepcopy(self.data))
        for k in keys:
            inds = np.isnan(self[k])
            if len(inds) == 0:
                return self
            if len(inds.nonzero()[0]) == len(self[k]):
                logger.warning("Data column {} in {} data is full of nans. Dropping this column.".format(k, self.source))
                self.vars.remove(k)
                self[k] = 0.
            else:
                self[k][inds] = np.interp(inds.nonzero()[0], (~inds).nonzero()[0], self[k][~inds])
        if not return_masked_array:
            return self
        else:
            masked_array = np.ma.masked_where(np.isnan(orig_data), self.data)
            satdata_masked = copy.deepcopy(self)
            satdata_masked.data = masked_array
            return self, satdata_masked


    def interp_to_time(self, tarray, keys=None):
        """Linearly interpolates over nans in array.

        Parameters
        ==========
        tarray : np.ndarray
            Array containing new timesteps in number format.
        keys : list (default=None)
            Provide list of keys (str) to be interpolated over, otherwise all.
        """

        if keys == None:
            keys = self.vars

        # Create new time array
        data_dict = {'time': tarray}
        for k in keys:
            na = np.interp(tarray, self['time'], self[k])
            data_dict[k] = na

        # Create new data opject:
        newData = SatData(data_dict, header=copy.deepcopy(self.h), source=copy.deepcopy(self.source))
        newData.h['SamplingRate'] = tarray[1] - tarray[0]
        # Interpolate position data:
        if self.pos != None:
            newPos = self.pos.interp_to_time(self['time'], tarray)
            newData.pos = newPos

        return newData
    
    
    
    def make_hourly_data(self):
        """Takes data with minute resolution and interpolates to hour.
        Uses .interp_to_time(). See that function for more usability options.

        Parameters
        ==========
        None

        Returns
        =======
        Data_h : new SatData obj
            New array with hourly interpolated data. Header is copied from original.
        """

        # Round to nearest hour
        stime = self['time'][0] - self['time'][0]%(1./24.)
        # Roundabout way to get time_h ensures timings with full hours:
        nhours = (num2date(self['time'][-1])-num2date(stime)).total_seconds()/60./60.
        # Create new time array
        time_h = np.array(stime + np.arange(1, nhours)*(1./24.))
        Data_h = self.interp_to_time(time_h)

        return Data_h


    
    
    def shift_wind_to_L1(self, L1Pos=[]):
        """Corrects for differences in B and density values due to solar wind
        expansion at different radii.
        
        Exponents taken from Kivelson and Russell, Introduction to Space Physics (Ch. 4.3.2).
        Magnetic field components are scaled according to values in
        Hanneson et al. (2020) JGR Space Physics paper.
        https://doi.org/10.1029/2019JA027139

        Parameters
        ==========
        None

        Returns
        =======
        self
        """

        dttime = [num2date(t).replace(tzinfo=None) for t in self['time']]
        if len(L1Pos) == 0:
            L1Pos = get_l1_position(dttime, units=self.pos.h['Units'], refframe=self.pos.h['ReferenceFrame'])
        r_ratio = L1Pos['r']/self.pos['r']
        
        #if 'density' in self.vars:
         #   self['density'] = self['density'] * (r_ratio)**(-2)
        
        #if 'btot' in self.vars:
        #    self['btot'] = self['btot'] * (r_ratio)**(-1.55) #original: 1.49
        #    #btot_upper = self['btot'] * (r_ratio)**(-1.)
        #    btot_lower = self['btot'] * (r_ratio)**(-2.)
        #    self['btot_err'] = np.abs(self['btot'] - btot_lower)
        print('scaling factor= ', -1.66)
        shift_vars_r = ['br', 'bx'] # radial component
        shift_vars = [v for v in shift_vars_r if v in self.vars]      # behave according to 1/r
        for var in shift_vars:
            self[var] = np.mean([self[var], self[var] * (r_ratio)**(-1.66)], axis=0) #original: 1.94
            self[var+'_err'] = np.std([self[var], self[var] * (r_ratio)**(-1.66)], axis=0)
        
        shift_vars_t = ['bt', 'by'] # tangential component
        shift_vars = [v for v in shift_vars_t if v in self.vars]      # behave according to 1/r
        for var in shift_vars:
            self[var] = np.mean([self[var], self[var] * (r_ratio)**(-1.66)], axis=0) #original: 1.26
            self[var+'_err'] = np.std([self[var], self[var] * (r_ratio)**(-1.66)], axis=0)

        shift_vars_n = ['bn', 'bz'] # normal component
        shift_vars = [v for v in shift_vars_n if v in self.vars]      # behave according to 1/r
        for var in shift_vars:
            self[var] = np.mean([self[var], self[var] * (r_ratio)**(-1.66)], axis=0) #original: 1.34
            self[var+'_err'] = np.std([self[var], self[var] * (r_ratio)**(-1.66)], axis=0)
            
        if 'btot' in self.vars:
            self['btot'] = np.sqrt(self['bx']**2 + self['by']**2 + self['bz']**2)
            #btot_lower = self['btot'] * (r_ratio)**(-2.)
            self['btot_err'] = np.sqrt((self['bx']/np.sqrt(self['bx']**2 + self['by']**2 + self['bz']**2)*self['bx_err'])**2+(self['by']/np.sqrt(self['bx']**2 + self['by']**2 + self['bz']**2)*self['by_err'])**2+(self['bz']/np.sqrt(self['bx']**2 + self['by']**2 + self['bz']**2)*self['bz_err'])**2)


        logger.info("shift_wind_to_L1: Scaled B and density values to L1 distance")

        return self
    
    
    
    def make_dst_prediction(self, dst1, dst2, dst3, method='temerin_li', t_correction=False, minute_res=True):
        """Makes prediction with data in array.

        Parameters
        ==========
        method : str
            Options = ['burton', 'obrien', 'temerin_li', 'temerin_li_2006']
        t_correction : bool
            For TL-2006 method only. Add a time-dependent linear correction to
            Dst values (required for anything beyond 2002).

        Returns
        =======
        dstData : new SatData obj
            New object containing predicted Dst data.
        """

        if method.lower() == 'temerin_li':
            if 'speedx' in self.vars:
                vx = self['speedx']
            else:
                vx = self['speed']
            logger.info("Calculating Dst for {} using Temerin-Li model 2002 version (updated parameters)".format(self.source))
            dst_pred = calc_dst_temerin_li(self['time'], self['btot'], self['bx'], self['by'], self['bz'], self['speed'], vx, self['density'], version='2002n')
        elif method.lower() == 'temerin_li_2002':
            if 'speedx' in self.vars:
                vx = self['speedx']
            else:
                vx = self['speed']
            logger.info("Calculating Dst for {} using Temerin-Li model 2002 version".format(self.source))
            dst_pred = calc_dst_temerin_li(self['time'], self['btot'], self['bx'], self['by'], self['bz'], self['speed'], vx, self['density'], version='2002')
        elif method.lower() == 'temerin_li_2006':
            if 'speedx' in self.vars:
                vx = self['speedx']
            else:
                vx = self['speed']
            logger.info("Calculating Dst for {} using Temerin-Li model 2006 version".format(self.source))
            dst_pred = calc_dst_temerin_li(self['time'], self['btot'], self['bx'], self['by'], self['bz'], self['speed'], vx, self['density'], dst1, dst2, dst3,
                                           version='2006', linear_t_correction=t_correction, minute_res=minute_res)
        elif method.lower() == 'obrien':
            logger.info("Calculating Dst for {} using OBrien model".format(self.source))
            dst_pred = calc_dst_obrien(self['time'], self['bz'], self['speed'], self['density'])
        elif method.lower() == 'burton':
            logger.info("Calculating Dst for {} using Burton model".format(self.source))
            dst_pred = calc_dst_burton(self['time'], self['bz'], self['speed'], self['density'])

        dstData = SatData({'time': copy.deepcopy(self['time']), 'dst': dst_pred})
        dstData.h['DataSource'] = "Dst prediction from {} data using {} method".format(self.source, method)
        dstData.h['SamplingRate'] = 1./24.

        return dstData


# =======================================================================================
# ----------------------------- III. CLASS PositionData() -------------------------------
# =======================================================================================
    

class PositionData():
    """Data object containing satellite position data.

    Init Parameters
    ===============
    --> PositionData(input_dict, source=None, header=None)
    posdata : list(x,y,z) or list(r,lon,lat)
        Dict containing the input data in the form of key: data (in array or list)
        Example: {'time': timearray, 'bx': bxarray}. The available keys for data input
        can be accessed in SatData.default_keys.
    header : dict(headerkey: value)
        Dict containing metadata on the data array provided. Useful data headers are
        provided in SatData.empty_header but this can be expanded as needed.
    source : str
        Provide quick-access name of satellite/data type for source.

    Attributes
    ==========
    .positions : np.ndarray
        Array containing position information. Best accessed using SatData[key].
    .h : dict
        Dict of metadata as defined by input header.

    Methods
    =======
    ...

    Examples
    ========
    """

    empty_header = {'ReferenceFrame': '',
                    'CoordinateSystem': '',
                    'Units': '',
                    'Object': '',
                    }
    
    # -----------------------------------------------------------------------------------
    # Internal methods for class PositionData()
    # -----------------------------------------------------------------------------------

    def __init__(self, posdata, postype, header=None):
        """Create new instance of class."""

        if not postype.lower() in ['xyz', 'rlonlat']:
            raise Exception("PositionData __init__: postype must be either 'xyz' or 'rlonlat'!")
        self.positions = np.asarray(posdata)
        if header == None:               # Inititalise empty header
            self.h = copy.deepcopy(PositionData.empty_header)
        else:
            self.h = header
        self.h['CoordinateSystem'] = postype.lower()
        self.coors = ['x','y','z'] if postype == 'xyz' else ['r','lon','lat']


    def __getitem__(self, var):
        if isinstance(var, str):
            if var in self.coors:
                return self.positions[self.coors.index(var)]
            else:
                raise Exception("PositionData object does not contain data under the key '{}'!".format(var))
        return self.positions[:,var]


    def __setitem__(self, var, value):
        if isinstance(var, str):
            if var in self.coors:
                self.positions[self.coors.index(var)] = value
            else:
                raise Exception("PositionData object does not contain the key '{}'!".format(var))
        else:
            raise ValueError("Cannot interpret {} as index for __setitem__!".format(var))


    def __len__(self):
        return len(self.positions[0])


    def __str__(self):
        return self.positions.__str__()

    
    def interp_to_time(self, t_orig, t_new):
        """Linearly interpolates over nans in array.

        Parameters
        ==========
        t_orig : np.ndarray
            Array containing original timesteps.
        t_new : np.ndarray
            Array containing new timesteps.
        """

        na = []
        for k in self.coors:
            na.append(np.interp(t_new, t_orig, self[k]))

        # Create new data opject:
        newData = PositionData(na, copy.deepcopy(self.h['CoordinateSystem']), header=copy.deepcopy(self.h))

        return newData

 
# =======================================================================================
# ----------------------------- III. COORDINATE CONVERSIONS------------------------------
# =======================================================================================

    
def convert_RTN_to_GSE_sta_l1(sc_in):

    sc=copy.deepcopy(sc_in)

    print('conversion RTN to GSE') 

    heeq_bx=np.zeros(len(sc))
    heeq_by=np.zeros(len(sc))
    heeq_bz=np.zeros(len(sc))

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))


    ########## first RTN to HEEQ 

    #go through all data points
    for i in np.arange(0,len(sc)):
        #print(sc.time[i])
        #HEEQ vectors
        X_heeq=[1,0,0]
        Y_heeq=[0,1,0]
        Z_heeq=[0,0,1]

        #normalized X RTN vector
        Xrtn=[sc.x[i],sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i],sc.y[i],sc.z[i]])
        #solar rotation axis at 0, 0, 1 in HEEQ
        Yrtn=np.cross(Z_heeq,Xrtn)/np.linalg.norm(np.cross(Z_heeq,Xrtn))
        Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))

        #project into new system
        heeq_bx[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),X_heeq)
        heeq_by[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Y_heeq)
        heeq_bz[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Z_heeq)

        #then HEEQ to GSE
        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5))  #use modified julian date  

        #then lambda_sun
        T00=(mjd[i]-51544.5)/36525.0
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   
        LAMBDA=280.460+36000.772*T00+0.04107*UT
        M=357.528+35999.050*T00+0.04107*UT

        #lt2 is lambdasun in Hapgood, equation 5, here in rad
        lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180

        #note that some of these equations are repeated later for the GSE to GSM conversion
        S1=np.matrix([[np.cos(lt2+np.pi), np.sin(lt2+np.pi),  0], [-np.sin(lt2+np.pi) , np.cos(lt2+np.pi) , 0], [0,  0,  1]])
        #create S2 matrix with angles with reversed sign for transformation HEEQ to HAE
        omega_node=(73.6667+0.013958*((mjd[i]+3242)/365.25))*np.pi/180 #in rad
        S2_omega=np.matrix([[np.cos(-omega_node), np.sin(-omega_node),  0], [-np.sin(-omega_node) , np.cos(-omega_node) , 0], [0,  0,  1]])
        inclination_ecl=7.25*np.pi/180
        S2_incl=np.matrix([[1,0,0],[0,np.cos(-inclination_ecl), np.sin(-inclination_ecl)], [0, -np.sin(-inclination_ecl), np.cos(-inclination_ecl)]])
        #calculate theta
        theta_node=np.arctan(np.cos(inclination_ecl)*np.tan(lt2-omega_node)) 

        #quadrant of theta must be opposite lt2 - omega_node Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360   
        lambda_omega_deg=np.mod(lt2-omega_node,2*np.pi)*180/np.pi
        x = np.cos(np.deg2rad(lambda_omega_deg))
        y = np.sin(np.deg2rad(lambda_omega_deg))

        #get theta_node in deg
        theta_node_deg=theta_node*180/np.pi
        x_theta = np.cos(theta_node)
        y_theta = np.sin(theta_node)
        #if in same quadrant, then theta_node = theta_node +pi   
        #if abs(lambda_omega_deg-theta_node_deg) < 180: theta_node=theta_node+np.pi ------> diese Zeile mit den if-Schleifen drunter ersetzen.
        if (x>=0 and y>=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node - np.pi
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node - np.pi/2
            elif (x_theta<=0 and y_theta>=0): theta_node = np.pi+(theta_node-np.pi/2)

        elif (x<=0 and y<=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node + np.pi
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node + np.pi/2
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node-np.pi/2

        elif (x>=0 and y<=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node + np.pi/2
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node - np.pi/2
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node + np.pi
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node

        elif (x<0 and y>0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node - np.pi/2
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node + np.pi/2
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node -np.pi          


        S2_theta=np.matrix([[np.cos(-theta_node), np.sin(-theta_node),  0], [-np.sin(-theta_node) , np.cos(-theta_node) , 0], [0,  0,  1]])

        #make S2 matrix
        S2=np.dot(np.dot(S2_omega,S2_incl),S2_theta)
        #this is the matrix S2^-1 x S1
        HEEQ_to_HEE_matrix=np.dot(S1, S2)
        #convert HEEQ components to HEE
        HEEQ=np.matrix([[heeq_bx[i]],[heeq_by[i]],[heeq_bz[i]]]) 
        HEE=np.dot(HEEQ_to_HEE_matrix,HEEQ)
        #change of sign HEE X / Y to GSE is needed
        sc.bx[i]=-HEE[0]
        sc.by[i]=-HEE[1]
        sc.bz[i]=HEE[2]

    print('conversion RTN to GSE done') 

    return sc

def convert_GSE_to_GSM_new(sc_in):

    sc=copy.deepcopy(sc_in)

    print('conversion GSE to GSM')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))

    for i in np.arange(0,len(sc)):
        #get all dates right
        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5))  #use modified julian date 

        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours    

        #define position of geomagnetic pole in GEO coordinates
        pgeo=78.8+4.283*((mjd[i]-46066)/365.25)*0.01 #in degrees
        lgeo=289.1-1.413*((mjd[i]-46066)/365.25)*0.01 #in degrees

        #GEO vector
        Qg=[np.cos(pgeo*np.pi/180)*np.cos(lgeo*np.pi/180), np.cos(pgeo*np.pi/180)*np.sin(lgeo*np.pi/180), np.sin(pgeo*np.pi/180)]

        #CREATE T1, T00, UT is known from above
        zeta=(100.461+36000.770*T00+15.04107*UT)*np.pi/180
        T1=np.matrix([[np.cos(zeta), np.sin(zeta),  0], [-np.sin(zeta) , np.cos(zeta) , 0], [0,  0,  1]]) #angle for transpose

        #lambda_sun in Hapgood, equation 5, here in degrees
        LAMBDA=280.460+36000.772*T00+0.04107*UT
        M=357.528+35999.050*T00+0.04107*UT
        lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180

        #CREATE T2, LAMBDA, M, lt2 known from above
        ##################### lamdbda und Z
        t2z=np.matrix([[np.cos(lt2), np.sin(lt2),  0], [-np.sin(lt2) , np.cos(lt2) , 0], [0,  0,  1]])
        et2=(23.439-0.013*T00)*np.pi/180

        ###################### epsilon und x
        t2x=np.matrix([[1,0,0],[0,np.cos(et2), np.sin(et2)], [0, -np.sin(et2), np.cos(et2)]])
        T2=np.dot(t2z,t2x)  #equation 4 in Hapgood 1992
        #matrix multiplications   
        T2T1t=np.dot(T2,np.matrix.transpose(T1))
        ################
        Qe=np.dot(T2T1t,Qg) #Q=T2*T1^-1*Qq
        psigsm=np.arctan(Qe.item(1)/Qe.item(2)) #arctan(ye/ze) in between -pi/2 to +pi/2

        T3=np.matrix([[1,0,0],[0,np.cos(-psigsm), np.sin(-psigsm)], [0, -np.sin(-psigsm), np.cos(-psigsm)]])
        b_gse=np.matrix([[sc.bx[i]],[sc.by[i]],[sc.bz[i]]]) 
        b_gsm=np.dot(T3,b_gse)   #equation 6 in Hapgood
        sc.bx[i]=b_gsm[0]
        sc.by[i]=b_gsm[1]
        sc.bz[i]=b_gsm[2]

    print('conversion GSE to GSM done')   

    return sc


    
# -----------------------------------------------------------------------------------
# Other functions
# -----------------------------------------------------------------------------------    
    
    
def expand_icme(df_timeshifted, l1, arr_time, t_le, t_te, power=0.8):

    df_s = df_timeshifted.copy(deep=True)

    mo_mask = (df_timeshifted['time'] >= date2num(t_le)) & (df_timeshifted['time'] <= date2num(t_te))
    prior_mask = (df_timeshifted['time'] < date2num(t_le))
    post_mask = (df_timeshifted['time'] > date2num(t_te))

    D1 = (t_te - t_le).total_seconds()
    FR = df_timeshifted[mo_mask]
    r1 = FR['r'].mean()

    #idx = df_s.set_index('time').index.get_loc(t_le, method='nearest')
    #ts_le = df_s['time_shifted'].iloc[idx]

    #idx2 = def_ref_sc.set_index('time').index.get_loc(ts_le, method='nearest')
    dct=(arr_time.replace(tzinfo=None))-l1.time
    idx2 = np.argmin(np.abs(dct))
    r2 = l1.r[idx2]
    D2 = D1 * (r2/r1)**power

    expansion_delta = np.linspace(0, len(df_timeshifted[mo_mask])-1, len(df_timeshifted[mo_mask]))*60*(D2/D1)

    df_timeshifted['expansion_delta'] = np.nan
    df_mo = df_timeshifted[mo_mask].assign(expansion_delta=expansion_delta)
    df_prior = df_timeshifted[prior_mask].assign(expansion_delta=expansion_delta.min())
    df_post = df_timeshifted[post_mask].assign(expansion_delta=expansion_delta.max())

    stitched_df = pd.concat([df_prior, df_mo, df_post])

    t = []
    for i in range(len(stitched_df)):
        new_t = num2date(stitched_df['time_shifted'].iloc[i]) + timedelta(seconds=stitched_df['expansion_delta'].iloc[i])
        t = np.append(t, new_t)
    stitched_df['time_shifted_exp'] = t

    return stitched_df