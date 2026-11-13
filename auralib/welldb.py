# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:23:13 2021

@author: whamlyn1
"""

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import os
import datetime
import pandas as pd
from scipy.interpolate import interp1d, PchipInterpolator


class WellDB():
    """
    """
    
    def __init__(self, dbpath):
        """
        Constructor method 
        """
        
        self.dbpath = dbpath
    
    
    def create_new_welldb(self):
        """
        Create a new well database
        """
        
        try:
            os.makedirs(self.dbpath, exist_ok=False)
    
        except Exception as e:
            print(e)
            
        with open(os.path.join(self.dbpath, 'description.txt'), 'w') as fd:
            fd.write('Aura WellDB Version 1.0\n')
            fd.write('Created on: %s' % datetime.date.today())
    
    
    def create_new_well(self, wellid):
        """
        Initialize a new well in the current database
        """
        
        wellpath = os.path.join(self.dbpath, wellid)
        
        try:
            os.makedirs(wellpath)
        except Exception as e:
            print(e)
        
        try:
            for subdir in ['log', 'point', 'td', 'dev']:
                os.makedirs(os.path.join(wellpath, subdir))
            
            with open(os.path.join(wellpath, 'well_header.txt'), 'w') as fd:
                fd.write('WELLID:%s\n' % wellid)
                fd.write('WELLNAME:\n')
                fd.write('GL:\n')
                fd.write('KB:\n')
                fd.write('X_TOP:\n')
                fd.write('X_BOT:\n')
            
            with open(os.path.join(wellpath, 'markers.csv'), 'w') as fd:
                fd.write('WELLID,NAME,ZREF,ZTYPE\n')
                
        except Exception as e:
            print(e)
    
    
    def create_new_log_set(self, setid, wellid):
        """
        Add a new log set to a well
        """
        
        setpath = os.path.join(self.dbpath, wellid, 'log', setid)
        
        try:
            os.makedirs(setpath)
        except Exception as e:
            print(e)
    
    
    def add_log_to_set(self, log, setid, wellid, overwrite=False):
        """"
        Add a new log to a set
        """
        
        logpath = os.path.join(self.dbpath, wellid, 'log', setid, log.name)
        
        if (os.path.exists(logpath)==True) & (overwrite==False):
            print('Log already exists: %s' % logpath)

        else:
            #print('Writing log...')
            #print('%s' % logpath)
            
            try:
                with open(logpath, 'w') as fd:
                    fd.write('mnem:%s\n' % log.name)
                    fd.write('type:%s\n' % log.type)
                    fd.write('unit:%s\n' % log.units)
                    fd.write('ztype:%s\n' % log.ztype)
                    fd.write('start:%f\n' % log.start)
                    fd.write('stop:%f\n' % log.stop)
                    fd.write('step:%f\n' % log.step)
                    fd.write('nsamp:%i\n' % log.nsamp)
                    
                    for sample in log.data:
                        if np.isnan(sample):
                            fd.write('%s\n' % 'nan')
                        else:
                            fd.write('%f\n' % sample)
                        
            except Exception as e:
                print(e)
    
    
    def get_log_from_set(self, logid, setid, wellid, returnz=True):
        """
        Read well log from set
        """
        
        logpath = os.path.join(self.dbpath, wellid, 'log', setid, logid)
        
        log = aura.welldb.AuraLog()
        
        try:
            with open(logpath, 'r') as fd:
                buf = fd.readlines()
            
            log.name = buf[0].split(':')[1].strip()
            log.type = buf[1].split(':')[1].strip()
            log.units = buf[2].split(':')[1].strip()
            log.ztype = buf[3].split(':')[1].strip()
            log.start = float(buf[4].split(':')[1].strip())
            log.stop =  float(buf[5].split(':')[1].strip())
            log.step =  float(buf[6].split(':')[1].strip())
            log.nsamp = int(buf[7].split(':')[1].strip())
            
            log.data = np.array(buf[8:], dtype=float)
            nanidx = np.nonzero(log.data==-999.25)
            log.data[nanidx] = np.nan
            
            log.zref = np.arange(log.nsamp)*log.step + log.start
            
            return log
        
        except Exception as e:
            print(e)
            
    
    def add_marker_to_well(self, wellid, mkr_name, mkr_zref, mkr_ztype):
        """
        """
        
        mkrpath = os.path.join(self.dbpath, wellid, 'markers.csv')
        outstr = '%s,%s,%f,%s\n' % (wellid, mkr_name, mkr_zref, mkr_ztype)

        try:
            with open(mkrpath, 'r+') as fd:
                fd.seek(0, 2)
                fd.write(outstr)
        
        except Exception as e:
            print(e)
            
    
    def get_markers_for_well(self, wellid):
        """
        """
        
        mkrpath = os.path.join(self.dbpath, wellid, 'markers.csv')
        
        df = pd.read_csv(mkrpath)
        df = df.set_index('NAME')
        
        mkrs = df.to_dict(orient='index')
        
        return mkrs

    
    
    def get_well_list(self):
        """
        Returns a list of wellids in the well database
        """
        
        wells = []
        for item in os.listdir(self.dbpath):
            itempath = os.path.join(self.dbpath, item)
            if os.path.isdir(itempath):
                wells.append(item)
        
        return wells
    
    def get_log_names(self, wellid, setid=None):
        """
        """

        wellpath = os.path.join(self.dbpath, wellid)
        logpath = os.path.join(wellpath, 'log')
        
        log_sets = os.listdir(logpath)
        
        logs = {}
        for cur_set in log_sets:
            log_names = os.listdir(os.path.join(logpath, cur_set))
            logs[cur_set] = log_names
        
        return logs
        



class AuraLog():
    """
    Aura Well Log Class
    """
    
    def __init__(self, name=None, ltype=None, data=None, units=None, ztype='MD',
                 start=None, step=None):

        self.name = name
        self.type = ltype
        self.data = data
        self.units = units
        self.ztype = ztype
        self.start = start
        self.step = step

        try:
            self.nsamp = len(data)
            self.stop = self.nsamp*step + self.start
            self.strip_nans()
            self.zref = np.arange(self.nsamp)*self.step + self.start
            
        except:
            self.nsamp = 0
            self.stop = self.start
            self.zref = None
        
    
    def strip_nans(self):
        """
        Strip nans from start and end of log. Internal nans will not be removed
        """
        
        z = self.start + np.arange(self.nsamp)*self.step
        
        mask = np.isnan(self.data)
        
        idx1 = mask.argmin()
        idx2 = mask.size - mask[::-1].argmin()
        
        data = self.data[idx1:idx2]
        start = z[idx1]
        stop = z[idx2-1]
        nsamp = len(data)
        
        self.data = data
        self.start = start
        self.stop = stop
        self.nsamp = nsamp
        
    def write_log(self, outpath, overwrite=False):
        """
        """
        
        if (os.path.exists(outpath)==True) & (overwrite==False):
            print('Log already exists: %s' % outpath)

        else:
            
            try:
                with open(outpath, 'w') as fd:
                    fd.write('mnem:%s\n'  % self.name)
                    fd.write('type:%s\n'  % self.type)
                    fd.write('unit:%s\n'  % self.units)
                    fd.write('ztype:%s\n' % self.ztype)
                    fd.write('start:%f\n' % self.start)
                    fd.write('stop:%f\n'  % self.stop)
                    fd.write('step:%f\n'  % self.step)
                    fd.write('nsamp:%i\n' % self.nsamp)
                    
                    for sample in self.data:
                        if np.isnan(sample):
                            fd.write('%s\n' % 'nan')
                        else:
                            fd.write('%f\n' % sample)
                        
            except Exception as e:
                print(e)
    
    def read_log(self, inpath):
        """
        """
        
        try:
            with open(inpath, 'r') as fd:
                buf = fd.readlines()
            
            self.name = buf[0].split(':')[1].strip()
            self.type = buf[1].split(':')[1].strip()
            self.units = buf[2].split(':')[1].strip()
            self.ztype = buf[3].split(':')[1].strip()
            self.start = float(buf[4].split(':')[1].strip())
            self.stop =  float(buf[5].split(':')[1].strip())
            self.step =  float(buf[6].split(':')[1].strip())
            self.nsamp = int(buf[7].split(':')[1].strip())
            
            self.data = np.array(buf[8:], dtype=float)
            nanidx = np.nonzero(self.data==-999.25)
            self.data[nanidx] = np.nan
            
        except Exception as e:
            print(e)
    
    
    def resample_log(self, new_zref):
        """
        Resample well log data from the original zref to a new zref.
        """
        
        start = new_zref[0]
        stop = new_zref[-1]
        nsamp = len(new_zref)
        step = (stop - start)/(nsamp - 1)
        
        f = interp1d(self.zref, self.data, kind='linear', bounds_error=False,
                     fill_value=(self.data[0], self.data[-1]))
        
        #f = PchipInterpolator(self.zref, self.data, extrapolate=False)
        
        data = f(new_zref)
        
        log = AuraLog(name=self.name, ltype=self.type, data=data, 
                      units=self.units, ztype=self.ztype, 
                      start=start, step=step)
        
        return log
    

def read_petrel_tops(infile):
    """
    Load tops from Petrel
    """
    
    with open(infile, 'r') as fd:
        buf = fd.readlines()
    
    fields = []
    tops = []
    markers = {}
    
    count = -1
    in_header = False
    in_data = False
    for line in buf:
        
        count += 1
        
        line = line.strip()
        
        if line[0] == '#':
            print('Comment line. Skipping...')
            continue
        
        if line == 'BEGIN HEADER':
            print('Reading Header Fields')
            in_header = True
            continue
    
        if line == 'END HEADER':
            print('Finished reading header')
            in_header = False
            in_data = True
            continue
            
        if in_header == True:
            fields.append(line)
            
        if in_data == True:
            tmp = line.split('"')
            tmp.remove(' ')
            tmp.remove(' ')
            tmp0 = tmp[0].strip().split()
            tmp4 = tmp[4].strip().split()
            data = []
            data.extend(tmp0)
            data.append(tmp[1])
            data.append(tmp[2])
            data.append(tmp[3])
            data.extend(tmp4)
            
            tops.append(data)
    
    all_wells = []
    for each in tops:
        all_wells.append(each[fields.index('Well')])
        
    for well in list(set(all_wells)):
            markers[well] = {}
            
    for i in range(len(tops)):
        curtop = tops[i]
        WELL = curtop[fields.index('Well')]
        MD = float(curtop[fields.index('MD')])
        TVDss = float(curtop[fields.index('Z')])*-1
        NAME = curtop[fields.index('Surface')]
        
        markers[WELL][NAME] = {'MD': MD, 'TVDss': TVDss}
    
    return markers