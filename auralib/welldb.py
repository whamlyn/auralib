# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 14:46:29 2020

@author: wesha
"""

#%%

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import os
import shutil


def create_well(welldb, wellid, overwrite=False):
    """
    Function to create a new well
    """
    curwellpath = os.path.join(welldb, wellid)

    if os.path.exists(curwellpath):
        if overwrite == False:
            print('Path already exists:\n%s' % curwellpath)
            ans = input('Overwrite? (Y/[N]): ')

            if (len(ans)==0) | (ans.lower()=='n'):
                print('Skipping well %s' % curwellpath)
            else:
                overwrite = True

        elif overwrite==True:
            shutil.rmtree(curwellpath)

    os.mkdir(curwellpath)
    os.mkdir(os.path.join(curwellpath, 'logs'))
    os.mkdir(os.path.join(curwellpath, 'dev'))
    os.mkdir(os.path.join(curwellpath, 'td'))
    os.mkdir(os.path.join(curwellpath, 'markers'))
    os.mkdir(os.path.join(curwellpath, 'synth'))
    os.mkdir(os.path.join(curwellpath, 'gathers'))

    well_header = {'WELLID': wellid, 'TOPX': 0.0, 'TOPY': 0.0,
                   'KB': 0.0, 'GL': 0.0}
    with open(os.path.join(curwellpath, 'well_header.csv'), 'w') as fd:
        for key in well_header:
            line = '%s,%s\n' % (key, well_header[key])
            fd.write(line)
    
    # initialize markers file
    with open(os.path.join(curwellpath, 'markers', 'markers.csv'), 'w') as fd:
        fd.write('WELLID,%s\n' % wellid)
        fd.write('#NAME,MD\n')



def create_las_logs(welldb, lasfile, setname='LAS', ztype='md', wellidfield='WELL',
                    create_missing_well=True):
    """
    Function to load logs from LAS to a database well
    """

    # create LASReader object
    buf = aura.las.LASReader(lasfile)

    # get some important constants
    wellid = buf.well_info[wellidfield]['data']
    strt = buf.well_info['STRT']['data']
    stop = buf.well_info['STOP']['data']
    step = buf.well_info['STEP']['data']
    zunit = buf.well_info['STRT']['unit'].lower()
    ztype = ztype.lower()

    # If the well doesn't exist in the database, create it now
    if os.path.exists(os.path.join(welldb, wellid)) == False:
        if create_missing_well == True:
            print('Well %s does not exist in DB. Creating now...\n' % wellid)
            create_well(welldb, wellid)
        else:
            print('Well %s does not exist in DB. Skipping...\n' % wellid)

    # If the log SET directory doesn't exist, create it now too
    logsdir = os.path.join(welldb, wellid, 'logs')
    setdir = os.path.join(logsdir, setname)
    if os.path.exists(setdir) == False:
        os.mkdir(setdir)

    # Write logs containg in LAS file to the log SET directory
    for logname in buf.curves.keys():
        logunit = buf.curve_info[logname]['unit']
        logdesc = buf.curve_info[logname]['desc']
        logfile = os.path.join(setdir, logname+'.csv')
        with open(logfile, 'w') as fd:
            fd.write('#WELL,%s\n' % wellid)
            fd.write('#ZTYPE,%s\n' % ztype)
            fd.write('#ZUNIT,%s\n' % zunit)
            fd.write('#STRT,%f\n' % strt)
            fd.write('#STOP,%f\n' % stop)
            fd.write('#STEP,%f\n' % step)
            fd.write('#UNIT,%s\n' % logunit)
            fd.write('#DESC,%s\n' % logdesc)
            fd.write('~DATA\n')
            zlog = np.arange(strt, stop+step, step)
            for z, data in zip(zlog, buf.curves[logname]):
                fd.write('%f,%f\n' % (z, data))


def add_top(welldb, wellid, z, name):
    """
    Add well markers for well.
    """
    
    well_exists = os.path.exists(os.path.join(welldb, wellid))
    
    if well_exists:
        
        tops_file = os.path.join(welldb, wellid, 'markers', 'markers.csv')
    
        with open(tops_file, 'a') as fd:
            
            fd.seek(0, 2)
            
            if type(name) == str:
                line = '%s,%f\n' % (name, z)
                fd.write(line)
        
            else:
                for i in range(len(z)):
                    line = '%s,%f\n' % (name[i], z[i])
                    fd.write(line)
    else:
        print('Well %s does not exist in database, skipping...\n' % wellid)
    
        
def update_well_header(welldb, wellid, hdr_name, hdr_value):
    """
    Update and/or add new well header fields
    """
    
    curwellpath = os.path.join(welldb, wellid)
    well_exists = os.path.exists(curwellpath)
    
    if well_exists:
        with open(os.path.join(curwellpath, 'well_header.csv'), 'a') as fd:
            buf = fd.readlines()
        
        header = {}
        for line in buf:
            line = line.split()
            header[line[0]] = line[1]
        
        for i in range(len(hdr_name)):
            if type(hdr_value[i]) == str:
                header[hdr_name[i]] = hdr_value[i]
            else:
                header[hdr_name[i]] = '%f' % hdr_value[i]
            
        with open(os.path.join(curwellpath, 'well_header.csv'), 'w') as fd:
            for key in header.keys():
                line = '%s,%s\n' % (key, header[key])
                fd.write(line)
    
    
    
    