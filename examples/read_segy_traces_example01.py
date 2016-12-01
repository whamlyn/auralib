"""
Example script to illustrate how auralib may be used to read trace data from
SEG-Y format data files.

Written by: Wes Hamlyn
Created:    1-Dec-2016
"""

import auralib as aura

infile = r'C:\Users\whamlyn\Dropbox\data\QSI\data\Project_data\cdps_line2.sgy'

def_bhead = {'samp_rate':{'bpos':17, 'fmt':'h', 'nbyte':2},
             'num_samp':{'bpos':21, 'fmt':'h', 'nbyte':2},
             'samp_fmt':{'bpos':25, 'fmt':'h', 'nbyte':2}
             }

def_thead = {'il':{'bpos':189,  'fmt':'l', 'nbyte':4},
             'xl':{'bpos':193, 'fmt':'l', 'nbyte':4},
             'cmpx':{'bpos':181, 'fmt':'l', 'nbyte':4},
             'cmpy':{'bpos':185, 'fmt':'l', 'nbyte':4},
             'offset':{'bpos':37, 'fmt':'l', 'nbyte':4}
             }
			 
segybuf = aura.segy.Segy(infile, def_bhead, def_thead)

print(segybuf.ebcdic)