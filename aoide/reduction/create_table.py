#!/usr/bin/env python
import os
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import sys
import numpy
import optparse


parser = optparse.OptionParser(usage="%prog -d dir -t table")
parser.add_option("-d", "--dir", dest="dir_in", type="string", help="directory with all the raw files")
parser.add_option("-t", "--table", dest="table_out", type="string", help="file nname of the output fits file")

(options, args) = parser.parse_args()
dir_in = options.dir_in
table_out = options.table_out


files = os.listdir(dir_in)
MUSE_files = []
MUSE_name = []
MUSE_type = []
OB_ID = []
EXPTIME = [] 
NEXP =[]
EXPNO=[]
time =[]
mode=[]
mjd=[]

for f in files:
 if 'MUSE' in f and '.fits.fz'  in f:
     header = pyfits.getheader("%s%s"%(dir_in,f),0)
     MUSE_files.append("%s/%s"%(dir_in,f))
     MUSE_name.append(header['OBJECT'])
     MUSE_type.append(header['ESO DPR TYPE'])
     NEXP.append(int(header['ESO TPL NEXP']))
     EXPNO.append(int(header['ESO TPL EXPNO']))
     time.append(header['ESO TPL START'])
     OB_ID.append(int(header['ESO OBS ID']))
     mode.append(header['ESO INS MODE'])
     mjd.append(float(header['MJD-OBS']))
     EXPTIME.append(float(header['EXPTIME']))

MUSE_files = numpy.array(MUSE_files)
MUSE_name = numpy.array(MUSE_name)
MUSE_type = numpy.array(MUSE_type)
NEXP = numpy.array(NEXP)
EXPNO = numpy.array(EXPNO)
time = numpy.array(time)
OB_ID = numpy.array(OB_ID)
mode = numpy.array(mode)
mjd = numpy.array(mjd)
EXPTIME = numpy.array(EXPTIME)
idxsort = numpy.argsort(MUSE_files)

columns=[]
columns.append(pyfits.Column(name='file',format='100A',array=MUSE_files[idxsort]))
columns.append(pyfits.Column(name='name',format='20A',array=MUSE_name[idxsort]))
columns.append(pyfits.Column(name='type',format='20A',array=MUSE_type[idxsort]))
columns.append(pyfits.Column(name='time',format='30A',array=time[idxsort]))
columns.append(pyfits.Column(name='mjd',format='E',array=mjd[idxsort]))
columns.append(pyfits.Column(name='EXPTIME',format='E',array=EXPTIME[idxsort]))
columns.append(pyfits.Column(name='mode',format='10A',array=mode[idxsort]))
columns.append(pyfits.Column(name='OB_ID',format='K',array=OB_ID[idxsort]))
columns.append(pyfits.Column(name='NEXP',format='J',array=NEXP[idxsort]))
columns.append(pyfits.Column(name='EXPNO',format='J',array=EXPNO[idxsort]))

table = pyfits.new_table(columns)
table.writeto(table_out,clobber=True)




     
