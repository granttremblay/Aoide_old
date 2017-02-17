#!/usr/bin/env python
import os
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import sys
import numpy
import muse
import optparse


parser = optparse.OptionParser(usage="%prog -o dir_out -t table -f flag")
parser.add_option("-o", "--dir_out", dest="dir_out", type="string", help="directory with all the computed Master calibration and reduced files")
parser.add_option("-t", "--table", dest="raw_table", type="string", help="fits file with all the raw data infos")
parser.add_option("-f", "--flag", dest="flag", type="string", help="Name of the object to be reduced.")


(options, args) = parser.parse_args()
raw_table=options.raw_table
dir_out=options.dir_out
name=options.flag


hdu = pyfits.open(raw_table)
tab = hdu[1].data

MUSE_files = tab.field('file')
MUSE_name = tab.field('name')
MUSE_type = tab.field('type')
NEXP = tab.field('NEXP')
EXPNO = tab.field('EXPNO')
time = tab.field('time')
mjd=tab.field('mjd')


select_obj = (MUSE_type=='OBJECT') &  (MUSE_name==name)


m=1
times = numpy.unique(time[select_obj])
for i in range(len(times)):
    select_sky = (MUSE_type=='SKY') &  (MUSE_name==name) & (time==times[i])
    select_obj_OB = (MUSE_type=='OBJECT') &  (MUSE_name==name) & (time==times[i])
    mjd_sky = mjd[select_sky]
    expno_obj = EXPNO[select_obj_OB]
    expno_sky = EXPNO[select_sky]
    SKY_in = ['%s/%s/SKY_SPECTRUM_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno_sky[s]) for s in range(len(expno_sky))]
    
    for o in range(len(expno_obj)):
        print expno_obj[o]
        hdr = pyfits.getheader(MUSE_files[select_obj_OB][o],0)
        muse.interpolate_sky(SKY_in,mjd_sky,'%s/%s/SKY_SPECTRUM_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno_obj[o]),mjd[select_obj_OB][o],header_out=hdr)
 
    

	

     





     
