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


parser = optparse.OptionParser(usage="%prog -d dir_temp -o dir_out -t table -f flag -c cores -n nifu")
parser.add_option("-d", "--dir_temp", dest="dir_temp", type="string", help="directory for temporary execution")
parser.add_option("-o", "--dir_out", dest="dir_out", type="string", help="directory with all the computed Master calibration and reduced files")
parser.add_option("-t", "--table", dest="raw_table", type="string", help="fits file with all the raw data infos")
parser.add_option("-f", "--flag", dest="flag", type="string", help="Name of the object to be reduced.")
parser.add_option("-c", "--cores", dest="cores", type="int", help="number of cores to be used")


(options, args) = parser.parse_args()
raw_table=options.raw_table
dir_out=options.dir_out
dir_temp=options.dir_temp
name=options.flag
cores=options.cores


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
select_sky = (MUSE_type=='SKY') &  (MUSE_name==name)

os.system('mkdir %s/%s'%(dir_out,name.replace(' ','')))

FOVs=[]
PIXTABS=[]
for i in range(len(time[select_obj])):
   FOVs.append('%s/%s/IMAGE_FOV_%s_%02d.fits'%(dir_out,name.replace(' ',''),time[select_obj][i],EXPNO[select_obj][i]))
   PIXTABS.append('%s/%s/PIXTABLE_REDUCED_%s_%02d.fits'%(dir_out,name.replace(' ',''),time[select_obj][i],EXPNO[select_obj][i]))
    

cwd=os.getcwd()
os.chdir(dir_temp)   
tags={'IMAGE_FOV':FOVs}
muse.create_sof_align('align.sof',tags)
muse.compute_align('align.sof')
tags={'PIXTABLE_REDUCED':PIXTABS}
muse.create_sof_cube('cube.sof',tags)
muse.compute_cube('cube.sof',cores)

if muse.nocache:
    os.system('nocache cp DATACUBE_FINAL.fits %s/%s/'%(dir_out,name.replace(' ','')))
    os.system('nocache cp IMAGE_FOV_00*.fits %s/%s/'%(dir_out,name.replace(' ',''))) 
else:
    os.system('cp DATACUBE_FINAL.fits %s/%s/'%(dir_out,name.replace(' ','')))
    os.system('cp IMAGE_FOV_00*.fits %s/%s/'%(dir_out,name.replace(' ',''))) 
os.system('rm IMAGE_FOV_*.fits')
os.system('rm DATACUBE_FINAL.fits')

PIXTABS=[]
for i in range(len(time[select_sky])):
   PIXTABS.append('%s/%s/PIXTABLE_REDUCED_%s_%02d.fits'%(dir_out,name.replace(' ',''),time[select_sky][i],EXPNO[select_sky][i]))

tags={'PIXTABLE_REDUCED':PIXTABS}
muse.create_sof_cube('cube_sky.sof',tags,offset=False)
muse.compute_cube('cube_sky.sof',cores,filter_FOV='white')

if muse.nocache:
    os.system('nocache cp DATACUBE_FINAL.fits %s/%s/DATACUBE_SKY.fits'%(dir_out,name.replace(' ','')))
    os.system('nocache cp IMAGE_FOV_0001.fits %s/%s/IMAGE_FOV_SKY.fits'%(dir_out,name.replace(' ',''))) 
else:
     os.system('cp DATACUBE_FINAL.fits %s/%s/DATACUBE_SKY.fits'%(dir_out,name.replace(' ','')))
     os.system('cp IMAGE_FOV_0001.fits %s/%s/IMAGE_FOV_SKY.fits'%(dir_out,name.replace(' ',''))) 
os.system('rm IMAGE_FOV_*.fits')
os.system('rm DATACUBE_FINAL.fits')    





     
