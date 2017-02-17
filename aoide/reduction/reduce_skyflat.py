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
parser.add_option("-f", "--flag", dest="flag", type="string", help="Either ALL or the date of the template to be reduced")
parser.add_option("-c", "--cores", dest="cores", type="int", help="number of cores to be used")


(options, args) = parser.parse_args()
raw_table=options.raw_table
dir_out=options.dir_out
dir_temp=options.dir_temp
flag=options.flag
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

select_bias = (MUSE_type=='BIAS') &  (MUSE_name=='BIAS') & (EXPNO==1)
select_flat = (MUSE_type=='FLAT,LAMP') &  (MUSE_name=='FLAT,LAMP') & (EXPNO==1)
select_illum = (MUSE_type=='FLAT,LAMP,ILLUM') &  (MUSE_name=='FLAT,LAMP,ILLUM') & (EXPNO==1)
select_wave = (MUSE_type=='WAVE') &  (MUSE_name=='WAVE') & (EXPNO==1)
select_skyflat = (MUSE_type=='FLAT,SKY') &  (MUSE_name=='FLAT,SKY')

if flag=='ALL':
    times = numpy.unique(time[select_skyflat])
    for i in range(len(times)):
        select_TPL = (MUSE_type=='FLAT,SKY') &  (MUSE_name=='FLAT,SKY') & (time==times[i])
        mjd_sky = mjd[select_TPL][0]
        if NEXP[select_TPL][0]==numpy.sum(select_TPL):
            tags={'SKYFLAT':MUSE_files[select_TPL]}
            mjd_biases =  mjd[select_bias]
            time_biases = time[select_bias]
            idx_bias = numpy.argmin((mjd_biases-mjd_sky)**2)
            tags['MASTER_BIAS']=['%s/BIAS/MASTER_BIAS_%s.fits'%(dir_out,time_biases[idx_bias])]
            mjd_illum = mjd[select_illum]
            time_illum = time[select_illum]
            idx_illum = numpy.argmin((mjd_illum-mjd_sky)**2)
            tags['ILLUM']=[MUSE_files[select_illum][idx_illum]]
            mjd_wave =  mjd[select_wave]
            time_wave = time[select_wave]
            idx_wave = numpy.argmin((mjd_wave-mjd_sky)**2)
            tags['WAVECAL_TABLE']=['%s/WAVE/WAVECAL_TABLE_%s.fits'%(dir_out,time_wave[idx_wave])]
            mjd_flat =  mjd[select_flat]
            time_flat = time[select_flat]
            idx_flat = numpy.argmin((mjd_flat-mjd_sky)**2)
            tags['TRACE_TABLE']=['%s/TRACE/TRACE_TABLE_%s.fits'%(dir_out,time_flat[idx_flat])]
            tags['MASTER_FLAT']=['%s/FLAT/MASTER_FLAT_%s.fits'%(dir_out,time_flat[idx_flat])]
            
            muse.create_sof_twilight('%s/twilight.sof'%(dir_temp),tags,time[i].split('T')[0])
        
        cwd=os.getcwd()
        os.chdir(dir_temp)
        muse.compute_twilight('twilight.sof',cores)
        if muse.nocache:
            os.system('nocache cp TWILIGHT_CUBE.fits %s/SKYFLAT/TWILIGHT_CUBE_%s.fits'%(dir_out,times[i]))
        else:
            os.system('cp TWILIGHT_CUBE.fits %s/SKYFLAT/TWILIGHT_CUBE_%s.fits'%(dir_out,times[i]))
        os.system('rm DATACUBE_SKYFLAT*.fits')
        os.system('rm TWILIGHT_CUBE*.fits')
        os.chdir(cwd)
else:
    select_TPL = (MUSE_type=='FLAT,SKY') &  (MUSE_name=='FLAT,SKY') & (time==flag)
    mjd_sky = mjd[select_TPL][0]
#    if NEXP[select_TPL][0]==numpy.sum(select_TPL)-1:
    try:
        tags={'SKYFLAT':MUSE_files[select_TPL]}
        mjd_biases =  mjd[select_bias]
        time_biases = time[select_bias]
        idx_bias = numpy.argmin((mjd_biases-mjd_sky)**2)
        tags['MASTER_BIAS']=['%s/BIAS/MASTER_BIAS_%s.fits'%(dir_out,time_biases[idx_bias])]
        mjd_illum = mjd[select_illum]
        time_illum = time[select_illum]
        idx_illum = numpy.argmin((mjd_illum-mjd_sky)**2)
        tags['ILLUM']=[MUSE_files[select_illum][idx_illum]]
        mjd_wave =  mjd[select_wave]
        time_wave = time[select_wave]
        idx_wave = numpy.argmin((mjd_wave-mjd_sky)**2)
        tags['WAVECAL_TABLE']=['%s/WAVE/WAVECAL_TABLE_%s.fits'%(dir_out,time_wave[idx_wave])]
        mjd_flat =  mjd[select_flat]
        time_flat = time[select_flat]
        idx_flat = numpy.argmin((mjd_flat-mjd_sky)**2)
        tags['TRACE_TABLE']=['%s/TRACE/TRACE_TABLE_%s.fits'%(dir_out,time_flat[idx_flat])]
        tags['MASTER_FLAT']=['%s/FLAT/MASTER_FLAT_%s.fits'%(dir_out,time_flat[idx_flat])]
	print dir_temp
        muse.create_sof_twilight('%s/twilight.sof'%(dir_temp),tags,flag.split('T')[0])
    except IOError:
	pass     
    cwd=os.getcwd()
    os.chdir(dir_temp)
    muse.compute_twilight('twilight.sof',cores)
    if muse.nocache:
        os.system('nocache cp TWILIGHT_CUBE.fits %s/SKYFLAT/TWILIGHT_CUBE_%s.fits'%(dir_out,flag))
    else:
        os.system('cp TWILIGHT_CUBE.fits %s/SKYFLAT/TWILIGHT_CUBE_%s.fits'%(dir_out,flag))
    os.system('rm DATACUBE_SKYFLAT*.fits')
    os.system('rm TWILIGHT_CUBE*.fits')
    os.chdir(cwd)
     





     
