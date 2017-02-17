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
parser.add_option("-s", "--save", dest="save", type="string", help="save option for muse_scipost", default="cube,individual")
parser.add_option("-f", "--flag", dest="flag", type="string", help="Name of the object to be reduced.")
parser.add_option("-c", "--cores", dest="cores", type="int", help="number of cores to be used")
parser.add_option("-n", "--nifu", dest="nifu", type="int", help="mode of reduction per ifu")


(options, args) = parser.parse_args()
raw_table=options.raw_table
dir_out=options.dir_out
dir_temp=options.dir_temp
name=options.flag
cores=options.cores
nifu=options.nifu
save=options.save


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
select_skyflat = (MUSE_type=='FLAT,SKY') &  (MUSE_name=='FLAT,SKY') & (EXPNO==1)
select_obj = (MUSE_type=='OBJECT') &  (MUSE_name==name)
select_std = (MUSE_type=='STD') & (MUSE_name=='STD')

os.system('mkdir %s/%s'%(dir_out,name.replace(' ','')))

m=1
times = numpy.unique(time[select_obj])
for i in range(len(times)):
    select_TPL = (MUSE_type=='OBJECT') &  (MUSE_name==name) & (time==times[i])
    expno = EXPNO[select_TPL]
    mjd_obj = mjd[select_TPL][0]
    tags={'OBJECT':MUSE_files[select_TPL]}
    mjd_biases =  mjd[select_bias]
    time_biases = time[select_bias]
    idx_bias = numpy.argmin((mjd_biases-mjd_obj)**2)
    tags['MASTER_BIAS']=['%s/BIAS/MASTER_BIAS_%s.fits'%(dir_out,time_biases[idx_bias])]
    mjd_illum = mjd[select_illum]
    time_illum = time[select_illum]
    idx_illum = numpy.argmin((mjd_illum-mjd_obj)**2)
    tags['ILLUM']=[MUSE_files[select_illum][idx_illum]]
    mjd_wave =  mjd[select_wave]
    time_wave = time[select_wave]
    idx_wave = numpy.argmin((mjd_wave-mjd_obj)**2)
    tags['WAVECAL_TABLE']=['%s/WAVE/WAVECAL_TABLE_%s.fits'%(dir_out,time_wave[idx_wave])]
    mjd_flat =  mjd[select_flat]
    time_flat = time[select_flat]
    idx_flat = numpy.argmin((mjd_flat-mjd_obj)**2)
    tags['TRACE_TABLE']=['%s/TRACE/TRACE_TABLE_%s.fits'%(dir_out,time_flat[idx_flat])]
    tags['MASTER_FLAT']=['%s/FLAT/MASTER_FLAT_%s.fits'%(dir_out,time_flat[idx_flat])]

    mjd_skyflat = mjd[select_skyflat]
    time_skyflat = time[select_skyflat]
    idx_skyflat = numpy.argmin((mjd_skyflat-mjd_obj)**2)
    tags['TWILIGHT_CUBE']=['%s/SKYFLAT/TWILIGHT_CUBE_%s.fits'%(dir_out,time_skyflat[idx_skyflat])]
          
    muse.create_sof_scibasic('%s/scibasic.sof'%(dir_temp),tags,geometry='',time=time[i].split('T')[0])

    mjd_std =  mjd[select_std]
    time_std = time[select_std]
    idx_std = numpy.argmin((mjd_std-mjd_obj)**2)
        
        
    cwd=os.getcwd()
    os.chdir(dir_temp)
    muse.compute_scibasic('scibasic.sof',cores,nifu,'TRUE')
    for j in range(len(MUSE_files[select_TPL])):
        tags={'PIXTABLE_OBJECT':[ 'PIXTABLE_OBJECT_%04d-%02d.fits'%(j+1,n+1) for n in range(24)]}
        tags['STD_RESPONSE']=['%s/STD/STD_RESPONSE_%s.fits'%(dir_out,time_std[idx_std])]
        tags['STD_TELLURIC']=['%s/STD/STD_TELLURIXS_%s.fits'%(dir_out,time_std[idx_std])]
        #tags['LSF_PROFILE']=['%s/WAVE/LSF_PROFILE_%s.fits'%(dir_out,time_wave[idx_wave])]
        
        muse.create_sof_scipost('scipost_temp.sof',tags,astrometry='',time=time[i].split('T')[0])
        muse.compute_scipost('scipost_temp.sof',cores,save=save,skymethod='none',darcheck='none',skymodel_frac=0.00)
        
        muse.subtract_sky('PIXTABLE_REDUCED_0001.fits','%s/%s/SKY_SPECTRUM_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno[j]),'%s/%s/PIXTABLE_REDUCED_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno[j]))
        if muse.nocache:
            os.system('nocache cp IMAGE_FOV_0001.fits %s/%s/IMAGE_FOV_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno[j])) 
        else:
            os.system('cp IMAGE_FOV_0001.fits %s/%s/IMAGE_FOV_%s_%02d.fits'%(dir_out,name.replace(' ',''),times[i],expno[j])) 
        #break
        os.system('rm IMAGE_FOV_*.fits')
        os.system('rm DATACUBE_FINAL.fits')
        m+=1
    os.system('rm PIXTABLE_*.fits')   
    os.chdir(cwd)

     





     
