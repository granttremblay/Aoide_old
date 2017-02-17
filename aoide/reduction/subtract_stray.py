import pyfits
import numpy
from scipy import ndimage


def remove_stray_light_raw(input_file,output_file):
    hdu = pyfits.open(input_file)
    hdu_out=[None]*25
    hdu_out[0] = pyfits.PrimaryHDU(hdu[0].data)
    hdu_out[0].header = hdu[0].header


    for i in range(1,25):
        img = hdu[i].data
        hdr = hdu[i].header
        print    img.shape
        bias1=numpy.median(img[22:2119,1:31].flatten())
        bias2=numpy.median(img[2152:4204,1:31].flatten())
        bias3=numpy.median(img[22:2119,4192:-1].flatten())
        bias4=numpy.median(img[2120:4204,4192:-1].flatten())
        print bias1,bias2,bias3,bias4
        stray1=numpy.median(ndimage.filters.gaussian_filter1d(img[32:2087,32:39]-bias1,5,axis=0),1)
        stray2=numpy.median(ndimage.filters.gaussian_filter1d(img[2151:4207,32:39]-bias2,5,axis=0),1)
        stray3=numpy.median(ndimage.filters.gaussian_filter1d(img[32:2087,4183:4190]-bias3,5,axis=0),1)
        stray4=numpy.median(ndimage.filters.gaussian_filter1d(img[2151:4207,4183:4190]-bias4,5,axis=0),1)
        img[32:2087,32:2079] = img[32:2087,32:2079]-stray1[:,numpy.newaxis]
        img[2151:4207,32:2079] = img[2151:4207,32:2079]-stray2[:,numpy.newaxis]
        img[32:2087,2144:4191] = img[32:2087,2144:4191]-stray3[:,numpy.newaxis]
        img[2151:4207,2144:4191] = img[2151:4207,2144:4191]-stray4[:,numpy.newaxis]
    
        hdu_out[i] = pyfits.ImageHDU(img)
        hdu_out[i].header=hdu[i].header
    hdu =pyfits.HDUList(hdu_out)
    hdu.writeto(output_file,clobber=True)


for i in range(7,12):
     remove_stray_light_raw('../MUSE_WFM_FLAT039_%04d.fits'%(i),'MUSE_WFM_FLAT039_%04d_stray.fits'%(i))
