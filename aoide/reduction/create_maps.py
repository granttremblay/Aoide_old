import pyfits
import numpy

hdu = pyfits.open('/home/grant/Projects/muse3c346/reduction/11a_science_modelsky_scipost/IMAGE_FOV_0001.fits')
hdr = hdu[1].header
hdu.close()

hdu = pyfits.open('3C346.eline_table.fits')
tab = hdu[1].data

x_cor = tab.field('x_cor')
y_cor = tab.field('y_cor')

Ha_flux = tab.field('Halpha_flux')
Ha_eflux = tab.field('Halpha_flux_err')
Ha_vel = tab.field('Halpha_vel')
Ha_evel = tab.field('Halpha_vel_err')
Ha_fwhm = tab.field('Halpha_fwhm')
Ha_efwhm = tab.field('Halpha_fwhm_err')


select = (Ha_flux/Ha_eflux>3.0) & (Ha_evel<10.0) & (Ha_efwhm<40.0)

Ha_flux_map = numpy.zeros((numpy.max(y_cor),numpy.max(x_cor)))
Ha_vel_map = numpy.zeros((numpy.max(y_cor),numpy.max(x_cor)))
Ha_fwhm_map = numpy.zeros((numpy.max(y_cor),numpy.max(x_cor)))

Ha_flux_map[y_cor[select],x_cor[select]] = Ha_flux[select]
Ha_vel_map[y_cor[select],x_cor[select]] = Ha_vel[select]
Ha_fwhm_map[y_cor[select],x_cor[select]] = Ha_fwhm[select]

hdu = pyfits.PrimaryHDU(Ha_flux_map)
hdu.header = hdr
hdu.writeto('Ha_flux_map.fits',clobber=True,output_verify='fix')

hdu = pyfits.PrimaryHDU(Ha_fwhm_map)
hdu.header = hdr
hdu.writeto('Ha_fwhm_map.fits',clobber=True,output_verify='fix')

hdu = pyfits.PrimaryHDU(Ha_vel_map)
hdu.header = hdr
hdu.writeto('Ha_vel_map.fits',clobber=True,output_verify='fix')