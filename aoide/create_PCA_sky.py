import argparse 
import muse

parser = argparse.ArgumentParser(description="Script to create a PCA library of spectra from an datacube given a mask")
parser.add_argument('input_cube',metavar='CUBE',type=str,nargs='?',help='Input FITS datacube from which sky spectra are selected')
parser.add_argument('output_pca',metavar='PCA',type=str,nargs='?',help='Output spectra of PCA components as FITS file')
parser.add_argument('sky_mask',metavar='MASK',type=str,nargs='?',help='FITS mask for the selection of spectra taken into acount')
parser.add_argument('-f','--filter',type=int,nargs='?',default=30, help='Size of median filter in wavelength direction to remove  continuum signal from the spectra')
parser.add_argument('-s','--spectra',type=int,nargs='?',default=20000, help='Maximum number of spectra considered for PCA analysis. The number of selected by the MASK is always the absolute maximum and can only be reduced.')
parser.add_argument('-p','--parallel',type=str,nargs='?',default='1', help='Number of cores used for computation.')

args = parser.parse_args()

muse.create_PCA_sky(args.input_cube,args.output_pca,args.sky_mask,args.filter,args.spectra,args.parallel)

