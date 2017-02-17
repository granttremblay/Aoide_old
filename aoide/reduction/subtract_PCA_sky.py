import argparse 
import muse
import copy_reg
from types import *

parser = argparse.ArgumentParser(description="Script to subtract sky residuals from a datacube based on a PCA spectral library")
parser.add_argument('input_cube',metavar='CUBE_IN',type=str,nargs='?',help='Input FITS datacube from which sky spectra are selected')
parser.add_argument('output_cube',metavar='CUBE_OUT',type=str,nargs='?',help='Output FITS datacube with sky residuals subtracted')
parser.add_argument('pca_sky',metavar='PCA',type=str,nargs='?',help='Input FITS file with PCA components')
parser.add_argument('-c','--components',type=int,nargs='?',default=150, help='Number of PCA components to be used')
parser.add_argument('-f','--filter',type=int,nargs='?',default=30, help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction')
parser.add_argument('-p','--parallel',type=str,nargs='?',default='1', help='Number of cores used for computation.')

args = parser.parse_args()

muse.subtract_PCA_sky(args.input_cube,args.output_cube,args.pca_sky,args.components,args.filter,args.parallel)

def _pickle_method(method):
      func_name = method.im_func.__name__
      obj = method.im_self
      cls = method.im_class
      if func_name.startswith('__') and not func_name.endswith('__'):
	cls_name = cls.__name__.lstrip('_')
	if cls_name: func_name = '_' + cls_name + func_name
      return _unpickle_method, (func_name, obj, cls)
    
def _unpickle_method(func_name, obj, cls):
      for cls in cls.mro():
	try:
	  func = cls.__dict__[func_name]
	except KeyError:
	  pass
	else:
	  break
      return func.__get__(obj, cls)
copy_reg.pickle(MethodType,_pickle_method, _unpickle_method)
