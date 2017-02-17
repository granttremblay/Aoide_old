import os

objects=['NGC5248']

cwd = os.getcwd()
for q in objects:
      print q
      os.chdir('reduced/%s'%(q))
      os.system('nocache python ../../create_PCA_sky.py DATACUBE_SKY.fits PCA_sky.fits SKY_MASK.fits -f 25 -s 20000 -p 24')
      os.system('nocache python ../../subtract_PCA_sky.py DATACUBE_FINAL.fits DATACUBE_FINAL_clean.fits PCA_sky.fits -f 25 -c 50 -p 24')
      os.chdir(cwd)
