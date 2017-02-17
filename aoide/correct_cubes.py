import muse


objects=['NGC3351']
#objects=['HE0045-2145']
star_center={'NGC3351':[]}
A_V={'NGC3351':0.092}

offset={'NGC3351':(0,0)}
binning={'NGC3351':False}

for name in objects:
    print name
    muse.correct_cube2('reduced/%s/DATACUBE_FINAL_clean.fits'%(name),'final/%s.final.fits'%(name),A_V[name],star_center[name],offset[name],binning[name])
