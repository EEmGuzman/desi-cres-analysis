#/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def getdata(fiberrow):
    """
    Get wavelength and flux data from frame file and return arrays of data.
    """
    with fits.open(str(sys.argv[1])) as inputfile:
        xdataname = inputfile['WAVELENGTH'].data
        yintermed = inputfile['FLUX'].data
        ydataname = yintermed[fiberrow]
        return (xdataname, ydataname)

with fits.open(str(sys.argv[1])) as infile:
    fibermap = infile[5].data
    fiberids = fibermap['FIBER']
    camera = infile[0].header['CAMERA']
    expid = infile[0].header['EXPID']
    flavor = infile[0].header['FLAVOR']


idealx, idealy = getdata(int(sys.argv[2]))
near1bpx, near1bpy = getdata(int(sys.argv[3]))
worstx, worsty = getdata(int(sys.argv[4]))
FibIDideal = fiberids[int(sys.argv[2])]
FibIDnearby = fiberids[int(sys.argv[3])]
FibIDworst = fiberids[int(sys.argv[4])]

plt.figure()
plt.plot(idealx, idealy, c='k' ,label='Ideal Spectrum (Fiber ID {})'.format(FibIDideal))
plt.plot(near1bpx, near1bpy, c='lawngreen' ,label='Nearby Spectrum 1bp (Fiber ID {})'.format(FibIDnearby))
plt.plot(worstx, worsty, c='b',label='Worst Spectrum (Fiber ID {})'.format(FibIDworst))
plt.legend()
plt.title('Camera: {}, EXPID: {}, Flavor: {}'.format(camera, expid, flavor))
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux')
plt.savefig(str(sys.argv[5]))
plt.close()
