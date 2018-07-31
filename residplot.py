#!/usr/bin/env python3
#./residplot.py framefile idealfiberrow worstfiberrow outputfname

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def getdata(fiberrow):
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
worstx, worsty = getdata(int(sys.argv[3]))
FibIDideal = fiberids[int(sys.argv[2])]
FibIDworst = fiberids[int(sys.argv[3])]

worstydiff = worsty - idealy

fig1 = plt.figure(1)
frame1 = fig1.add_axes((.1,.3,.8,.6))
plt.plot(idealx, idealy/1000, c='k',label='Ideal Spectrum (Fiber ID {})'.format(FibIDideal))
plt.plot(worstx, worsty/1000, c='b',label='Worst Spectrum (Fiber ID {})'.format(FibIDworst))
plt.title('Camera: {}, EXPID: {}, Flavor: {}'.format(camera, expid, flavor))
plt.ylabel('Flux ($10^{3}$)')
plt.autoscale(enable=True, axis='both', tight=False)
plt.legend()
frame2 = fig1.add_axes((.1,.1,.8,.2))
plt.scatter(worstx, worstydiff/1000, s=8, marker='o')
plt.axhline(y=0, linestyle='--', color='red')
frame2.set_ylabel('Residual')
plt.xlabel('Wavelength (Angstroms)')
plt.savefig('residp-{}-{}{}.eps'.format(camera,str(expid),flavor))
plt.close(fig1)
