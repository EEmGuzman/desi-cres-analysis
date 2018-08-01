#!/usr/bin/env python3
#./residplot.py framefile idealfiberrow worstfiberrow

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
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
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

# Subplot 1
ax0 = plt.subplot(gs[0])
ax0.plot(idealx, idealy/1000, c='k',label='Ideal Spectrum (Fiber ID {})'.format(FibIDideal))
ax0.plot(worstx, worsty/1000, c='b',label='Worst Spectrum (Fiber ID {})'.format(FibIDworst))
ax0.legend()
ax0.set_ylabel('Flux ($10^{3}$)')
ax0.set_title('Camera: {}, EXPID: {}, Flavor: {}'.format(camera, expid, flavor))

# Subplot 2
ax1 = plt.subplot(gs[1], sharex = ax0)
ax1.scatter(worstx, worstydiff/1000, s=8, marker='o')
ax1.axhline(y=0, linestyle='--', color='red')
plt.setp(ax0.get_xticklabels(), visible=False)
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
ax1.set_ylabel('Residual')
plt.xlabel('Wavelength (Angstroms)')

plt.autoscale(enable=True, axis='both', tight=False)
plt.subplots_adjust(hspace=.0)
plt.savefig('residp-{}-{}{}.eps'.format(camera,str(expid),flavor))
plt.close(fig1)
