#!/usr/bin/env python3

import sys
import csv
import numpy as np
from astropy.io import fits

# Importing the [500, 3] WSIGMA array.
with fits.open(str(sys.argv[1])) as file_psf:
    wsigma_array = file_psf['WSIGMA'].data

# Importing frame.fibermap['FIBERS']
with fits.open(str(sys.argv[2])) as infile:
    fibermap = infile[5].data
    fiberids = fibermap['FIBER']

# Variables for plot
p0 = wsigma_array[:, 0:1].flatten()
p1 = wsigma_array[:, 1:2].flatten()
p2 = wsigma_array[:, 2:3].flatten()
arrlist = [p0, p1, p2]
# Calculating the medians by column. Coeffs to be used as 'Model'
medlegpolcoef = np.median(wsigma_array,axis=0)
print("The wsigma_medians are {} ".format(medlegpolcoef))

# Calcualting Standard Deviation (New Version)
wsigma_rms = np.sqrt(np.mean((wsigma_array - medlegpolcoef)**2,axis=0))
print("The wsigma_rms are {} ".format(wsigma_rms))

# Calculating the distance each value is away from the median.
diffp0 = p0 - medlegpolcoef[0]
diffp1 = p1 - medlegpolcoef[1]
diffp2 = p2 - medlegpolcoef[2]

# Checking how many of each parameter are outside of +- 2 RMS of the median.
toperror = np.array([medlegpolcoef[val] + 2*wsigma_rms[val] for val in [0,1,2]])
bottomerror = np.array([medlegpolcoef[val] - 2*wsigma_rms[val] for val in [0,1,2]])

rms1toperr = np.array([medlegpolcoef[val] + 1*wsigma_rms[val] for val in [0,1,2]])
rms1boterr = np.array([medlegpolcoef[val] - 1*wsigma_rms[val] for val in [0,1,2]])

badparamrnum0 = list(np.where(np.logical_or(p0>toperror[0], p0<bottomerror[0]))[0])
badparamrnum1 = list(np.where(np.logical_or(p1>toperror[1], p1<bottomerror[1]))[0])
badparamrnum2 = list(np.where(np.logical_or(p2>toperror[2], p2<bottomerror[2]))[0])
gparamrnum0 = list(np.where(np.logical_and(p0<rms1toperr[0], p0>rms1boterr[0]))[0])
gparamrnum1 = list(np.where(np.logical_and(p1<rms1toperr[1], p1>rms1boterr[1]))[0])
gparamrnum2 = list(np.where(np.logical_and(p2<rms1toperr[2], p2>rms1boterr[2]))[0])


nbadparam = np.array([len(badparamrnum0), len(badparamrnum1), len(badparamrnum2)])
print('The number of bad coefficients are {} '.format(nbadparam))
print('The rows with a bad param0 are {} '.format(badparamrnum0))
print('The rows with a bad param1 are {} '.format(badparamrnum1))
print('The rows with a bad param2 are {} '.format(badparamrnum2))

# Finding ideal fibers (row number) that have all params <1RMS
goodfibrnum = sorted(list(set(gparamrnum0) & set(gparamrnum1) & set(gparamrnum2)))
print('The Ideal Fiber Rows are {},\nThe number of Ideal Fibers are {}'.format(goodfibrnum,len(goodfibrnum)))

# Determining which fibers have 2 or more parameters >2*RMS away from median
rnum3bparam = sorted(list(set(badparamrnum0) & set(badparamrnum1) & set(badparamrnum2)))

rnum2bparams1 = set(badparamrnum0) & set(badparamrnum1)
rnum2bparams2 = set(badparamrnum0) & set(badparamrnum2)
rnum2bparams3 = set(badparamrnum1) & set(badparamrnum2)
rnum2bparamint = rnum2bparams1|rnum2bparams2|rnum2bparams3
rnum2bparam = sorted(list(rnum2bparamint.difference(set(rnum3bparam))))

#print('Row with 3 bad params {} '.format(rnum3bparam))
#print('Row with 2 bad params {} '.format(rnum2bparam))
bp2fibids = []
bp3fibids = []
if len(rnum3bparam) != 0:
    for value in range(len(rnum3bparam)):
        bp3fibids.append(fiberids[rnum3bparam[value]])
    print('Fiber IDs with 3 bad params {}'.format(sorted(bp3fibids)))
if len(rnum2bparam) != 0:
    for value in range(len(rnum2bparam)):
        bp2fibids.append(fiberids[rnum2bparam[value]])
    print('Fiber IDs with 2 bad params {}'.format(sorted(bp2fibids)))
print('The number of odd fibers are {}'.format(len(rnum2bparam + rnum3bparam)))

# Writing information about bad fibers to csv file
finfibid = []
findiffp0 = []
findiffp1 = []
findiffp2 = []
finfibid_3bp = []
findiffp0_3bp = []
findiffp1_3bp = []
findiffp2_3bp = []
if len(rnum2bparam) != 0:
    counter = 0
    for value in rnum2bparam:
        finfibid.append(bp2fibids[counter])
        findiffp0.append(diffp0[value])
        findiffp1.append(diffp1[value])
        findiffp2.append(diffp2[value])
        counter += 1
if len(rnum3bparam) != 0:
    counter = 0
    for value in rnum3bparam:
        finfibid_3bp.append(bp3fibids[counter])
        findiffp0_3bp.append(diffp0[value])
        findiffp1_3bp.append(diffp1[value])
        findiffp2_3bp.append(diffp2[value])
        counter += 1

bfibinfo = zip(rnum2bparam, finfibid, findiffp0, findiffp1, findiffp2)
bfibinfo_3bp = zip(rnum3bparam, finfibid_3bp, findiffp0_3bp, findiffp1_3bp, findiffp2_3bp)
with open(str(sys.argv[3]), 'w') as inputfile:
    dwrite = csv.writer(inputfile)
    dwrite.writerow(["RNum","Fiber_ID", 'LPC_p0_diff', 'LPC_p1_diff','LPC_p2_diff'])
    for row in bfibinfo:
        dwrite.writerow(row)
    dwrite.writerow(["3 Bad Coefficients"])
    for row in bfibinfo_3bp:
        dwrite.writerow(row)
