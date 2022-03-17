#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import timeit
import pyMelt as m
import pyMelt_multinest as mi
from datetime import datetime


# Enter the name of the model ran here:
run_name = 'BVCClusterMPITestPassivePxGtSp'

# This function provides a date for the model label.
now = datetime.now()
date = now.strftime('%Y%m%d')

# These are functions to keep track of how long the code has been running for:
start = timeit.default_timer()

# Define the lithology objects to be used
lz = m.lithologies.matthews.klb1()
px = m.lithologies.matthews.kg1()
hz = m.lithologies.shorttle.harzburgite()

# These are the parameters the model should match. First item is the parameter value, the second
# is its 1 s.d. uncertainty.
data = {
        'Tcrys': [1426, 26],
        'tc': [4.5, 1.5],
        'La_Yb': [9.70, 0.85],
        'Dy_Yb': [2.17, 0.07]
        }

# These are the parameters you don't want the inversion to find. For MORB P_lith and P_cryst must
# be set to zero to use the calculated crustal thickness.
knowns = {
        'DeltaS': 300
        }

# These are the parameters for the inversion to find. The first item in the list is the type of
# prior, in this case all the priors are uniform. The next item is another list. For the uniform
# prior this list is the range of values the inversion should consider. If it were a normally
# distributed prior, use ['norm',[val,1sd]].
unknowns = {
    'Tp': ['uni', [1250, 1750]],
    'F_px': ['uni', [0, 1]],
    'F_hz': ['uni', [0, 1]],
    'P_lith': ['norm', [60*28.0/1000, 20*28.0/1000]],
    'ambientTp': ['norm', [1364, 15/2]],
    'ambientPx': ['norm', [0.021, 0.012]],
    'ambientHz': ['norm', [0.42, 0.15]],
    'P_cryst': ['norm', [28*28.0/1000, 3*28.0/1000]],
    'La_lz': ['norm', [0.648, 0.648*0.1]],
    'Dy_lz': ['norm', [0.674, 0.674*0.1]],
    'Yb_lz': ['norm', [0.441, 0.441*0.1]],
    'La_px': ['norm', [2.701, 2.701*0.1]],
    'Dy_px': ['norm', [3.2925, 3.2925*0.1]],
    'Yb_px': ['norm', [1.9975, 1.9975*0.1]],
#    'lambda': ['uni', [0.0, 1.0]],
#    'amplitude': ['uni', [0.0, 1.0]]
            }

# This is the inversion object. Make sure that for each of the knowns and unknowns the correct
# combinations of arguments has been entered.
inv = mi.inversion(
    lithologies=[lz, px, hz],
    deltaP=0.004,
    data=data,
    knowns=knowns,
    unknowns=unknowns,
    SpreadingCentre=False,
    ContinentalRift=True,
    Passive=True,
    Traces=True,
    MORBmelts=False,
    TcrysShallow=True,
    buoyancy=True,
    buoyancyPx='kg1',
    resume=False,
    DensityFile='LithDensity_80kbar.csv',
    livepoints=400,
    name=date+'_'+run_name+'_'+str(data['Tcrys'][0])+'_'+str(data['Tcrys'][1])
    )

# Run the inversion
inv.run_multinest()

# Stop the timer at the end of the inversion
stop = timeit.default_timer()
print('Time in minutes: ', (stop - start)/60)
