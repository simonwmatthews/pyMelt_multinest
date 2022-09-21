#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pyMelt as m
import pyMelt_multinest as mi

# Define the lithology objects to be used
lz = m.lithologies.matthews.klb1()
px = m.lithologies.matthews.kg1()
hz = m.lithologies.shorttle.harzburgite()

# These are the parameters the model should match. First item is the parameter value, the second is its 1 s.d. uncertainty.
data = {'tc':[5.74,0.27],
       'Tcrys':[1280,20],
       'Fpx':[0.175,0.1]}

# These are the parameters you don't want the inversion to find. For MORB P_lith and P_cryst must be set to zero to
# use the calculated crustal thickness.
knowns = {
       'DeltaS':300,
       'P_lith':0.0,
       'P_cryst':0.0,
       }

# These are the parameters for the inversion to find. The first item in the list is the type of prior, in this case
# all the priors are uniform. The next item is another list. For the uniform prior this list is the range of values
# the inversion should consider. If it were a normally distributed prior, use ['norm',[val,1sd]].
unknowns = {
       'Tp':['uni',[1250,1600]],
       'F_px':['uni',[0.0,1.0]],
       'F_hz':['uni',[0.0,0.99]]
       }

# Create the inversion object
inv = mi.inversion([lz,px,hz],data,knowns,unknowns,name='MORB')

# Run the inversion
inv.run_multinest()
