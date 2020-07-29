#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyMelt as m
import numpy as np
from pymultinest.run import run
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.interpolate import interp1d
import json
import os
import pandas as pd

class inversion:
    """
    Inversion object. Contains all the methods required by pyMultiNest to run.

    The parameters that must be defined collectively between knowns and unknowns are:
    Tp (mantle Tp in degC), DeltaS (Entropy change on melting), P_lith (pressure at base of lithosphere in GPa),
    P_cryst (pressure at which crystallisation takes place in GPa), F_px (proportion of pyroxenite in mantle),
    F_hz (proportion of harzburgite in the mantle).

    When buoyancy calculations are being used, the ambient Tp, F_px, and F_hz must also be included in either
    knowns or unknowns. If the data to be matched includes the melt flux (Q), knowns or unknowns must have r
    (the plume conduit radius in m), and mu (the plume viscosity in Pa s-1).

    Parameters
    ----------
    lithologies     list
        List of the pyMelt lithology objects in the order lherzolite, pyroxenite, harzburgite.
    data    dict
        Dictionary of the parameters and values for the inversion to match. Can include:
        'tc' (crustal thickness), 'Tcrys' (crystallisation temperature in degC), 'Fpx' (fraction
        of aggregate melts derived from pyroxenite), 'Q' (melt flux in m3s-1).
    knowns  dict
        Parameters required by the melting model that are to be set as fixed. See above. Keys are
        parameter names, values are the parameter values.
    unknowns dict
        Parameters to be found by the inversion. Keys are the parameter name. Dictionary values are
        a list of structure [PRIORTYPE,[low bound or mean, high bound or 1sd]. The prior type can be
        uniform ('uni'), log-uniform ('log-uni'), normal ('norm'), or log-normal ('lognorm').
    SpreadingCenter bool
        Model as a spreading center, or not? Default is True.
    TcrysShallow    bool
        Use the shallow Tcrys endmember, as opposed to the deep endmember. Default is True
    bouyancy    bool
        Require the solutions to be buoyant with respect to the defined ambient mantle properties.
        Default is False
    buoyancyPx  str
        Which column in the DensityFile to use for calculaitng pyroxenite density in buoyancy
        calculations. Default is 'kg1'.
    resume  bool
        Whether to resume an incomplete run of MultiNest, if files are present. Default is False.
    DensityFile str
        The file name for the csv file containing lithology densities as a function of temperature.
        Default is the file included in the repository.
    livepoints  int
        The number of livepoints that MultiNest should use.
    name    str
        The name to call the inversion, and the name of the folder to store the results in.

    """
    def __init__(self,lithologies,data,knowns,unknowns,SpreadingCenter=True,
                 TcrysShallow=True,buoyancy=False,buoyancyPx='kg1',resume=False,
                 DensityFile='LithDensity_80kbar.csv',livepoints=400,name='default'):

        self.livepoints = livepoints
        self.name = name
        self.unknowns = unknowns
        self.knowns = knowns
        self.parameters = len(self.unknowns)
        self.data = data
        self.lithologies = lithologies
        self.lithology_names = ['lz','px','hz']
        self.SpreadingCenter = SpreadingCenter
        self.TcrysShallow = TcrysShallow
        self.buoyancy = buoyancy
        self.buoyancyPx = buoyancyPx
        self.DensityFile = DensityFile
        self.resume = resume

        if self.buoyancy == False:
            self.var_list = ['Tp','DeltaS','P_lith','P_cryst','F_px','F_hz']
        else:
            if 'Q' in self.data.keys():
                self.var_list = ['Tp','DeltaS','P_lith','P_cryst','F_px','F_hz',
                             'ambientTp','ambientPx','ambientHz','r','mu']
            else:
                self.var_list = ['Tp','DeltaS','P_lith','P_cryst','F_px','F_hz',
                                 'ambientTp','ambientPx','ambientHz']
            self.rho = pd.read_csv(DensityFile,sep='\t')
            self.rhoLz = interp1d(self.rho.Tp,self.rho.klb1,fill_value='extrapolate')
            self.rhoHz = interp1d(self.rho.Tp,self.rho.hz,fill_value='extrapolate')
            self.rhoPx = interp1d(self.rho.Tp,self.rho[self.buoyancyPx],fill_value='extrapolate')


        self.unknowns_list = list()
        for variable in self.var_list:
            if variable in self.unknowns == False and variable in self.knowns == False:
                print('Please provide a value for '+variable+'.')
            if variable in list(self.unknowns.keys()):
                self.unknowns_list.append(variable)


    def prior(self,cube,ndim,nparams):
        i = 0
        for var in self.var_list:
            if var in self.unknowns:
                if self.unknowns[var][0] == 'uni':
                    cube[i] = self.unknowns[var][1][0] + cube[i] * (self.unknowns[var][1][1]-
                              self.unknowns[var][1][0])
                if self.unknowns[var][0] == 'norm':
                    g = norm(self.unknowns[var][1][0],self.unknowns[var][1][1])
                    cube[i] = g.ppf(cube[i])
                if self.unknowns[var][0] == 'loguni':
                    cube[i]= 10**(cube[i]*np.log10(self.unknowns[var][1][1]) + (1-cube[i])*np.log10(self.unknowns[var][1][0]))
                if self.unknowns[var][0] == 'lognorm':
                    g = lognorm(self.unknowns[var][1][0],self.unknowns[var][1][1])
                    cube[i] = g.ppf(cube[i])

                i = i+1
        return cube

    def loglike(self,cube,ndim,nparams):
        # List of values to pass to the forward model function
        x = list()
        # Position in n-dimensional cube
        i = 0
        for var in self.var_list:
            if var in self.unknowns:
                x.append(cube[i])
                i = i+1
            else:
                x.append(self.knowns[var])


        for i in range(len(self.lithologies)):
            self.lithologies[i].DeltaS = x[1]

        run_model = True
        if x[4] + x[5] > 1.0:
            run_model = False
            likelihood = -1e10*np.exp(1+x[4]+x[5])

        if run_model == True:
            proportions = [(1.0-x[5]-x[4]),x[4],x[5]]

            mantle = m.mantle(self.lithologies,proportions,self.lithology_names)
            if self.SpreadingCenter == False:
                results = mantle.AdiabaticMelt_1D(x[0],Pend=x[2],Pstart=10,ReportSSS=False)
            else:
                results = mantle.AdiabaticMelt_1D(x[0],Pstart=10,ReportSSS=False)
            if self.SpreadingCenter == True:
                results.integrate_tri()

            likelihood = 0

            if self.buoyancy == True or 'Q' in self.data.keys():
                ambient_rho = ((1.0-x[7]-x[8])*self.rhoLz(x[6])+
                                   x[7]*self.rhoPx(x[6])+
                                   x[8]*self.rhoHz(x[6])
                                    )

                model_rho = ((1.0-x[4]-x[5])*self.rhoLz(x[0])+
                                   x[4]*self.rhoPx(x[0])+
                                   x[5]*self.rhoHz(x[0])
                                    )

                buoyancy = ambient_rho - model_rho

            if 'Q' in self.data.keys():
                Q = np.pi/8 * (buoyancy * 9.81 * x[9]**4)/x[10] * results.F_total.max()
                likelihood = (likelihood + (-0.5*(np.log(2*np.pi*self.data['Q'][1]**2))-
                            (self.data['Q'][0]-Q)**2/(2*self.data['Q'][1]**2)))

            if self.buoyancy == True:
                if buoyancy < 0:
                    likelihood = likelihood - np.exp(-buoyancy) - 1


            if 'Tcrys' in self.data.keys():
                if self.SpreadingCenter == True:
                    TcrysMin, TcrysMax = results.MeltCrystallisationT(ShallowMeltP=False,MeltStorageP=False)
                else:
                    TcrysMin, TcrysMax = results.MeltCrystallisationT(ShallowMeltP=x[2],MeltStorageP=x[3])

                if self.TcrysShallow == True:
                    likelihood = (likelihood + (-0.5*(np.log(2*np.pi*self.data['Tcrys'][1]**2))-
                            (self.data['Tcrys'][0]-TcrysMin)**2/(2*self.data['Tcrys'][1]**2)))
                else:
                    likelihood = (likelihood + (-0.5*(np.log(2*np.pi*self.data['Tcrys'][1]**2))-
                            (self.data['Tcrys'][0]-TcrysMax)**2/(2*self.data['Tcrys'][1]**2)))

            if 'tc' in self.data.keys():
                likelihood = (likelihood + (-0.5*(np.log(2*np.pi*self.data['tc'][1]**2))-
                            (self.data['tc'][0]-results.tc)**2/(2*self.data['tc'][1]**2)))

            if 'Fpx' in self.data.keys():
                if self.SpreadingCenter == True:
                    likelihood = (likelihood + (-0.5*(np.log(2*np.pi*self.data['Fpx'][1]**2))-
                            (self.data['Fpx'][0]-results.tc_lithology_contributions[1])**2/(2*self.data['Fpx'][1]**2)))

        return likelihood

    def run_multinest(self,verbose=True):
        if not os.path.exists(self.name): os.mkdir(self.name)
        result = run(LogLikelihood=self.loglike,Prior=self.prior,n_dims=self.parameters,
                        outputfiles_basename=self.name+'/',resume=self.resume,
                        verbose=verbose,n_live_points=self.livepoints)
        json.dump(self.unknowns_list, open(self.name+'/params.json','w'))

        return result
