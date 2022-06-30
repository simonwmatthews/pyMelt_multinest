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
import sys


class inversion:
    """
    Inversion object. Contains all the methods required by pyMultiNest to run.

    The parameters that must be defined collectively between knowns and unknowns are:
    Tp (mantle Tp in degC),
    DeltaS (Entropy change on melting),
    P_lith (pressure at base of lithosphere in GPa),
    P_cryst (pressure at which crystallisation takes place in GPa),
    F_px (proportion of pyroxenite in mantle),
    F_hz (proportion of harzburgite in the mantle).

    When buoyancy calculations are being used, the following must also be included in
    knowns and unknowns:
    ambientTp (ambient mantle Tp in degC),
    ambientPx (proportion of pyroxenite in ambient mantle)
    ambientHz (proportion of harzburgite in ambient mantle)

    When plume flux calculations are being used, the following must also be included in
    knowns and unknowns:
    r (plume radius in m)
    mu (plume viscosity in Pa s-1)

    When trace element calculations are being used, the following must also be included in
    knowns and
    unknowns:
    La_lz, Dy_lz, Yb_lz (concentrations of La, Dy, and Yb in mantle lherzolite)
    La_px, Dy_px, Yb_px (concentrations of La, Dy, and Yb in mantle pyroxenite
                         if not using MORBmelts)
    MORBmelts (melt fraction of a MORB batch melt if using MORBmelts)

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
        Parameters to be found by the inversion. Keys are the parameter name. Dictionary values
        are a list of structure [PRIORTYPE,[low bound or mean, high bound or 1sd]. The prior type
        can be uniform ('uni'), log-uniform ('log-uni'), normal ('norm'), or log-normal
        ('lognorm').
    DeltaP  float
        Fixed integration pressure step in GPa.
    SpreadingCentre bool
        Model as a spreading center, or not? Default is True.
    ContinentalRift bool
        Model as a magma-assisted continental rift, or not? Default is False.
    Passive    bool
        Model with passive upwelling, or not? Default is True.
    Traces     bool
        Model with trace element constraints or not? Default is False.
    MORBmelts  bool
        Model pyroxenite traces with melts of MORB or not? Default is False
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

    def __init__(self, lithologies, data, knowns, unknowns, DeltaP=0.004, SpreadingCentre=True,
                 ContinentalRift=False, Passive=True, Traces=False, MORBmelts=False,
                 TcrysShallow=True, buoyancy=False, buoyancyPx='kg1', resume=False,
                 DensityFile='LithDensity_80kbar.csv', livepoints=400, name='default'):

        self.livepoints = livepoints
        self.name = name
        self.unknowns = unknowns
        self.knowns = knowns
        self.parameters = len(self.unknowns)
        self.data = data
        self.lithologies = lithologies
        self.lithology_names = ['lz', 'px', 'hz']
        self.DeltaP = DeltaP
        self.SpreadingCentre = SpreadingCentre
        self.ContinentalRift = ContinentalRift
        self.Passive = Passive
        self.Traces = Traces
        self.MORBmelts = MORBmelts
        self.TcrysShallow = TcrysShallow
        self.buoyancy = buoyancy
        self.buoyancyPx = buoyancyPx
        self.DensityFile = DensityFile
        self.resume = resume

        if (self.SpreadingCentre and self.ContinentalRift) is True:
            sys.exit('The system cannot simultaneously be an oceanic spreading center '
                     + 'and a continental rift.')

        if self.Traces is False and self.MORBmelts is True:
            sys.exit('Cannot model pyroxenite traces without considering all traces.')

        variables = ['Tp', 'DeltaS', 'P_lith', 'P_cryst', 'F_px', 'F_hz']

        if self.buoyancy is True:
            variables.extend(['ambientTp', 'ambientPx', 'ambientHz'])

            if 'Qv' in self.data.keys() or 'Qb' in self.data.keys() or 'Qm' in self.data.keys():
                variables.extend(['r', 'mu'])

            self.rho = pd.read_csv(DensityFile)
            self.rhoLz = interp1d(self.rho['Tp'],
                                  self.rho['klb1'],
                                  fill_value='extrapolate')
            self.rhoHz = interp1d(self.rho['Tp'],
                                  self.rho['hz'],
                                  fill_value='extrapolate')
            self.rhoPx = interp1d(self.rho['Tp'],
                                  self.rho[self.buoyancyPx],
                                  fill_value='extrapolate')

        if self.Passive is False:
            variables.extend(['lambda', 'amplitude'])

        if self.Traces is True:
            if ('La_Yb' and 'Dy_Yb') in self.data.keys():
                variables.extend(['La_lz', 'Dy_lz', 'Yb_lz'])
                if self.MORBmelts is True:
                    variables.extend(['MORBmelts'])
                else:
                    variables.extend(['La_px', 'Dy_px', 'Yb_px'])
            else:
                sys.exit('Please enter \'La/Yb\', and \'Dy/Yb\' into the data and the '
                         + 'prior source compositions.'
                         )

        self.var_list = variables
        self.unknowns_list = list()
        for variable in self.var_list:
            if variable in self.unknowns is False and variable in self.knowns is False:
                print('Please provide a value for '+variable+'.')

            if variable in list(self.unknowns.keys()):
                self.unknowns_list.append(variable)

    def prior(self, cube, ndim, nparams):
        i = 0
        for var in self.var_list:
            if var in self.unknowns:
                if self.unknowns[var][0] == 'uni':
                    cube[i] = self.unknowns[var][1][0] + cube[i] * (self.unknowns[var][1][1] -
                                                                    self.unknowns[var][1][0])

                if self.unknowns[var][0] == 'norm':
                    g = norm(self.unknowns[var][1][0], self.unknowns[var][1][1])
                    cube[i] = g.ppf(cube[i])

                if self.unknowns[var][0] == 'loguni':
                    cube[i] = 10**(cube[i] * np.log10(self.unknowns[var][1][1]) +
                                   (1 - cube[i]) * np.log10(self.unknowns[var][1][0]))

                if self.unknowns[var][0] == 'lognorm':
                    g = lognorm(self.unknowns[var][1][0], self.unknowns[var][1][1])
                    cube[i] = g.ppf(cube[i])

                i = i+1
        return cube

    def loglike(self, cube, ndim, nparams):
        x = list()
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
        if x[self.var_list.index('F_px')] + x[self.var_list.index('F_hz')] > 1.0:
            run_model = False
            likelihood = -1e10 * np.exp(1 + x[self.var_list.index('F_px')] +
                                        x[self.var_list.index('F_hz')])

        elif x[self.var_list.index('P_lith')] < x[self.var_list.index('P_cryst')]:
            run_model = False
            likelihood = -1e10 * np.exp(1 + x[self.var_list.index('P_lith')] +
                                        x[self.var_list.index('P_cryst')])

        else:
            proportions = [(1.0 - x[self.var_list.index('F_hz')] - x[self.var_list.index('F_px')]),
                           x[self.var_list.index('F_px')], x[self.var_list.index('F_hz')]]
            mantle = m.mantle(self.lithologies, proportions, self.lithology_names)
            SolidusPressures = mantle.solidusIntersection(x[self.var_list.index('Tp')])

            SolidusPressureCheck = np.isnan(SolidusPressures).all()
            if SolidusPressureCheck:
                run_model = False
                likelihood = -1e12
            else:
                self.SolidusIntersectP = np.nanmax(SolidusPressures)
                if self.SolidusIntersectP < x[self.var_list.index('P_lith')]:
                    run_model = False
                    likelihood = -1e10 * np.exp(1 + x[self.var_list.index('P_lith')] -
                                                self.SolidusIntersectP)
                elif abs(self.SolidusIntersectP - x[self.var_list.index('P_lith')]) < self.DeltaP:
                    run_model = False
                    likelihood = -1e10 * np.exp(1 + x[self.var_list.index('P_lith')] -
                                                self.SolidusIntersectP)

        if run_model is True:
            likelihood = 0

            if self.SpreadingCentre is True:
                results = mantle.adiabaticMelt(Tp=x[self.var_list.index('Tp')],
                                               Pstart=self.SolidusIntersectP,
                                               dP=(-1)*self.DeltaP,
                                               ReportSSS=False)
            else:
                results = mantle.adiabaticMelt(Tp=x[self.var_list.index('Tp')],
                                               Pstart=self.SolidusIntersectP,
                                               Pend=x[self.var_list.index('P_lith')],
                                               dP=(-1)*self.DeltaP,
                                               ReportSSS=False)

            if self.Traces is True:
                if self.MORBmelts is True:
                    La_px = 0.192 / (x[self.var_list.index('MORBmelts')] * (1 - 0.01) + 0.01)
                    Dy_px = 0.505 / (x[self.var_list.index('MORBmelts')] * (1 - 0.079) + 0.079)
                    Yb_px = 0.365 / (x[self.var_list.index('MORBmelts')] * (1 - 0.115) + 0.115)
                    results.calculateChemistry(
                        elements={'lz': {'La': x[self.var_list.index('La_lz')],
                                         'Dy': x[self.var_list.index('Dy_lz')],
                                         'Yb': x[self.var_list.index('Yb_lz')]},
                                  'px': {'La': La_px,
                                         'Dy': Dy_px,
                                         'Yb': Yb_px},
                                  'hz': m.chemistry.workman05_dmm},
                        cpxExhaustion={'lz': 0.18,
                                       'px': 0.70,
                                       'hz': 0.10},
                        garnetOut={
                            'lz': m.chemistry.mineralTransition_linear(
                                {'gradient': 1/666.7, 'intercept': 400/666.7}),
                            'px': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 1.5}),
                            'hz': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 1.5}),
                                     },
                        spinelIn={
                            'lz': m.chemistry.mineralTransition_linear(
                                {'gradient': 1/666.7, 'intercept': 533/666.7}),
                            'px': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 2.5}),
                            'hz': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 2.5}),
                                     },
                        mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                            'px': m.chemistry.kg1_MineralProportions,
                                            'hz': m.chemistry.klb1_MineralProportions}
                        )

                else:
                    results.calculateChemistry(
                        elements={'lz': {'La': x[self.var_list.index('La_lz')],
                                         'Dy': x[self.var_list.index('Dy_lz')],
                                         'Yb': x[self.var_list.index('Yb_lz')]},
                                  'px': {'La': x[self.var_list.index('La_px')],
                                         'Dy': x[self.var_list.index('Dy_px')],
                                         'Yb': x[self.var_list.index('Yb_px')]},
                                  'hz': m.chemistry.workman05_dmm},
                        cpxExhaustion={'lz': 0.18,
                                       'px': 0.70,
                                       'hz': 0.10},
                        garnetOut={
                            'lz': m.chemistry.mineralTransition_linear(
                                {'gradient': 1/666.7, 'intercept': 400/666.7}),
                            'px': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 1.5}),
                            'hz': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 1.5}),
                                     },
                        spinelIn={
                            'lz': m.chemistry.mineralTransition_linear(
                                {'gradient': 1/666.7, 'intercept': 533/666.7}),
                            'px': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 2.5}),
                            'hz': m.chemistry.mineralTransition_isobaric(
                                {'transition_pressure': 2.5}),
                                     },
                        mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                            'px': m.chemistry.kg1_MineralProportions,
                                            'hz': m.chemistry.klb1_MineralProportions}
                        )

            if self.Passive is True:
                if self.SpreadingCentre is True:
                    geosetting = m.geosettings.spreadingCentre(results)
                elif self.ContinentalRift is True:
                    geosetting = m.geosettings.spreadingCentre(
                        results, P_lithosphere=x[self.var_list.index('P_lith')],
                        extract_melt=True)
                else:
                    geosetting = m.geosettings.intraPlate(
                        results, P_lithosphere=x[self.var_list.index('P_lith')])
            else:
                if self.SpreadingCentre is True:
                    geosetting = m.geosettings.spreadingCentre(
                        results, weightingFunction=m.geosettings.weighting_expdecay,
                        weighting_wavelength=x[self.var_list.index('lambda')],
                        weighting_amplitude=x[self.var_list.index('amplitude')])
                elif self.ContinentalRift is True:
                    geosetting = m.geosettings.spreadingCentre(
                        results, P_lithosphere=x[self.var_list.index('P_lith')],
                        extract_melt=True,
                        weightingFunction=m.geosettings.weighting_expdecay,
                        weighting_wavelength=x[self.var_list.index('lambda')],
                        weighting_amplitude=x[self.var_list.index('amplitude')])
                else:
                    geosetting = m.geosettings.intraPlate(
                        results, P_lithosphere=x[self.var_list.index('P_lith')],
                        weightingFunction=m.geosettings.weighting_expdecay,
                        weighting_wavelength=x[self.var_list.index('lambda')],
                        weighting_amplitude=x[self.var_list.index('amplitude')])

            if self.buoyancy is True:
                ambient_rho = ((1.0 - x[self.var_list.index('ambientPx')] -
                                x[self.var_list.index('ambientHz')]) *
                               self.rhoLz(x[self.var_list.index('ambientTp')]) +
                               x[self.var_list.index('ambientPx')] *
                               self.rhoPx(x[self.var_list.index('ambientTp')]) +
                               x[self.var_list.index('ambientHz')] *
                               self.rhoHz(x[self.var_list.index('ambientTp')])
                               )
                model_rho = ((1.0 - x[self.var_list.index('F_px')] -
                              x[self.var_list.index('F_hz')]) *
                             self.rhoLz(x[self.var_list.index('Tp')]) +
                             x[self.var_list.index('F_px')] *
                             self.rhoPx(x[self.var_list.index('Tp')]) +
                             x[self.var_list.index('F_hz')] *
                             self.rhoHz(x[self.var_list.index('Tp')]))
                buoyancy = ambient_rho - model_rho

            if 'Qv' in self.data.keys():
                Qv = np.pi/8 * (buoyancy * 9.81 *
                                x[self.var_list.index('r')]**4)/x[self.var_list.index('mu')]
                likelihood = (likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Qv'][1]**2)) -
                                            (self.data['Qv'][0] - Qv)**2 /
                                            (2*self.data['Qv'][1]**2)))
            elif 'Qb' in self.data.keys():
                Qv = np.pi/8 * (buoyancy * 9.81 *
                                x[self.var_list.index('r')]**4)/x[self.var_list.index('mu')]
                Qb = (Qv/1e3) * model_rho * 30e-6 * (x[0]-x[6])
                likelihood = (
                    likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Qb'][1]**2)) -
                                  (self.data['Qb'][0] - Qb)**2 / (2 * self.data['Qb'][1]**2)))
            elif 'Qm' in self.data.keys():
                Qv = np.pi/8 * (buoyancy * 9.81 *
                                x[self.var_list.index('r')]**4)/x[self.var_list.index('mu')]
                Qm = Qv * results.F_total.max()
                likelihood = (
                    likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Qm'][1]**2)) -
                                  (self.data['Qm'][0] - Qm)**2 / (2*self.data['Qm'][1]**2)))

            if self.buoyancy is True:
                if buoyancy < 0:
                    likelihood = likelihood - np.exp(-buoyancy) - 1

            if 'Tcrys' in self.data.keys():
                if self.SpreadingCentre is True:
                    TcrysMin, TcrysMax = geosetting.meltCrystallisationT(
                        ShallowMeltP=None,
                        MeltStorageP=None)
                else:
                    TcrysMin, TcrysMax = geosetting.meltCrystallisationT(
                        ShallowMeltP=x[self.var_list.index('P_lith')],
                        MeltStorageP=x[self.var_list.index('P_cryst')])

                if self.TcrysShallow is True:
                    likelihood = (
                        likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Tcrys'][1]**2)) -
                                      (self.data['Tcrys'][0] - TcrysMin)**2 /
                                      (2 * self.data['Tcrys'][1]**2)))
                else:
                    likelihood = (
                        likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Tcrys'][1]**2)) -
                                      (self.data['Tcrys'][0] - TcrysMax)**2 /
                                      (2 * self.data['Tcrys'][1]**2)))

            if 'tc' in self.data.keys():
                if (self.SpreadingCentre or self.ContinentalRift) is True:
                    likelihood = (
                        likelihood + (-0.5 * (np.log(2 * np.pi * self.data['tc'][1]**2)) -
                                      (self.data['tc'][0] - geosetting.tc)**2 /
                                      (2 * self.data['tc'][1]**2)))

            if 'Fpx' in self.data.keys():
                if (self.SpreadingCentre or self.ContinentalRift) is True:
                    likelihood = (
                        likelihood + (-0.5*(np.log(2 * np.pi*self.data['Fpx'][1]**2)) -
                                      (self.data['Fpx'][0] -
                                       geosetting.lithology_contributions[1])**2 /
                                      (2*self.data['Fpx'][1]**2)))

            if self.Traces is True:
                La_Yb_ratio = geosetting.chemistry['La'] / geosetting.chemistry['Yb']
                Dy_Yb_ratio = geosetting.chemistry['Dy'] / geosetting.chemistry['Yb']

                likelihood = (
                    likelihood + (-0.5 * (np.log(2 * np.pi * self.data['La_Yb'][1]**2)) -
                                  (self.data['La_Yb'][0] - La_Yb_ratio)**2 /
                                  (2 * self.data['La_Yb'][1]**2)))
                likelihood = (
                    likelihood + (-0.5 * (np.log(2 * np.pi * self.data['Dy_Yb'][1]**2)) -
                                  (self.data['Dy_Yb'][0] - Dy_Yb_ratio)**2 /
                                  (2 * self.data['Dy_Yb'][1]**2)))
        return likelihood

    def run_multinest(self, verbose=True):
        if not os.path.exists(self.name):
            os.makedirs(self.name, exist_ok=True)
        result = run(LogLikelihood=self.loglike, Prior=self.prior, n_dims=self.parameters,
                     outputfiles_basename=self.name+'/', resume=self.resume, verbose=verbose,
                     n_live_points=self.livepoints)
        json.dump(self.unknowns_list, open(self.name+'/params.json', 'w'))
        return result
