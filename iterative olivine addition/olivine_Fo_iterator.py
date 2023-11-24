# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 11:31:54 2021

@author: kevwo
"""

import numpy as np
from scipy.optimize import fsolve


def D_olLiqMg(P, T, H2O, totAlkalis):
    """
    Determine partitioning coefficients of Mg between melt and olivine (Putirka, 2007).

    Parameters
    ----------
    P : float
        Pressure (GPa).
    T : float
        Temperature (degC).
    H2O : float
        Water content of host magma (wt%).
    totAlkalis : float
        Total alkali content (Na2O + K2O) of host magma (wt%).

    Returns
    -------
    D_olLiqMg : float
        float, olivine-melt Mg partition coefficient

    """
    lnD_olLiqMg_1 = -2.158 + 55.09 * (P/T)
    lnD_olLiqMg_2 = -6.213e-2 * H2O
    lnD_olLiqMg_3 = (4430/T) + 5.115e-2 * totAlkalis
    lnD_olLiqMg = lnD_olLiqMg_1 + lnD_olLiqMg_2 + lnD_olLiqMg_3
    return np.exp(lnD_olLiqMg)


def D_olLiqFe(P, T, H2O, totAlkalis, SiO2):
    """
    Determine partitioning coefficients of Fe between melt and olivine (Putirka, 2007).

    Parameters
    ----------
    P : float
        Pressure (GPa).
    T : float
        Temperature (degC).
    H2O : float
        Water content of host magma (wt%).
    totAlkalis : float
        Total alkali content (Na2O + K2O) of host magma (wt%).
    SiO2 : float
        Silica content of host magma (wt%).

    Returns
    -------
    D_olLiqFe : float
        float, olivine-melt Fe partition coefficient

    """
    lnD_olLiqFe_1 = -3.300 + 47.57 * (P/T)
    lnD_olLiqFe_2 = -5.192e-2 * H2O
    lnD_olLiqFe_3 = (3344/T) + 5.595e-2 * totAlkalis
    lnD_olLiqFe_4 = 1.633e-2 * SiO2
    lnD_olLiqFe = lnD_olLiqFe_1 + lnD_olLiqFe_2 + lnD_olLiqFe_3 + lnD_olLiqFe_4
    return np.exp(lnD_olLiqFe)


def olivineStoich(Fo, D_olLiqMg, D_olLiqFe):
    """
    Function to convert Fo to Fe and Mg cation fractions in olivine and liquid via
    olivine stoichiometry.

    Parameters
    ----------
    Fo : float
        Olivine Fo (mol%).
    D_olLiqMg : float
        Partition function for Mg between olivine and liquid.
    D_olLiqFe : float
        Partition function for Fe between olivine and liquid.

    Returns
    -------
    list of proportions
        list, as follows: molar proportion of: Fe in olivine, Mg in olivine, Fe in liquid,
        Mg in liquid.

    """
    X_olMg = 0.667 * Fo / 100
    X_olFe = 0.667 - X_olMg
    X_liqFe = X_olFe / D_olLiqFe
    X_liqMg = X_olMg / D_olLiqMg
    return [X_olFe, X_olMg, X_liqFe, X_liqMg]


def FindXMgFe(Fo, D_olLiqMg, D_olLiqFe, KD=0.3):
    """
    Function to convert Fo to Fe and Mg cation fractions via the olivine-liquid Fe-Mg distribution
    coefficient and olivine stoichiometry.

    Parameters
    ----------
    Fo : float
        Olivine forsterite (mol%).
    D_olliq_Mg : float
        Partition function for Mg between olivine and liquid.
    D_olliq_Fe : float
        Partition function for Fe between olivine and liquid.
    KD : TYPE, optional
        Olivine-liquid Fe-Mg distribution coefficient. The default is 0.3.

    Returns
    -------
    list of proportions
        list, as follows: molar proportion of: Fe in liquid, Mg in liquid

    """
    MgFeLiq = KD / ((100/Fo) - 1)
    X_liqFe = 0.667 / (MgFeLiq * D_olLiqMg + D_olLiqFe)
    X_liqMg = MgFeLiq * X_liqFe
    return [X_liqFe, X_liqMg]


def iteration(Fo_0, T_crys_0, P, fraction=0.0001, H2O=0.0, totAlkalis=5.0, SiO2=48.0,
              targetFo=91.0):
    """
    Numerical iterator to find best-fitting T for given liquid Mg and Fe in equilibrium with a
    specific Fo, until a target Fo is achieved.

    Parameters
    ----------
    Fo_0 : float
        Olivine forsterite (mol%).
    T_crys_0 : float
        Crystallisation temperature (degC) of olivine with forsterite of Fo_0.
    P : float
        Pressure at which to perform iteration (GPa).
    fraction : float
        Fraction of olivine to add at every iteration step. The default is 0.0001.
    H2O : float, optional
        Water content of host magma (wt%). The default is 0.0.
    totAlkalis : float, optional
        Total alkali content (Na2O + K2O) of host magma (wt%). The default is 5.0.
    SiO2 : float, optional
        Silica content of host magma (wt%). The default is 48.0.
    targetFo : float, optional
        Target olivine forsterite (mol%) to terminate iteration at. The default is 91.0.

    Returns
    -------
    data : numpy.array
        numpy.array containing two rows. data[0] is olivine Fo, data[1] is crystallisation
        temperature for the corresponding olivine Fo.

    """

    Fo_list = [Fo_0]
    T_list = [T_crys_0]

    D_olLiqFeStart = D_olLiqFe(P, T_crys_0, H2O, totAlkalis, SiO2)
    D_olLiqMgStart = D_olLiqMg(P, T_crys_0, H2O, totAlkalis)

    Fo_crys = Fo_0
    T_crys = T_crys_0

    while Fo_crys < targetFo:
        cationFractions = FindXMgFe(Fo_crys, D_olLiqMgStart, D_olLiqFeStart)

        X_FeLiqNew = (cationFractions[0] + fraction * (1 - (Fo_crys/100))*2/3) / (1 + fraction)
        X_MgLiqNew = (cationFractions[1] + fraction * Fo_crys * 2/3/100) / (1 + fraction)

        def func(T):
            func = (X_MgLiqNew * D_olLiqMg(P, T, H2O, totAlkalis) +
                    X_FeLiqNew * D_olLiqFe(P, T, H2O, totAlkalis, SiO2))
            return 0.667 - func

        T_crys = fsolve(func, T_crys, full_output=0)[0]
        D_olLiqFeStart = D_olLiqFe(P, T_crys, H2O, totAlkalis, SiO2)
        D_olLiqMgStart = D_olLiqMg(P, T_crys, H2O, totAlkalis)

        Fo_crys = 100 * ((X_MgLiqNew/X_FeLiqNew/0.3)**-1 + 1)**-1

        Fo_list.append(Fo_crys)
        T_list.append(T_crys)

    data = np.zeros((2, len(Fo_list)))
    data[0] = Fo_list
    data[1] = T_list
    return data
