## pyMelt_MultiNest (pyMelt v2.0)
pyMelt can be used in conjunction with the MultiNest algorithm (Feroz and Hobson, 2008; Feroz et al., 2009, 2013) via its Python frontend, pyMultiNest (Buchner et al., 2014). This permits the inversion of measured data (e.g., crystallisation temperature, crustal thickness, or rare-earth element concentration ratios) to obtain unknowns (e.g., potential temperature) via Bayesian inference. More details of the inversion methods are provided in Matthews et al. (2021).

For pyMelt_MultiNest to work, MultiNest and pyMultiNest must be installed, as well as ```pyMelt```, ```numpy```, ```scipy```, and ```pandas```. The user is directed to the [pyMultinest installation instructions](https://johannesbuchner.github.io/PyMultiNest/) for further guidance.

pyMelt_MultiNest is installed through placing the pyMelt_multinest.py file in the working directory: 


```python
import pyMelt as m
import pyMelt_multinest as mi
```

If operating on a Linux installation, the latest versions of pyMelt and pyMelt_MultiNest can be downloaded from the [code repository](https://github.com/simonwmatthews/pyMelt) and installed locally:

```shell
$ git clone https://github.com/simonwmatthews/pyMelt
$ git clone https://github.com/simonwmatthews/pyMelt_multinest
```

First define the lithology objects to be used in the inversion:


```python
lz = m.lithologies.matthews.klb1()
px = m.lithologies.matthews.kg1()
hz = m.lithologies.shorttle.harzburgite()
```

Next the data, knowns, and unknowns must be specified.

```data``` is a dictionary of the initial parameters that the model should aim to match. The data currently supported by pyMelt_MultiNest are:

* ```'Tcrys'```: melt liquidus/crystallisation temperature (&deg;C)
* ```'tc'```: crustal thickness (km)
* ```'Fpx'```: fraction of aggregate melts derived from pyroxenite
* ```'Qm'```, ```'Qb'```, or ```'Qv'```: mantle plume melt flux (m<sup>3</sup> s<sup>-1</sup>), buoyancy flux (Mg s<sup>-1</sup>), and volume flux (m<sup>3</sup> s<sup>-1</sup>) respectively
* ```'La_Yb'``` and ```'Dy_Yb'```: La/Yb and Dy/Yb rare-earth element concentration ratios in generated basalts.

```data``` keys are the names of the parameters listed above, and ```data``` values are presented as a list of two values; the first being the parameter value; the second being its 1 standard deviation uncertainty.

Below is an example ```data``` dictionary for mid-ocean ridge basalt from Matthews et al. (2021):


```python
data = {
       'tc': [5.74, 0.27],
       'Tcrys': [1280, 20],
       'Fpx': [0.175, 0.1]
       }
```

```knowns``` is a dictionary of the parameters that the inversion should not find. These are parameters required by the melting model that are set to be fixed during the inversion. ```knowns``` keys are the parameter names, ```knowns``` values are the parameter values.

```unknowns``` is a dictionary of the parameters to be fond by the inversion. ```unknowns``` keys are the parameter names. ```unknowns``` values are lists comprising of two items: a prior type, and a prior-dependent nested list of either defined lower and upper bounds or a mean and standard deviation (see below):

The parameters that must be defined collectively between ```knowns``` and ```unknowns``` are:

* ```'Tp'```: mantle potential temperature (&deg;C)
* ```'P_lith'```: pressure at the base of the lithosphere (GPa)
* ```'P_cryst'```: pressure at which crystallisation takes place (GPa)
* ```'F_px'```: proportion of pyroxenite in the mantle
* ```'F_hz'```: proportion of harzburgite in the mantle
* ```'DeltaS'```: entropy change on melting (J K<sup>-1</sup>)

If trace element concentrations of a multi-lithology mantle (lherzolite, pyroxenite, harzburgite) are considered then the following parameters must also be defined:

* ```'La_lz'```: lherzolite source lanthanum concentration (ppm)
* ```'Dy_lz'```: lherzolite source dysprosium concentration (ppm)
* ```'Yb_lz'```: lherzolite source ytterbium concentration (ppm)
* ```'La_px'```: pyroxenite source lanthanum concentration (ppm)
* ```'Dy_px'```: pyroxenite source dysprosium concentration (ppm)
* ```'Yb_px'```: pyroxenite source ytterbium concentration (ppm)

If pyroxenite trace compositions are to be approximated by batch melts of DMM composition mantle then ```MORBmelts```, the batch melt fraction of DMM, must be defined.

If the buoyancy of a multi-lithology mantle (lherzolite, pyroxenite, harzburgite) is considered then the following parameters must also be defined:

* ```'ambientTp'```: temperature of ambient mid-ocean ridge mantle (&deg;C)
* ```'ambientPx'```: proportion of pyroxenite in ambient mantle
* ```'ambientHz'```: proportion of harzburgite in ambient mantle

If mantle plume fluxes (either ```'Qm'```, ```'Qb'```, or ```'Qv'```) are considered then the following parameters must also be defined:

* ```'r'```: radius of the mantle plume (m)
* ```'mu'```: plume viscosity (Pa)

If exponentially decaying active upwelling is considered then the following parameters must also be defined:
* ```'lambda'```: wavelength of decay for active upwelling
* ```'amplitude'```: amplitude of decay for active upwelling

The prior types that can be used are:

* ```'uni'```: uniform distribution between a lower and upper bound
* ```'loguni'```: log-uniform distribution between a lower and upper bound
* ```'norm'```: normal distribution defined by a mean and standard deviation
* ```'lognorm'```: log-normal distribution defined by a mean and standard deviation

The following example code is also from the mid-ocean ridge basalt inversion of Matthews et al. (2021):


```python
knowns = {
       'DeltaS': 300,
       'P_lith': 0.0,
       'P_cryst': 0.0,
       }

unknowns = {
       'Tp':['uni', [1250,1600]],
       'F_px':['uni', [0.0,1.0]],
       'F_hz':['uni', [0.0,0.99]]
       }
```

To run MultiNest, an ```inversion``` object must be created:


```python
inv = mi.inversion(
    lithologies=[lz, px, hz],
    data=data,
    knowns=knowns,
    unknowns=unknowns,
    DeltaP=0.004,
    SpreadingCentre=True,
    ContinentalRift=False,
    Passive=True,
    SuperSolidus=True
    Traces=False,
    MORBmelts=False,
    TcrysShallow=True,
    buoyancy=True,
    buoyancyPx='kg1',
    resume=False,
    DensityFile='LithDensity_80kbar.csv',
    livepoints=400,
    name='test_inversion'
    )

inv.run_multinest()
```

This inversion object comprises the following parameters:

* ```lithologies```: list comprising the pyMelt lithology objects in the order lherzolite, pyroxenite, harzburgite.
* ```data```: dictionary of parameters and values for the inversion to match as described above.
* ```knowns```: dictionary of parameters and values that are to be set as fixed as described above.
* ```unknowns```: dictionary of parameters to be found by the inversion as described above.
* ```DeltaP```: a float defining the pressure step of integration during adiabatic decompression melting.
* ```SpreadingCentre```: option to model as a spreading centre, set to ```True``` by default.
* ```ContinentalRift```: option to model as a continental rift, set to ```False``` by default.
* ```Passive```: option to model as passive upwelling, set to ```True``` by default. If set to ```False```, exponentially decaying active upwelling will be considered.
* ```SuperSolidus```: option to include supersolidus melts (```True```), or restrict inversion to subsolidus starts (```False```).
* ```Traces```: option to use La/Yb and Dy/Yb ratios as a constraint, set to ```False``` by default.
* ```MORBmelts```: option to use batch melts of DMM mantle for La/Yb and Dy/Yb compositions of pyroxenite, set to ```False``` by default.
* ```TcrysShallow```: option to use the shallow Tcrys endmember as opposed to the deep endmember, set to ```True``` by default.
* ```buoyancy```: option to necessitate solutions to be buoyant with respect to ambient mantle, set to ```False``` by default.
* ```buoyancyPx```: selects the column in the ```DensityFile``` to use for calculating pyroxenite density in buoyancy calculations, set to ```'kg1'``` of ```'LithDensity_80kbar.csv'``` by default.
* ```resume```: option to resume an incomplete run of MultiNest if files are present, set to ```False``` by default.
* ```DensityFile```: file name for the .csv file containing lithology densities as a function of temperature, set to the included tab-separated repository file ```'LithDensity_80kbar.csv'``` by default.
* ```livepoints```: number of livepoints used by MultiNest, set to ```400``` by default.
* ```name```: the name to call the inversion, and the name of the folder to store results in, set to ```'default'``` by default.

MultiNest may take several hours to run; the run time can be altered by changing the size of the pressure steps undergone during the pyMelt ```adiabaticMelt``` fourth-order Runge-Kutta integration. Output files will be saved to a folder named with the ```name``` string given during the creation of the inversion object. Output files generated include:

* .txt: inversion data presented as a tab-delimited text file which can be imported to Excel (using the 'Get data' function). The first two columns are MultiNest likelihood outputs (sample probability, -2\*loglikelihood); the remaining columns are individual inversion points presented in the order called by pyMelt_MultiNest (this order is also found within the params.json output file).
* stats.dat: summary statistics for the inverted parameters presented in the order called by pyMelt_MultiNest. These statistics are: 
    * Nested sampling global log-evidence for the inversion
    * Nested importance sampling global log-evidence for the inversion
    * Mean and standard deviation of inversion parameters
    * Maximum likelihood of inversion parameters
    * Maximum a posteriori (MAP) estimates of inversion parameters
* stats.json: summary statistics of inversion parameters presented in the order called by pyMelt_MultiNest.
* post_equal_weights.dat: all weighted accepted solutions to the inversion problem.

Upon completion the output files can be interpreted using the multinest_marginals.py script that is packaged with pyMultiNest. The shell command for running multinest_marginals.py is:

```shell
$ python multinest_marginals.py OUTPUT_FOLDER_NAME/
```

Results will be returned as a median value and a standard deviation, with more detailed data provided in the stats.json file generated. Marginal plots are also generated in .pdf and .png format.

## Citing pyMelt and pyMelt_MultiNest

Please cite the pyMelt manuscript in Volcanica. The latest release is v2.0.

Matthews, S., Wong, K. and Gleeson, M. (2022) “pyMelt: An extensible Python engine for mantle melting calculations”, Volcanica, 5(2), pp. 469–475. doi: 10.30909/vol.05.02.469475.

Examples of pyMelt_MultiNest in use (with pyMelt v.1.960 and v1.0 respectively) can be found in the following manuscripts:

Wong, K., Ferguson, D., Matthews, S., Morgan, D., Tadesse, A. Z., Sinetebeb, Y., & Yirgu, G. (2022). Exploring rift geodynamics in Ethiopia through olivine-spinel Al-exchange thermometry and rare-earth element distributions. Earth and Planetary Science Letters, 597, 117820. https://doi.org/10.1016/j.epsl.2022.117820

Matthews, S., Wong, K., Shorttle, O., Edmonds, M., & Maclennan, J. (2021). Do olivine crystallization temperatures faithfully record mantle temperature variability?. Geochemistry, Geophysics, Geosystems, 22(4), e2020GC009157. https://doi.org/10.1029/2020GC009157

If pyMelt is used with pyMultinest, please additionally cite the following papers which detail the development and introduction of MultiNest and pyMultinest:

Feroz, F., & Hobson, M. P. (2008). Multimodal nested sampling: An efficient and robust alternative to Markov Chain Monte Carlo methods for astronomical data analyses. Monthly Notices of the Royal Astronomical Society, 384(2), 449–463. 
https://doi.org/10.1111/j.1365-2966.2007.12353.x

Feroz, F., Hobson, M. P., & Bridges, M. (2009). MultiNest: An efficient and robust Bayesian inference tool for cosmology and particle physics. Monthly Notices of the Royal Astronomical Society, 398(4), 1601–1614. 
https://doi.org/10.1111/j.1365-2966.2009.14548.x

Feroz, F., Hobson, M. P., Cameron, E., & Pettitt, A. N. (2019). Importance Nested Sampling and the MultiNest Algorithm. The Open Journal of Astrophysics, 2(1), 10.21105/astro.1306.2144. 
https://doi.org/10.21105/astro.1306.2144

Buchner, J., Georgakakis, A., Nandra, K., Hsu, L., Rangel, C., Brightman, M., Merloni, A., Salvato, M., Donley, J., & Kocevski, D. (2014). X-ray spectral modelling of the AGN obscuring region in the CDFS: Bayesian model selection and catalogue. Astronomy & Astrophysics, 564, A125.https://doi.org/10.1051/0004-6361/201322971
