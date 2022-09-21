## pyMelt_MultiNest
The latest update to the repository should make it compatible with pyMelt v2.0, but it has not been seriously tested. There are also many possibilities for adding new constraints and streamlining the inversions which have not been taken advantage of here. It is probably best to proceed by making a bespoke version of the module (starting with this template) to solve your particular problem. Feel free to reach out to me (simonm@hi.is) for help with doing this.

pyMelt can be used in conjunction with the MultiNest algorithm (Feroz and Hobson, 2008; Feroz et al., 2009, 2013) via its python frontend, pyMultinest (Buchner et al., 2014). This permits the inversion of measured data (e.g. crystallisation temperature, crustal thickness) to obtain unknowns (e.g. potential temperature) via Bayesian inference. More details of the inversion methods are provided in Matthews et al. (in review).

For pyMelt_MultiNest to work, MultiNest and pyMultinest must be installed. The user is directed to the [pyMultinest installation instructions](https://johannesbuchner.github.io/PyMultiNest/) for further guidance.

pyMelt_MultiNest is installed by placing the pyMelt_multinest.py file in the working directory: 

```python
import pyMelt_multinest as mi
```

If operating on a Linux installation, the latest version can be downloaded from the [code repository](https://github.com/simonwmatthews?tab=repositories) and installed locally:

```shell
$ git clone https://github.com/simonwmatthews/pyMelt_multinest
```

First define the lithology objects to be used in the inversion:


```python
lz = m.lithologies.matthews.klb1()
px = m.lithologies.matthews.kg1()
hz = m.lithologies.shorttle.harzburgite()
```

Next the data, knowns, and unknowns must be specified.

```data``` is a dictionary of the initial parameters that the model should aim to match. This data includes:

* ```'Tcrys'```: melt liquidus/crystallisation temperature (&deg;C)
* ```'tc'```: crustal thickness (km)
* ```'Fpx'```: fraction of aggregate melts derived from pyroxenite
* ```'Q'```: melt flux (m<sup>3</sup> s<sup>-1</sup>)

```data``` keys are the names of the parameters listed above, and ```data``` values are presented as a list of two values; the first being the parameter value; the second being its 1 standard deviation uncertainty.

Below is an example ```data``` dictionary for mid-ocean ridge basalt from Matthews et al. (in review):


```python
data = {
       'tc': [5.74, 0.27],
       'Tcrys': [1280, 20],
       'Fpx': [0.175, 0.1]
       }
```

```knowns``` is a dictionary of the parameters that the inversion should not find. These are parameters required by the melting model that are set to be fixed during the inversion. ```knowns``` keys are the parameter names, ```knowns``` values are the parameter values.

```unknowns``` is a dictionary of the parameters to be fond by the inversion. ```unknowns``` keys are the parameter names. ```unknowns``` values are lists comprising of two items: a prior type, and a prior-dependent nested list of either defined lower and upper bounds or a mean and standard deviation (see below):

```python
unknowns = {
        'Tp': ['uni',[lower_bound, upper_bound]],
        'P_lith': ['norm',[value_mean, value_stdev]]
        }
```

The parameters that must be defined collectively between the knowns and the unknowns are:

* ```'Tp'```: mantle potential temperature (&deg;C)
* ```'P_lith'```: pressure at the base of the lithosphere (GPa)
* ```'P_cryst'```: pressure at which crystallisation takes place (GPa)
* ```'F_px'```: proportion of pyroxenite in the mantle
* ```'F_hz'```: proportion of harzburgite in the mantle
* ```'DeltaS'```: entropy change on melting (J K<sup>-1</sup>)

Additionally, the following conditions must be met.

* If mantle buoyancy is being considered, the ambient mantle potential temperature and the fractions of pyroxenite and harzburgite in the mantle must be included as knowns or unknowns. 
* If the data to be matched includes mantle plume melt flux (```'Q'``` in m<sup>3</sup> s<sup>-1</sup>), mantle plume radius (```'r'``` in m) and the plume viscosity (```'mu'``` in Pa s<sup>-1</sup>) must be included as knowns or unknowns.

The prior types that can be used are:

* ```'uni'```: uniform distribution between a lower and upper bound
* ```'loguni'```: log-uniform distribution between a lower and upper bound
* ```'norm'```: normal distribution defined by a mean and standard deviation
* ```'lognorm'```: log-normal distribution defined by a mean and standard deviation

The following example code is also from the MORB values of Matthews et al. (in review). In this case all priors are uniform:


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

To run MultiNest, an ```inversion``` object must be created.

```python
inversion(self, lithologies, data, knowns, unknowns, SpreadingCenter=True, TCrysShallow=True, buoyancy=False, buoyancyPx='kg1', resume=False, DensityFile'LithDensity_80kbar.csv', livepoints=400, name='default')
```

This inversion object comprises the following parameters:

* ```lithologies```: list comprising the pyMelt lithology objects in the order lherzolite, pyroxenite, harzburgite.
* ```data```: dictionary of parameters and values for the inversion to match as described above.
* ```knowns```: dictionary of parameters and values that are to be set as fixed as described above.
* ```unknowns```: dictionary of parameters to be found by the inversion as described above.
* ```SpreadingCenter```: option to model as a spreading center, set to ```True``` by default.
* ```TcrysShallow```: option to use the shallow Tcrys endmember as opposed to the deep endmember, set to ```True``` by default.
* ```buoyancy```: option to necessitate solutions to be buoyant with respect to ambient mantle, set to ```False``` by default.
* ```buoyancyPx```: selects the column in the ```DensityFile``` to use for calculating pyroxenite density in buoyancy calculations, set to ```'kg1'``` of ```'LithDensity_80kbar.csv'``` by default.
* ```resume```: option to resume an incomplete run of MultiNest if files are present, set to ```False``` by default.
* ```DensityFile```: file name for the .csv file containing lithology densities as a function of temperature, set to the included tab-separated repository file ```'LithDensity_80kbar.csv'``` by default.
* ```livepoints```: number of livepoints used by MultiNest, set to ```400``` by default.
* ```name```: the name to call the inversion, and the name of the folder to store results in, set to ```'default'``` by default.

An example ```inversion``` object is created as follows:


```python
inv = mi.inversion([lz,px,hz],data,knowns,unknowns,name='MORB')
```

MultiNest can then be run upon creation of the ```inversion``` object:


```python
inv.run_multinest()
```

A commented example pyMelt_MultiNest code snippet is provided below for the new Hawaii data published in Matthews et al. (in review):


```python
# Values in km for lithospheric thickness and crystallisation depth
tlith = 75
tlith_sd = 5
tc = 18
tc_sd = 1

# Densities in g/cm3 for lithosphere and crust
rho_lith = 2.8
rho_c = 2.6

# Gravity in m/s2
g_c = 10.0

# Observed data:
data = {
        'Tcrys': [1460, 20],
        'Q': [16, 2]
       }

# Fixed known values
knowns = {
        'DeltaS': 300,
        'mu': 1e19
         }

# Unknowns to be found by the inversion
# Ambient values are taken from the MORB inversion
unknowns = {
            'Tp': ['uni',[1250, 2000]],
            'F_px': ['uni',[0, 1]],
            'F_hz': ['uni',[0, 1]],
            'P_lith': ['norm',[tlith*rho_lith*g_c/1000, tlith_sd*rho_lith*g_c/1000]],
            'P_cryst': ['norm',[tc*rho_c*g_c/1000, tc_sd*rho_c*g_c/1000]],
            'r': ['loguni', [50*1e3, 300*1e3]],
            'ambientTp': ['norm',[1364, 15/2]],
            'ambientPx': ['norm',[0.021, 0.012]],
            'ambientHz': ['norm',[0.42, 0.15]]
           }

# Create inversion object
inv = mi.inversion([lz,px,hz],data,knowns,unknowns,name='Hawaii',
                               buoyancy=True,buoyancyPx='kg1',SpreadingCenter=False)

# Run MultiNest
inv.run_multinest()
```

MultiNest may take several hours to run; the run time can be altered by changing the number of pressure steps undergone during the ```AdiabaticMelt_1D``` fourth-order Runge-Kutta integration. Output files will be saved to a folder named with the string given during the creation of the inversion object. Output files generated include:

* .txt: inversion data presented as a tab-delimited text file which can be imported to Excel (using the 'Get data' function). The first two columns are MultiNest likelihood outputs (sample probability, -2\*loglikelihood); the remaining columns are individual inversion points presented in the order called by pyMelt_MultiNest (the order ```'Tp','DeltaS','P_lith','P_cryst','F_px','F_hz', 'ambientTp','ambientPx','ambientHz','r','mu'``` if all known/unknown variables are called; this order is also found within the params.json output file).
* stats.dat: summary statistics for the inverted parameters presented in the order called by pyMelt_MultiNest. These statistics are: 
    * Nested sampling global log-evidence for the inversion
    * Nested importance sampling global log-evidence for the inversion
    * Mean and standard deviation of inversion parameters
    * Maximum likelihood of inversion parameters
    * Maximum a posteriori (MAP) estimates of inversion parameters
* stats.json: summary statistics of inversion parameters presented in the order called by pyMelt_MultiNest.

Upon completion the output files can be interpreted using the multinest_marginals.py script that is packaged with pyMultiNest. The shell command for running multinest_marginals.py is:

```shell
$ python multinest_marginals.py OUTPUT_FOLDER_NAME/
```

Results will be returned as a median value and a standard deviation, with more detailed data provided in the stats.json file generated. Marginal plots are also generated in .pdf and .png format.

## Citing pyMelt

If pyMelt is used with pyMultinest, please additionally cite the following papers which detail the development and introduction of MultiNest and pyMultinest:

Feroz, F., & Hobson, M. P. (2008). Multimodal nested sampling: An efficient and robust alternative to Markov Chain Monte Carlo methods for astronomical data analyses. Monthly Notices of the Royal Astronomical Society, 384(2), 449–463. 
https://doi.org/10.1111/j.1365-2966.2007.12353.x

Feroz, F., Hobson, M. P., & Bridges, M. (2009). MultiNest: An efficient and robust Bayesian inference tool for cosmology and particle physics. Monthly Notices of the Royal Astronomical Society, 398(4), 1601–1614. 
https://doi.org/10.1111/j.1365-2966.2009.14548.x

Feroz, F., Hobson, M. P., Cameron, E., & Pettitt, A. N. (2019). Importance Nested Sampling and the MultiNest Algorithm. The Open Journal of Astrophysics, 2(1), 10.21105/astro.1306.2144. 
https://doi.org/10.21105/astro.1306.2144

Buchner, J., Georgakakis, A., Nandra, K., Hsu, L., Rangel, C., Brightman, M., Merloni, A., Salvato, M., Donley, J., & Kocevski, D. (2014). X-ray spectral modelling of the AGN obscuring region in the CDFS: Bayesian model selection and catalogue. Astronomy & Astrophysics, 564, A125.https://doi.org/10.1051/0004-6361/201322971
