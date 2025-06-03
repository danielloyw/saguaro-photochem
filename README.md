# Saguaro Atmospheric Photochemical Model
_Saguaro_ is a Fortran-based program for modeling photochemistry in planetary atmospheres. Currently, the implementation for Mars is available in this repository.

**Features:**
* 1-dimensional
* Coupled ion-neutral photochemistry with >60 neutral and ionic species, linked by >800 reactions including unimolecular, association, photodissociation, photoionization, electron impact, and dissociative recombination reactions
* Calculation of suprathermal electron population, including cascade through secondary collisions and ionization events
* Calculation of CO photodissociation with high spectral resolution to account for pre-dissociating states
* Transport by molecular and eddy diffusion
* Spans the full atmosphere from surface to exobase, or a more limited altitude range of interest
* Flexible definition of initial atmosphere, solar ultraviolet flux, escape, and surface exchange
* Calculation of both time-dependent and steady-state reaction rates and number densities

The _Saguaro_ codebase and this repository is maintained by Daniel Lo (https://github.com/danielloyw). Please direct any queries to danielloyw@gmail.com.


## Brief History
_Saguaro_ started out as a series of neutral and ion models focused on minor species in the ionosphere of Titan<sup>[1-4]</sup>. These initial models were then combined into a self-consistent ion-neutral chemistry model<sup>[5]</sup>, with a significant expansion in the N chemistry<sup>[6]</sup> and incorporation of molecular and eddy diffusion<sup>[7]</sup>. Subsequent developments include a dedicated module for detailed calculations of suprathermal electrons, and the use of high-resolution cross sections to account for pre-dissociating states in N<sub>2</sub> photodissociation<sup>[8]</sup>. This Titan implementation has been used to study a large variety of hydrocarbon and other organic species in the Titan atmosphere<sup>[9-11]</sup>. 

_Saguaro-Mars_ was adapted from the original Titan implementation. With a revised set of reactions, atmospheric parameters and initial atmospheric composition, _Saguaro-Mars_ has been used to study C photochemistry and escape<sup>[12-14]</sup>, H photochemistry and escape<sup>[15]</sup>, and O<sub>2</sub> variability<sup>[16]</sup> at Mars. 

A list of journal [references](#references) is available at the end of this README. 

# Model Inputs

# Model Outputs
Output data files produced by _Saguaro_ are located in the `/runs/<run-name>/output/` directory, with `<run-name>` being the name of the run set by the user when commencing the model run. These files use ASCII encoding, and can be opened by any text editor. Typically, each file begins with a header containing numbers describing the data structure. Each variable is described by a block with the variable name as header, and data entries (usually the variable as a function of altitude) within each block are split into lines of 10. 

## Description of Output Files
### General
| Filename         | Description  |
| :---------       | :----------- |
| `atm1D.out`      | Altitude profiles for key variables in the model, including the altitude grid, gravity, neutral and electron temperatures, pressure, mean molecular weight, eddy diffusion coefficient, and densities for each species. Same structure as input `atm1D.in` file. |
| `atm1D.csv`      | Same as `atm1D.out`, except in comma-separated format, and mixing ratios for each species are provided instead of their densities. |
| `solar_flux.out` | Solar flux as a function of wavelength and altitude. Units: photons cm<sup>-2</sup> s<sup>-1</sup> nm<sup>-1</sup>.                |
| `eflux.out`      | Electron fluxes as a function of energy and altitude. Units: electrons cm<sup>-2</sup> s<sup>-1</sup> eV<sup>-1</sup>.             |
| `balance.out`    | Overview of convergence for model species, with information about column number density (cm<sup>-2</sup>), production and loss rates (cm<sup>-2</sup> s<sup>-1</sup>), and fluxes at top and bottom (cm<sup>-2</sup> s<sup>-1</sup>). Balance (cm<sup>-2</sup> s<sup>-1</sup>) across the various contributing terms is calculated, and the stability time constant (s) is computed from dividing the balance with the column number density. |
| `mcolrates.out`  | Same as `balance.out` without some columns.                                                                                        |
| `diff.csv`       | Diffusion coefficients for each species with altitude. The eddy diffusion coeffient (which applies to all species) is also provided. Units: cm<sup>2</sup> s<sup>-1</sup>. |


### Reaction Rates
| Filename         | Description  |
| :---------       | :----------- |
| `ratecoeff.out`  | Reaction coefficients for each chemical reaction with altitude. Units: s<sup>-1</sup> for unimolecular reactions and cm<sup>3</sup> s<sup>-1</sup> for bimolecular reactions. |
| `chemrates.out`  | Reaction rates for each chemical reaction with altitude. Units: cm<sup>-3</sup> s<sup>-1</sup>. |
| `colrates.out`   | Table of column-integrated rates for each chemical reaction, sorted from highest to lowest. The columns of the table correspond to: (1) reaction index, (2) column-integrated rate down to , (3) column-integrate rate below , (4) total column-integrated rate, and (5) reaction. |
| `photorates.out` | Reaction rates for each photo reaction with altitude. Units: cm<sup>-3</sup> s<sup>-1</sup>. |
| `pcolrates.out`  |
| `elerates.out`   | Reaction rates for each electron reaction with altitude. Units: cm<sup>-3</sup> s<sup>-1</sup>. |
| `ecolrates.out`  |

### Species-Specific Files
Additonal information about each calculated species in the model is available in `/molecules/<species-name>.out`, where `<species-name>` represents the described species. 

| Column         | Units     | Description |
| :---------     | :-------- | :---------- | 
| alt     
| den        
| mole       
| flux       
| ext prd    
| pr_ph      
| ls_ph      
| pr_pe      
| ls_pe     
| pr_chem   
| ls_chem    
| net prod   
| net loss    
| condense   
| -div_flx   
| balance


## Reading Output Files
A variety of functions are available in `saguaro-photochem/mars/plot/saguaro_read.py` for Python users to import into their own code for reading the output files. The functions are written in Python 3, and return a dictionary of the data in the files stored as `numpy` arrays and `pandas` DataFrames. More information about the specific functions is available below.

### Python Dependencies
* NumPy
* pandas

### `saguaro_read.atm(filename)`
Reads in variables describing general atmospheric structure and number densities of all species from `atm1D.in` and `atm1D.out`.

#### Parameters
| Name           | Type    | Required | Default | Description                           |
| :---------     | :----   | :------- | :------ | :------------------------------------ |
| `filename`     | `str`   | Yes      | N/A     | Path of `atm1D.in` / `atm1D.out`.     |

#### Returns (Value Structure)
| Key            | Type               | Description                                     |
| :------------- | :----------------- | :---------------------------------------------- |
| `n_alt`        | `int`              | Number of altitude bins.                        |
| `n_mol`        | `int`              | Number of species.                              |
| `Species`      | `str list`         | List of species.                                |
| `Profiles`     | `pandas DataFrame` | Altitude profiles of atmospheric variables and number densities of all species. Columns are variables while rows are altitudes. |


### `saguaro_read.molecules(filename)`
Reads in parameters for modeling neutral (`/input/nmolecules.dat`) and ionic (`/input/imolecules.dat`) species.

#### Parameters
| Name           | Type    | Required | Default | Description                                      |
| :---------     | :----   | :------- | :------ | :----------------------------------------------- |
| `filename`     | `str`   | Yes      | N/A     | Path of `nmolecules.dat` / `imolecules.dat`.     |

#### Returns (Value Structure)
| Key            | Type               | Description                                                |
| :------------- | :----------------- | :--------------------------------------------------------- |
| `Species`      | `pandas DataFrame` | Table of species and how they are calculated in the model. |


### `saguaro_read.rates(filename)`
Reads in reaction rates for chemical (`chemrates.out`), photo (`photorates.out`), and electron (`elerates.out`) reactions. Also reads in rate coefficients for the chemical reactions (`ratecoeff.out`). 

#### Parameters
| Name           | Type    | Required | Default | Description                           |
| :---------     | :----   | :------- | :------ | :------------------------------------ |
| `filename`     | `str`   | Yes      | N/A     | Path of `chemrates.out` / `photorates.out` / `elerates.out` / `ratecoeff.out`. |

#### Returns (Value Structure)
| Key            | Type               | Description                                     |
| :------------- | :----------------- | :---------------------------------------------- |
| `Altitude`     | `float` array      | List of altitudes.                              |
| `Reaction`     | `pandas DataFrame` | List of reactions.                              |
| `Rates`        | `float` array      | Array of reaction rates / coefficients, with dimensions (altitude, reaction).            |


### `saguaro_read.eflux(filename)`
Reads in electron fluxes from `eflux.out`.

#### Parameters
| Name           | Type    | Required | Default | Description                           |
| :---------     | :----   | :------- | :------ | :------------------------------------ |
| `filename`     | `str`   | Yes      | N/A     | Path of `eflux.out`.                  |

#### Returns (Value Structure)
| Key            | Type               | Description                                     |
| :------------- | :----------------- | :---------------------------------------------- |
| `Altitude`     | `float` array      | List of altitudes.                              |
| `Energy`       | `float` array      | List of energies.                               |
| `Flux`         | `float` array      | Array of electron fluxes, with dimensions (altitude, energy). |


### `saguaro_read.summary(filename)`
Reads in variables associated with specific species (`*.out` files in `/output/molecules/`).

#### Parameters
| Name           | Type    | Required | Default | Description                           |
| :---------     | :----   | :------- | :------ | :------------------------------------ |
| `filename`     | `str`   | Yes      | N/A     | Path of input file.                   |

#### Returns (Value Structure)
| Key            | Type               | Description                                     |
| :------------- | :----------------- | :---------------------------------------------- |
| `Profiles`     | `pandas DataFrame` | Altitude profiles of variables associated with specified species. Columns are variables while rows are altitudes. |


## Plotting Model Outputs
A variety of functions are available in `saguaro-photochem/mars/plot/saguaro_read.py` for Python users to import into their own code for reading the output files. The functions are written in Python 3, and return a dictionary of the data in the files stored as `numpy` arrays and `pandas` DataFrames. More information about the specific functions is available below.

### Python Dependencies
* NumPy
* Matplotlib
* pandas
* pypdf
* sys
* os
* warnings


## References

<sup>[1]</sup> Vuitton V., et al. (2006). _The nitrogen chemistry of Titan’s upper atmosphere revealed_. The Astrophysical Journal. https://doi.org/10.1086/507467

<sup>[2]</sup> Vuitton V., et al. (2007). _Ion chemistry and N-containing molecules in Titan’s upper atmosphere_. Icarus. https://doi.org/10.1016/j.icarus.2007.06.023

<sup>[3]</sup> Vuitton V., et al. (2008). _Formation and distribution of benzene on Titan_. Journal of Geophysical Research: Planets. https://doi.org/10.1029/2007JE002997

<sup>[4]</sup> Vuitton V., et al. (2009). _Negative ion chemistry in Titan’s upper atmosphere_. Planetary and Space Science. https://doi.org/10.1016/j.pss.2009.04.004

<sup>[5]</sup> Yelle R., et al. (2010). _Formation of NH<sub>3</sub> and CH<sub>2</sub>NH in Titan’s upper atmosphere_. Faraday Discussions. https://doi.org/10.1039/C004787M

<sup>[6]</sup> Lavvas P., et al. (2008). _Coupling photochemistry with haze formation in Titan’s atmosphere. Part I: Model description_. Planetary and Space Science. https://doi.org/10.1016/j.pss.2007.05.027

<sup>[7]</sup> Hörst S., et al. (2008). _The origin of oxygen species in Titan’s atmosphere_. Journal of Geophysical Research: Planets. https://doi.org/10.1029/2008JE003135

<sup>[8]</sup> Lavvas P., et al. (2011). _Energy deposition and primary chemical products in Titan’s upper atmosphere_. Icarus. https://doi.org/10.1016/j.icarus.2011.03.001

<sup>[9]</sup> Vuitton V., et al. (2012). _Rapid association reactions at low pressure: impact on the formation of hydrocarbons on Titan_. The Astrophysical Journal. http://dx.doi.org/10.1088/0004-637X/744/1/11

<sup>[10]</sup> Lavvas P., et al. (2015). _N<sub>2</sub> state population in Titan’s atmosphere_. Icarus. https://doi.org/10.1016/j.icarus.2015.06.033

<sup>[11]</sup> Vuitton V., et al. (2018). _Simulating the density of organic species in the atmosphere of Titan with a coupled ion-neutral photochemical model_. Icarus. https://doi.org/10.1016/j.icarus.2018.06.013

<sup>[12]</sup> Lo D., et al. (2020). _Carbon photochemistry at Mars: Updates with recent data_. Icarus. https://doi.org/10.1016/j.icarus.2020.114001.

<sup>[13]</sup> Lo D., et al. (2021). _Carbon photochemical escape rates from the modern Mars atmosphere_.” Icarus. https://doi.org/10.1016/j.icarus.2021.114371.

<sup>[14]</sup> Lo D., et al. (2022). _MAVEN/IUVS observations of C I 156.1 nm and 165.7 nm dayglow: Direct detection of carbon and implications on photochemical escape_. Icarus. https://doi.org/10.1016/j.icarus.2021.114664.

<sup>[15]</sup> Stone S., et al. (2020). _Hydrogen escape from Mars is driven by seasonal and dust storm transport of water_. Science. https://doi.org/10.1126/science.aba5229.

<sup>[16]</sup> Lo D., et al. (2024). _Evaluating atmospheric and surface drivers for O<sub>2</sub> variations at Gale crater as observed by MSL SAM_. The Planetary Science Journal. https://doi.org/10.3847/PSJ/ad251b.
