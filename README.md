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

This repository is maintained by Daniel Lo (https://github.com/danielloyw). Please direct any queries to danielloyw@gmail.com.

# Brief History
_Saguaro_ started out as a series of neutral and ion models focused on minor species in the ionosphere of Titan<sup>[1-4]</sup>. These initial models were then combined into a self-consistent ion-neutral chemistry model<sup>[5]</sup>, with a significant expansion in the N chemistry<sup>[6]</sup> and incorporation of molecular and eddy diffusion<sup>[7]</sup>. Subsequent developments include a dedicated module for detailed calculations of suprathermal electrons, and the use of high-resolution cross sections to account for pre-dissociating states in N<sub>2</sub> photodissociation<sup>[8]</sup>. This Titan implementation has been used to study a large variety of hydrocarbon and other organic species in the Titan atmosphere<sup>[9-11]</sup>. 

_Saguaro-Mars_ was adapted from the original Titan implementation. With a revised set of reactions, atmospheric parameters and initial atmospheric composition, _Saguaro-Mars_ has been used to study C photochemistry and escape<sup>[12-14]</sup>, H photochemistry and escape<sup>[15]</sup>, and O<sub>2</sub> variability<sup>[16]</sup> at Mars. 

A list of [references](#references) is available at the end of this README. 

# Model Inputs

## Reading Output Files
Output data files produced by _Saguaro_, with extensions `.out` and `.csv`, are located in the `/runs/<run_name>/output/` directory, with `<run_name>` being the name of the run set by the user when commencing the model run. These files use ASCII encoding, and can be opened by any text editor. 

A variety of functions are available in `saguaro-photochem\mars\plot\saguaro_read.py` for Python users to import into their own code for reading the output files. The functions are written in Python 3, and return a dictionary of the data in the files stored as `numpy` arrays and `pandas` DataFrames. More information about the specific functions is available below:

### `saguaro_read.atm(filename)`
Reads in variables describing general atmospheric structure and number densities of all species from `atm1D.in` and `atm1D.out`

#### Parameters

| Name       | Type  | Required | Default | Description                                 |
| :--------- | :---- | :------- | :------ | :------------------------------------------ |
| `filename` | `str` | Yes      | N/A     | Path of `atm1D.in` / `atm1D.out`            |

#### Returns (Value Structure)

| Key            | Type    | Description                               |
| :------------- | :------ | :---------------------------------------- |
| `result`       | `float` | The calculated numerical output.          |
| `status`       | `str`   | Status of the operation (`'success'`). |



# Read files for specific species (*.out files in /molecules/)
def summary(filename):

# Read nmolecules.dat and imolecules.dat
def molecules(filename):

# Read chemrates.out, photorates.out, elerates.out, ratecoeff.out
def rates(filename):

# Read eflux.out
def eflux(filename):




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
