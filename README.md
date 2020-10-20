# porE
![GitHub Logo](/images/porE_logo_v1.png)
## PORosity Evaluation tool

[![license](https://img.shields.io/badge/license-APACHE2-green)](https://www.apache.org/licenses/LICENSE-2.0)
[![language](https://img.shields.io/badge/language-Python3-blue)](https://www.python.org/)
[![language](https://img.shields.io/badge/language-FORTRAN-red)](https://www.fortran.com/)
[![version](https://img.shields.io/badge/version-1.0.3-lightgrey)](https://github.com/kaitrepte/porE/blob/master/README.md)  
[![researchgate](https://img.shields.io/static/v1?label=researchgate&message=MOFs&style=social&logo=researchgate)](https://www.researchgate.net/project/Systematic-and-efficient-theoretical-investigations-of-metal-organic-frameworks-MOFs)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4075260.svg)](https://doi.org/10.5281/zenodo.4075260)


Main developer: 

* Kai Trepte (KT, kai.trepte1987@gmail.com)  
* Sebastian Schwalbe (SS, theonov13@gmail.com)  

## Installation (global pip installation)
```bash 
$ pip3 install porE
```
## Installation (local pip installation)
```bash 
$ git clone https://github.com/kaitrepte/porE
$ cd porE
$ pip3 install -e .
```

### Running the code
There are several calculations types one can choose from:

* HEA : He approach. Calculates the He void fraction using a cell list approach.
* OSA : Overlapping sphere approach. Calculates the porosity via two-body overlaps of spheres.
* GPA : Grid point approach. Void and accessible porosity are calculated using a grid in the unit cell.
* PSD : Pore size distribution. Using a Monte-Carlo scheme, the pore size distribution is calculated.

To run the code, go to the *examples* folder.     
The file run_porE.py contains examples for the execution on the non-GUI level.    
The file run_GUI.py will execute the graphical user interface, including the OSA, GPA and PSD functionalities.     
The HEA will be implemented soon.

Special thanks goes to Sebastian Schwalbe for writing the GUI and for many other useful discussions.

## Additional information
The basic outline of the porE code and its applications are summarized in 

- [K. Trepte and S. Schwalbe, Systematic Analysis of Porosities in Metal-Organic Frameworks](https://chemrxiv.org/articles/Systematic_Analysis_of_Porosities_in_Metal-Organic_Frameworks/10060331)
- ChemRxiv, DOI: 10.26434/chemrxiv.10060331.v1, 2019
- ChemRxiv, DOI: 10.26434/chemrxiv.10060331.v2, 2020
- ChemRxiv, DOI: 10.26434/chemrxiv.10060331.v3, 2020

Also, we create some YouTube videos to introduce and explain the porE code

- [Introduction](https://www.youtube.com/watch?v=yp4IgFnDf9E)

## Reference results
Several results from porE are summarized in the *porE/results* folder.
