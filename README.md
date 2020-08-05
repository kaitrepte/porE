# porE 
## PORosity Evaluation tool

[![license](https://img.shields.io/badge/license-APACHE2-green)](https://www.apache.org/licenses/LICENSE-2.0)
[![language](https://img.shields.io/badge/language-Python3-blue)](https://www.python.org/)
[![language](https://img.shields.io/badge/language-FORTRAN-blue)](https://www.fortran.com/)
[![version](https://img.shields.io/badge/version-1.0.1-lightgrey)](https://github.com/kaitrepte/porE/blob/master/README.md)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3972196.svg)](https://doi.org/10.5281/zenodo.3972196)

Main developer: 

* Kai Trepte 
* Sebastian Schwalbe 

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
There are three calculations types one can choose from:

* OSA : Overlapping sphere approach. Calculates the porosity via two-body overlaps of spheres.
* GPA : Grid point approach. Void and accessible porosity are calculated using a grid in the unit cell.
* PSD : Pore size distribution. Using a Monte-Carlo scheme, the pore size distribution is calculated.

To run the code, go to the *examples* folder.     
The file run_porE.py contains examples for the execution on the non-GUI level.    
The file run_GUI.py will execute the graphical user interface, including all functionalities.   

Special thanks goes to Sebastian Schwalbe for writing the GUI and for many other useful discussions.    

## Additional information
Limitations: Some elements regarding their vdW radii are not implemented yet.    
More will be done in the future.

The basic outline of the porE code and its applications are summarized in 

-  [K. Trepte and S. Schwalbe, Systematic Analysis of Porosities in Metal-Organic Frameworks](https://chemrxiv.org/articles/Systematic_Analysis_of_Porosities_in_Metal-Organic_Frameworks/10060331), ChemRxiv, DOI: 10.26434/chemrxiv.10060331.v1, 2019 



## Reference results
Several results from porE are summarized in the *porE/results* folder.
