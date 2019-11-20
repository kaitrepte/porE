# porE
PORosity Evaluation tool  

Main developer: 

* Kai Trepte 

Sidekick:  

* Sebastian Schwalbe 

Coding language: FORTRAN   

## Installation
To compile the code, do

	cd src/

	make

	cd ..

## Running the code
There are three calculations types one can choose from:

* OSA : Overlapping sphere approach. Calculates the porosity via two-body overlaps of spheres.
* GPA : Grid point approach. Void and accessible porosity are calculated using a grid in the unit cell.
* PSD : Pore size distribution. Using a Monte-Carlo scheme, the pore size distribution is calculated.

To run the code, go to the *run* folder. The file run.sh includes these options:

* porE: This is the standard code. Includes OSA, and the standard GPA
* porE_subgrid: A modified, much faster version of the GPA (GPA_sub-grid). Includes OSA as well.
* get_PSD: Analyze the PSD in any MOF. One needs to be careful with choosing proper numerical parameters

Just uncomment the corresponding line in run.sh and comment out the others.

Usage: First, modify the input file (input_porosity OR input_PSD) according to the system you want to investigate. 

Then, type
```should work with all shells
./run.sh
```

## Additional information
Limitations: Some elements regarding their vdW radii are not implemented yet. 
More will be done in the future.

The basic outline of the porE code and its applications are summarized in 

K. Trepte and S. Schwalbe, 'Systematic Analysis of Porosities in Metal-Organic Frameworks', 2019, at ChemRxiv, see 
https://chemrxiv.org/articles/Systematic_Analysis_of_Porosities_in_Metal-Organic_Frameworks/10060331


## Reference results
Several results from porE are summarized in the *results* folder.
