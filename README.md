# porE
PORosity Evaluation tool  

Main developer: 

* Kai Trepte 

Sidekick:  

* Sebastian Schwalbe 

Coding language: 

* FORTRAN (no GUI)
* FORTRAN, Python (with GUI)

## Installation to run in shell (no GUI)
To compile the code, do

	cd src/

	make

	cd ..

### Running the code
There are three calculations types one can choose from:

* OSA : Overlapping sphere approach. Calculates the porosity via two-body overlaps of spheres.
* GPA : Grid point approach. Void and accessible porosity are calculated using a grid in the unit cell.
* PSD : Pore size distribution. Using a Monte-Carlo scheme, the pore size distribution is calculated.

To run the code, go to the *examples* folder. The file run.sh includes these options:

* porE: This is the standard code. Includes OSA, and the standard GPA
* porE_subgrid: A modified, much faster version of the GPA (GPA_sub-grid). Includes OSA as well. Further, there is an evaluation of the pore windows, which uses the output of get_PSD (see below)
* get_PSD: Analyze the PSD in any MOF. One needs to be careful with choosing proper numerical parameters

Just uncomment the corresponding line in run.sh and comment out the others.

* First, modify the input file (input_porosity OR input_PSD) according to the system you want to investigate. 
* When pore windows shall be computed, one needs an output_PSD file. This file is generated by executing get_PSD. 
Beware that for any given structure, this only needs to be done once. For the given structures, one can copy the corresponding files from results/results_PSD.

Once all is set, type
```should work with all shells
./run.sh
```

## Installation to run with a Python GUI
In the main folder of porE, there is a directory called *pyporE*. It contains a /src directory with all source files, 
and an executable *pypore.py*. This executable allows using a graphical user interface (GUI) for porE. 

All functionalities as mentioned above are included in this GUI.

To compile this code, do

	cd pyporE/src/

	bash compile_python.sh


### Running the code
After compiling the code, you can simply execute

	python3 pypore.py

to start the GUI. Follow the tool-tips to run any desired calculation. This GUI also allows to load .cif files.
In addition, pyporE can be used as a python library. An example is provided in the folder /pypore_as_library.

Special thanks goes to Sebastian Schwalbe for writing this GUI and for many other useful discussions.

## Additional information
Limitations: Some elements regarding their vdW radii are not implemented yet. 
More will be done in the future.

The basic outline of the porE code and its applications are summarized in 

K. Trepte and S. Schwalbe, 'Systematic Analysis of Porosities in Metal-Organic Frameworks', 2019, at ChemRxiv, see 
https://chemrxiv.org/articles/Systematic_Analysis_of_Porosities_in_Metal-Organic_Frameworks/10060331


## Reference results
Several results from porE are summarized in the *results* folder.
