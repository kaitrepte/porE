# porE
PORosity Evaluation tool  

Main developer: 

* Kai Trepte 

Sidekick:  

* Sebastian Schwalbe 

Coding language: FORTRAN   

To compile the code, go to the *src* directory and type

        bash compile.sh

To run the code, type

	./porE


This tool provides a simple way to study porosities in e.g. metal-organic frameworks (MOFs). 
All that is needed is an xyz-file of the coordinates of the structure and the corresponding cell vectors (see structures_xyz for examples). 
In addition, there is the possibility to evaluate the pore size distribution (beta version).

Limitations: So far, only the vdW radii of the following elements are implemented: H, C, N, O, Ni, Cu, Zn, Zr

More will be done in the future.


# find_pores
An independent tool is provided, called find_pores.f90 . 
It can be used to determine the pore sizes in a structure and get an estimate of the pore size distribution.
This tool is in a developer version and should be used with care.
