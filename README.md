# porE
PORosity Evaluation tool  

Main developer: 

* Kai Trepte 

Sidekick:  

* Sebastian Schwalbe 

Coding language: FORTRAN   

To compile the code, go to the *src* directory and type

        bash compile.sh

To run the code, go to the *run* folder. The file run.sh includes three options

* porE: This is the standard code. Includes OSA, standard GPA and an evaluation of the PSD (beta version)
* porE_subgrid: A modified version of the GPA (GPA_sub-grid). Includes OSA as well, but not the PSD evaluation
* get_PSD: Analyze the PSD in any MOF. One needs to be careful with choosing proper numerical parameters

Just uncomment the corresponding line in run.sh and comment out the others.


Limitations: So far, only the vdW radii of the following elements are implemented: H, C, N, O, Ni, Cu, Zn, Zr

More will be done in the future.
