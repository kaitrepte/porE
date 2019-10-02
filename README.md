# porE
PORosity Evaluation tool  

Main developer: 

* Kai Trepte 

Sidekick:  

* Sebastian Schwalbe 

Coding language: FORTRAN   

To compile the code, go to the *src* directory and type

make

To run the code, go to the *run* folder. The file run.sh includes three options

* porE: This is the standard code. Includes OSA, and the standard GPA
* porE_subgrid: A modified version of the GPA (GPA_sub-grid). Includes OSA as well.
* get_PSD: Analyze the PSD in any MOF. One needs to be careful with choosing proper numerical parameters

Just uncomment the corresponding line in run.sh and comment out the others.

Usage: First, modify the input file (input_porosity OR input_PSD) according to the system you want to investigate. Then, type
```should work with all shells
./run.sh
```

Limitations: Some elements regarding their vdW radii are not implemented yet. 

More will be done in the future.
