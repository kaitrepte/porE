#!/bin/bash

# old porE 
f2py3 -c -m porE porE.f90
mv porE.*.so porE.so

# porE subgrid 
f2py3 -c -m porE_subgrid porE_subgrid.f90
mv porE_subgrid.*.so porE_subgrid.so

# get_PSD
f2py3 -c -m get_PSD get_PSD.f90
mv get_PSD.*.so get_PSD.so

# pore_window 
f2py3 -c -m porE_window porE_window.f90
mv porE_window.*.so porE_window.so

