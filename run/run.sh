#!/bin/sh

# Author:  K. Trepte
# Date: 24.01.2020
porE_dir="$(dirname "$(pwd)")"

## Execute porE 
## Functionalities: OSA, normal GPA
#$porE_dir/src/porE

## Execute porE_subgrid
## Functionalities: OSA, GPA with sub-grid division (much faster), pore window evaluation (needs an output_PSD file as input)
## Usage for pore window: First, generate an output_PSD file using the get_PSD exectuable. Then, you can use this executable to evaluate the pore windows
$porE_dir/src/porE_subgrid

## Execute get_PSD
## Functionalities: Calculate the PSD with any desired accuracy
#$porE_dir/src/get_PSD
