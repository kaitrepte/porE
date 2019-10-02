#!/bin/sh

# Author:  K. Trepte
# Date: 02.10.2019 
porE_dir="$(dirname "$(pwd)")"

## Execute porE 
## Functionalities: OSA, normal GPA
#$porE_dir/src/porE

## Execute porE_subgrid
## Functionalities: OSA, GPA with sub-grid division (much faster)
$porE_dir/src/porE_subgrid

## Execute get_PSD
## Functionalities: Calculate the PSD with any desired accuracy
#$porE_dir/src/get_PSD
