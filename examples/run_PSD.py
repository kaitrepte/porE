from porE.hea.HEA import HEA
from pore import porosity as p
from pore import psd 
from porE.io.ase2pore import *
# start from a cif file 
cif = 'structures/cif/porE8/mof210_vesta.cif'
ase2pore(cif)
structure = 'pypore.xyz'
#
# Execute an analyis of the pore size distribution (PSD)
# First number   -> Number of starting points
# Second number  -> Number of MC steps
#
print('-----------')
print('\nRun PSD\n')
print('-----------')
no_pores,tmp1,tmp2,tmp3,tmp4 = psd.get_psd(structure,200,10000)
pores           = tmp1[0:no_pores]
distr           = tmp2[0:no_pores]
pore_pos_cart   = tmp3[0:no_pores]
pore_pos_frac   = tmp4[0:no_pores]
