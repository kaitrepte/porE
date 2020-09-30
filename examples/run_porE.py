from porE.hea.HEA import HEA
from pore import porosity as p
from pore import psd 
from porE.io.ase2pore import *
# start from a cif file 
cif = 'structures/cif/uio66_vesta.cif'
#
# Execute HEA, using the cif file
#
Phi_HEA = HEA(cif)[0]
# convertes cif to porE xyz (i.e., pypore.xyz) 
ase2pore(cif)
structure = 'pypore.xyz'
#
# Execute porosity evaluation, using the overlapping sphere apporach (OSA)
#
print('-----------')
print('\nRun OSA\n')
print('-----------')
Phi, density, poreV, V_total, V_vdwSum, V_overlap = p.osa(structure)
#
# Execute an analyis of the pore size distribution (PSD)
# First number   -> Number of starting points
# Second number  -> Number of MC steps
#
print('-----------')
print('\nRun PSD\n')
print('-----------')
no_pores,tmp1,tmp2,tmp3,tmp4 = psd.get_psd(structure,200,1000)
pores           = tmp1[0:no_pores]
distr           = tmp2[0:no_pores]
pore_pos_cart   = tmp3[0:no_pores]
pore_pos_frac   = tmp4[0:no_pores]
#
# Execute porosity evaluation, using the grid point apporach (GPA)
# a probe radius for the accessible porosity needs to be provided
#
probe_R = 1.20
#
# Here, explicitly provide the full grid (grid points along each cell vector)
# FullGrid
grid_a = grid_b = grid_c = 30
print('-----------------------------------')
print('\nRun GPA: grid_a, grid_b, grid_c\n')
print('-----------------------------------')
Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_fullgrid(structure,probe_R,grid_a,grid_b,grid_c)
#
# Here, a grid point density per A is provided instead
# GridPerA
grid_density = 2.0
print('-------------------------')
print('\nRun GPA: grid_density\n')
print('-------------------------')
Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_gridpera(structure,probe_R,grid_density)
