import pore
# two modules: porosity and psd

# define some abbreviations
osa          = pore.porosity.osa
gpa_FullGrid = pore.porosity.gpa_fullgrid
gpa_GridPerA = pore.porosity.gpa_gridpera
get_PSD      = pore.psd.get_psd
#
# Structure of the xyz file:
#   number of atoms
#   cell vectors as   a_x a_y a_z  b_x b_y b_z  c_x c_y c_z
#   Element specifiers and coordinates
##
# see folder 'structures_xyz' as well
#
structure = 'structures_xyz/uio66.xyz'
#
# Execute porosity evaluation, using the overlapping sphere apporach (OSA)
#
print(' ')
osa(structure)
#
# Execute an analyis of the pore size distribution (PSD)
# First number   -> Number of starting points
# Second number  -> Number of MC steps
get_PSD(structure,100,1000)
#####
# Execute porosity evaluation, using the grid point apporach (GPA)
# a probe radius for the accessible porosity needs to be provided
probe_R = 1.20
#
# Here, explicitely provide the full grid (grid points along each cell vector)
# FullGrid
grid_a  = 30
grid_b  = 30
grid_c  = 30
print(' ')
gpa_FullGrid(structure,probe_R,grid_a,grid_b,grid_c)
#
# Here, a grid point density per A is provided instead
# GridPerA
grid_density = 2.0
print(' ')
gpa_GridPerA(structure,probe_R,grid_density)

