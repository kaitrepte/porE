from porE import pore
from ase.io import read 

# usage: python3 pore_run.py

def ase2pore(f_file):
    # convert ase structural information into pore input 
    struct = read(f_file)
    pos = struct.get_positions()
    sym = struct.get_chemical_symbols()
    cell = struct.get_cell()[:]
    f = open('pypore.xyz','w')
    f.write('%i \n' %(len(pos)))
    f.write('%0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f \n' %(cell[0,0],cell[0,1],cell[0,2],cell[1,0],cell[1,1],cell[1,2],cell[2,0],cell[2,1],cell[2,2]))
    for s in range(len(sym)):
        f.write('%s %0.9f %0.9f %0.9f\n' %(sym[s],pos[s][0],pos[s][1],pos[s][2]))
    f.close()

# input: Python 

# from ase to pore 
f_file = 'VOL_03_final.cif'
ase2pore(f_file) 

# ud -- user defined name must be pypore.xyz
system='ud'
method=2 
probe_r =1.2 
grid_a = 100 
grid_b = 100
grid_c = 100

# porE: Fortran call 
pore(system,method,probe_r,grid_a,grid_b,grid_c)
