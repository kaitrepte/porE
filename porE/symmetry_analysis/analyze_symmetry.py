import numpy as np 
from ase.io import read,write 
from ase.atoms import symbols2numbers
from spglib import *
from ase import Atoms 
from ase.visualize import view 

def analyze_symmetry(f_cif,symprec=1e-5):
    # input cell 
    print('Dataset: Input structure')
    ase_atoms = read(f_cif) 
    cf = ase_atoms.get_chemical_formula()
    print('Chemical formula: %s' %(cf))
    sym = ase_atoms.get_chemical_symbols()
    Natoms = len(sym)
    cellparams = ase_atoms.get_cell_lengths_and_angles() 
    print('Cell parameters = ',cellparams)
    print('Cell vectors    = ',ase_atoms.cell[0],ase_atoms.cell[1],ase_atoms.cell[2])
    print('Natoms = %i' % Natoms)
    spacegroup = get_spacegroup(ase_atoms, symprec=symprec)
    print('symprec = %0.5f -> %s' %(symprec,spacegroup))
    print('\n')

    # find primitive cell 
    print('Dataset: Primitive cell')
    lattice, scaled_positions, numbers = find_primitive(ase_atoms, symprec=symprec)
    prim_atoms = Atoms(numbers,scaled_positions=scaled_positions,cell=lattice)
    write(f_cif.split('.')[0]+'_prim.cif',prim_atoms,'cif')
    ase_atoms = prim_atoms
    cf = ase_atoms.get_chemical_formula()
    print(cf)
    sym = ase_atoms.get_chemical_symbols()
    Natoms_prim = len(sym)
    cellparams_prim = ase_atoms.get_cell_lengths_and_angles()
    print('Cell parameters = ',cellparams_prim)
    print('Cell vectors    = ',ase_atoms.cell[0],ase_atoms.cell[1],ase_atoms.cell[2])
    print('Natoms = %i' % Natoms_prim)
    spacegroup_prim = get_spacegroup(ase_atoms, symprec=symprec)
    print('symprec = %0.5f -> %s' %(symprec,spacegroup_prim))
    print('\n')

    # find conventionell cell 
    print('Dataset: Conventional cell')
    ase_atoms = read(f_cif)
    lattice, scaled_positions, numbers = standardize_cell(ase_atoms, to_primitive=False, no_idealize=False, symprec=1e-5)
    conv_atoms = Atoms(numbers,scaled_positions=scaled_positions,cell=lattice)
    write(f_cif.split('.')[0]+'_conv.cif',conv_atoms,'cif')
    ase_atoms = conv_atoms
    cf = ase_atoms.get_chemical_formula()
    print(cf)
    sym = ase_atoms.get_chemical_symbols()
    Natoms_conv = len(sym)
    cellparams_conv = ase_atoms.get_cell_lengths_and_angles()
    print('Cell parameters = ',cellparams_conv)
    print('Cell vectors    = ',ase_atoms.cell[0],ase_atoms.cell[1],ase_atoms.cell[2])
    print('Natoms = %i' % Natoms_conv)
    spacegroup_conv = get_spacegroup(ase_atoms, symprec=symprec)
    print('symprec = %0.5f -> %s' %(symprec,spacegroup_conv))
    print('\n')

    # check what kind of structure is the input structure 
    print('Analysis:')
    if np.isclose(Natoms,Natoms_prim) and np.isclose(cellparams,cellparams_prim).all():
        print('Input structure is probably a primitive cell')
    if np.isclose(Natoms,Natoms_conv) and np.isclose(cellparams,cellparams_conv).all():
        print('Input structure is probably a conventional cell')
    print('\n')

symprec = [1e-5,1e-4,1e-3,1e-2,1e-1]
for sp in symprec:
    analyze_symmetry(f_cif='uio68.cif',symprec=sp)
