from ase.io import read

def ase2pore(f_file):
    # convert ase structural information into pore input 
    struct = read(f_file)
    pos = struct.get_positions()
    sym = struct.get_chemical_symbols()
    cell = struct.get_cell()[:]
    name = f_file.split('/')[-1].replace('.cif','')
    f = open(name+'.xyz','w')
    f.write('%i \n' %(len(pos)))
    f.write('%0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f \n' %(cell[0,0],cell[0,1],cell[0,2],cell[1,0],cell[1,1],cell[1,2],cell[2,0],cell[2,1],cell[2,2]))
    for s in range(len(sym)):
        f.write('%s %0.9f %0.9f %0.9f\n' %(sym[s],pos[s][0],pos[s][1],pos[s][2]))
    f.close()


