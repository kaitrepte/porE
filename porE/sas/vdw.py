import numpy
from ase.io import read 
from sas import get_distance, get_radius 
from sas import timeit 
# Ref.: https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291096-987X%2819990130%2920%3A2%3C217%3A%3AAID-JCC4%3E3.0.CO%3B2-A?saml_referrer

#def get_radius(atom):
#    """ 
#        Get the radii for a given ase.atom 
#    """
#    radii_lcpo = {'H' : 1.20, 'C' : 1.70, 'N' : 1.65, 'O': 1.60, 'P' : 1.90, 'S': 1.90, 'Cl': 1.80}
#    return radii_lcpo[atom.symbol]

def get_Si(ri):
    """
       Get the surface area of a isolated sphere 
       Input:
           ri: radius of species i  
       Output:
           Si: area of isolated sphere 
    """
    return 4*numpy.pi*ri**2

def get_Aij(pi,pj,ri,rj,cell):
    """
       Get two-body overlap area 
       Input:
           pi: coordinate of species i 
           pj: coordinate of species j 
           ri: radius of species i 
           rj: radius of species j 
           cell: ase.cell 
       Output:
           Aij: overlap area 
    """
    dij = get_distance(pi,pj,cell)
    if abs(ri - rj) < dij and dij < 0.62*(ri + rj):
       Aij= 2*numpy.pi*ri*(ri-dij/2-(ri**2-rj**2)/(2*dij))
    else: Aij = 0
    return Aij

@timeit 
def vdw(atoms):
    """ 
        Calculate vdW surface area using 
        Linear Combinations of Pairwise Overlaps (LCPO). 
        Input:
            atoms: ase.atoms 
        Output:
            A_vdw: vdW surface area 
    """
    A_vdw = 0
    for i in range(len(atoms)):
        ri = get_radius(atoms[i])
        A_vdw += get_Si(ri)
        for j in range(len(atoms)):
            if i !=j:
            	rj = get_radius(atoms[j])
            	pi = atoms[i].position 
            	pj = atoms[j].position
            	A_vdw -= get_Aij(pi,pj,ri,rj,atoms.cell)
    print('A_vdw = {} AA**2'.format(A_vdw))
    return A_vdw 

if __name__ == '__main__':
    import glob 
    #f_files = glob.glob('../../examples/structures/cif/porE8/*.cif')
    #for f_file in f_files:
    #    print(f_file)
    #	#f_file ='../../examples/structures/cif/porE8/uio66_vesta.cif'
    #	#f_file = 'chb.sdf'
    #    mof = read(f_file)
    #    vdw(mof)
    f_file = '236.mol'
    mof = read(f_file) 
    vdw(mof)
