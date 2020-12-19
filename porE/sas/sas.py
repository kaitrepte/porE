import math 
import numpy 
from matplotlib.pyplot import * 
from ase.atoms import Atoms 
from ase.data import atomic_numbers, covalent_radii, vdw_radii
from ase.io import read, write 
from ase.geometry import find_mic
import time 

# Origional File 
# - https://github.com/boscoh/pdbremix/blob/master/pdbremix/asa.py
# - http://boscoh.com/protein/calculating-the-solvent-accessible-surface-area-asa.html
# SS modifications 
# - python2 to python3 
# - remove v3 dependencies 
# - add pbc 
# - add total SAS area 
# - add timing 

def timeit(f):
    """
        Decorator for timing 
    """
    def f0(*args, **kwargs):
        before = time.time()
        res = f(*args, **kwargs)
        after = time.time()
        print('elapsed time ({}) = {} s'.format(f.__qualname__, after - before))
        return res
    return f0

def get_distance(p1,p2,cell):
    """ 
        Get pbc distance between point p1 and p2 in cell
    """
    v = p1 - p2 
    v_min = find_mic(v,cell)[1]
    vlen = numpy.linalg.norm(p1-p2)
    return v_min 

def get_radius(atom): 
    """ 
        Get the radii for a given ase.atom 
    """
    #tmp = vdw_radii.tolist() 
    #tmp.append([40,2.36])
    #vdw_radii = numpy.array(tmp)
    tmp = vdw_radii.copy()
    tmp[40] = 2.36 
    #print(tmp[40])
    return tmp[atomic_numbers[atom.symbol]]
    #return covalent_radii[atomic_numbers[atom.symbol]]

def get_points_on_sphere(n):
    """
        Get points on a sphere 
        using the Golden-Section Spiral algorithm.
        Input: 
            n: number of points 
        Output: 
            points: coordinates for points on a sphere 
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return points


def find_neighbor_indices(atoms, r_probe, k):
    """
        Returns list of indices of atoms within probe distance to atom k. 
        Input: 
            atoms: ase.atoms 
            r_probe: probe radius in [AA]
            k: index
        Output: 
            neighbor_indices: indices 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = get_radius(atom_k) + r_probe + r_probe
    indices = [i for i in range(k)]
    indices.extend([i for i in range(k+1, len(atoms))])
    for i in indices:
      atom_i = atoms[i]
      dist = get_distance(atom_k.position, atom_i.position,atoms.cell)
      if dist < radius + get_radius(atom_i):
        neighbor_indices.append(i)
    return neighbor_indices

@timeit
def sas(atoms, r_probe=1.4, n_points=960, verbose=3, name='SAS'):
    """
        Returns the Solvent Accessible Surface (SAS), by rolling a
        ball with probe radius over the atoms with their radius
        defined.

        Input:
            atoms: ase.atoms
            r_probe: probe radius [AA]
            n_points: number of points on a sphere
            verbose: verbosity
        Output:
            res : vdW surface area [AA**2]

    """
    sphere_points = get_points_on_sphere(n_points)
    
    const = 4.0 * math.pi / len(sphere_points)
    areas = []
    for i, atom_i in enumerate(atoms):
      
      neighbor_indices = find_neighbor_indices(atoms, r_probe, i)
      n_neighbor = len(neighbor_indices)
      j_closest_neighbor = 0
      radius = r_probe + get_radius(atom_i)
      
      # radius defininition by PoreBlazer
      # radius *= 2**(1/6.)
    
      n_accessible_point = 0
      for point in sphere_points:
        is_accessible = True
        test_point = numpy.array(point)*radius + atom_i.position 
        cycled_indices = [i for i in range(j_closest_neighbor, n_neighbor)]
        cycled_indices.extend([i for i in range(j_closest_neighbor)])
        
        for j in cycled_indices:
          atom_j = atoms[neighbor_indices[j]]
          r = get_radius(atom_j) + r_probe
          diff = get_distance(atom_j.position, test_point,atoms.cell)
          if diff*diff < r*r:
            j_closest_neighbor = j
            is_accessible = False
            break
        if is_accessible:
          n_accessible_point += 1
      
      area = const*n_accessible_point*radius*radius 
      if verbose > 3: 
          print('area = {} n = {} radius = {}'.format(area,n_accessible_point,radius))
      areas.append(area)
    
    A = numpy.array(areas)
    # Sum all atomic ASA contributions 
    SAS = A.sum()
    if verbose == 3: 
        print('{}(r_probe = {} AA, n_points = {}) = {} AA**2'.format(name,r_probe,n_points,SAS))
    return SAS

def write_points_on_sphere_xyz(n_points=96,symbol='X'):
    """ 
        Write points on sphere as xyz file. 
        Input: 
            n_points : number of points on a sphere 
            symbol: chemical symbol/ species 
        Output: 
            xyz file: xyz file with n_sphere_point coordinates as symbol species 
    """
    Nsphere = get_points_on_sphere(n=n_points)
    ase_atoms = Atoms(len(Nsphere)*symbol,Nsphere)
    write('n_{}_points_on_sphere.xyz'.format(n_points),ase_atoms,'xyz') 

def vdw(atoms,n_points,verbose=3):
    """ 
        Calculate vdW surface using SAS algorithm with r_probe = 0 AA.
        Input: 
            atoms: ase.atoms 
            n_points: number of points on a sphere 
            verbose: verbosity 
        Output: 
            res : vdW surface area [AA**2]
    """
    res = sas(atoms, r_probe=0, n_points=n_points, verbose=verbose, name='vdW')
    return res 

if __name__ == '__main__':

    f_file ='../../examples/structures/cif/porE8/uio66_vesta.cif'
    r_probe = 1.2
    n_points = 9 
    mof = read(f_file)
    write_points_on_sphere_xyz(n_points=n_points)
    A_sas = sas(mof, r_probe=r_probe, n_points=n_points)
    A_vdw = vdw(mof, n_points=n_points)


