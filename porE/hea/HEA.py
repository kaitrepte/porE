import numpy
from ase.io import read
from ase.atoms import Atoms
from ase.visualize import view
# from pybff.uff_bonds import ATOM_data, ATOM_key2idx
import time

ATOM_key2idx = {
    'r': 0,
    'max_coord': 1
}

# N - He covalent radii currently from
#  - https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-26#Sec11
#  - https://www.webelements.com/
ATOM_data = {
    'H': (0.23, 1),
    # 'He': (0.93, 4),
    'He': (1.4, 4),  # vdw
    'Li': (0.68, 1),
    'Be': (0.35, 4),
    'B': (0.83, 4),
    'C': (0.68, 4),
    # 'C': (1.7,4), #vdw
    'N': (0.68, 4),
    'O': (0.68, 2),
    'F': (0.64, 1),
    'Ne': (1.12, 4),
    'Na': (0.97, 1),
    'Mg': (1.10, 4),
    'Al': (1.35, 4),
    'Si': (1.20, 4),
    'P': (1.05, 4),
    'S': (1.02, 4),
    'Cl': (0.99, 1),
    'Ar': (1.57, 4),
    'K': (1.33, 1),
    'Ca': (0.99, 6),
    'Sc': (1.44, 4),
    'Ti': (1.47, 6),
    'V': (1.33, 4),
    'Cr': (1.35, 6),
    'Mn': (1.35, 6),
    'Fe': (1.34, 6),
    'Co': (1.33, 6),
    'Ni': (1.50, 4),
    'Cu': (1.52, 4),
    'Zn': (1.45, 4),
    'Ga': (1.22, 4),
    'Ge': (1.17, 4),
    'As': (1.21, 4),
    'Se': (1.22, 4),
    'Br': (1.21, 1),
    'Kr': (1.91, 4),
    'Rb': (1.47, 1),
    'Sr': (1.12, 6),
    'Y': (1.78, 4),
    'Zr': (1.56, 4),
    'Nb': (1.48, 4),
    'Mo': (1.47, 6),
    'Tc': (1.35, 6),
    'Ru': (1.40, 6),
    'Rh': (1.45, 6),
    'Pd': (1.50, 4),
    'Ag': (1.59, 2),
    'Cd': (1.69, 4),
    'In': (1.63, 4),
    'Sn': (1.46, 4),
    'Sb': (1.46, 4),  # used the same value as for Sn
    'Te': (1.47, 4),
    'I': (1.40, 1),
    'Xe': (1.98, 4),
    'Cs': (1.67, 1),
    'Ba': (1.34, 6),
    'La': (1.87, 4),
    'Ce': (1.83, 6),
    'Pr': (1.82, 6),
    'Nd': (1.81, 6),
    'Pm': (1.80, 6),
    'Sm': (1.80, 6),
    'Eu': (1.99, 6),
    'Gd': (1.79, 6),
    'Tb': (1.76, 6),
    'Dy': (1.75, 6),
    'Ho': (1.74, 6),
    'Er': (1.73, 6),
    'Tm': (1.72, 6),
    'W': (1.62, 6),  # https://www.webelements.com/tungsten/atom_sizes.html
    'Pt': (1.36, 3),  # https://www.webelements.com/platinum/atom_sizes.html
    'Au': (1.36, 2),  # https://www.webelements.com/gold/atom_sizes.html
    'Hg': (1.32, 1)  # https://www.webelements.com/mercury/atom_sizes.html
}


def get_NN(p, rcut):
    """ Get neighbour list
        p    ... positions
        rcut ... cutoff radius
    """
    N = len(p)
    NN = []
    for i in range(N):
        nn = 0
        nlist = []
        for j in range(N):
            n = numpy.linalg.norm(p[i] - p[j])
            if n <= rcut:
                nn += 1
                nlist.append(j)
        NN.append(nlist)
    return NN


def get_subcells(cell, M, p, verbose=3):
    """ Divide a given cell in M subcells
        cell ... 3x3 matrix, cell vectors
        M    ... [Mx,My,Mz] number of subcells
        p    ... positions
    """
    Lx = numpy.linalg.norm(cell[0, :])
    Ly = numpy.linalg.norm(cell[1, :])
    Lz = numpy.linalg.norm(cell[2, :])
    print('L: {} {} {}'.format(Lx, Ly, Lz))
    Mx = M[0]
    My = M[1]
    Mz = M[2]
    N = len(p)
    Celllist = numpy.zeros((Mx, My, Mz, N))
    for i in range(0, N):
        n = numpy.linalg.norm(p[i])
        x, y, z = p[i]
        for mx in range(Mx):
            for my in range(My):
                for mz in range(Mz):
                    cellx_l = Lx / Mx * (mx)
                    celly_l = Ly / My * (my)
                    cellz_l = Lz / Mz * (mz)
                    cellx_h = Lx / Mx * (mx + 1)
                    celly_h = Ly / My * (my + 1)
                    cellz_h = Lz / Mz * (mz + 1)
                    if verbose > 3:
                        print('cell: {}/{} {}/{} {}/{}'.format(cellx_l, cellx_h, celly_l, celly_h, cellz_l, cellz_h))
                        print('m: {} {} {}'.format(mx, my, mz))
                        print('p: {} {} {}'.format(x, y, z))
                    if cellx_l <= x <= cellx_h and celly_l <= y <= celly_h and cellz_l <= z <= cellz_h:
                        if verbose > 3:
                            print('check', x, cellx_h, y, celly_h, z, cellz_h, n)
                        Celllist[mx, my, mz, i] = 1
    return Celllist


def calc_Voverlap(p_i, p_j, s_j, r_i):
    """ Eq. (8) quite expensive
        p_i, p_j : positions
        s_i, s_j: symbols
    """
    d = numpy.linalg.norm(p_j - p_i)
    r_j = ATOM_data[s_j][ATOM_key2idx['r']]
    if abs(r_i - r_j) < d and d < (r_i + r_j):
        res = numpy.pi * (d ** 4 - 6 * d ** 2 * (r_i ** 2 + r_j ** 2) + 8 * d * (r_i ** 3 + r_j ** 3) - 3 * (
                r_i ** 2 - r_j ** 2) ** 2) / (12 * d)
    elif abs(r_i - r_j) > d:
        """ smaller sphere is inside the bigger """
        r_min = numpy.array([r_i, r_j]).min()
        res = 4 / 3. * numpy.pi * r_min ** 3.
    else:
        res = 0
    return res


def get_COMs_celllist(cell, M):
    """ Calculate Center of Mass (COM) for celllist """
    Mx = M[0]
    My = M[1]
    Mz = M[2]
    Lx = numpy.linalg.norm(cell[0, :])
    Ly = numpy.linalg.norm(cell[1, :])
    Lz = numpy.linalg.norm(cell[2, :])
    NEW_P = []
    vec_x = cell[:][0, :].copy()
    n_x = numpy.linalg.norm(vec_x)
    vec_x /= n_x
    vec_y = cell[:][1, :].copy()
    n_y = numpy.linalg.norm(vec_y)
    vec_y /= n_y
    vec_z = cell[:][2, :].copy()
    n_z = numpy.linalg.norm(vec_z)
    vec_z /= n_z
    for mx in range(Mx):
        for my in range(My):
            for mz in range(Mz):
                cellx_l = Lx / Mx * mx
                celly_l = Ly / My * my
                cellz_l = Lz / Mz * mz
                cellx_h = Lx / Mx * (mx + 1)
                celly_h = Ly / My * (my + 1)
                cellz_h = Lz / Mz * (mz + 1)
                # COMs of the subcells for cubic
                # new_p = numpy.array([cellx_l + (cellx_h - cellx_l) / 2., celly_l + (celly_h - celly_l) / 2.,
                #                      cellz_l + (cellz_h - cellz_l) / 2.])
                # COMS of the subcells for lattice basis vectors
                newo = cellx_l * vec_x + celly_l * vec_y + cellz_l * vec_z
                newx = ((cellx_h - cellx_l) / 2.) * vec_x
                newy = ((celly_h - celly_l) / 2.) * vec_y
                newz = ((cellz_h - cellz_l) / 2.) * vec_z
                new_p2 = newo + newx + newy + newz
                NEW_P.append(new_p2)
    return NEW_P


def prune_points(p_framework, p_insert, s_framework, s_insert, r_i):
    """ Prune inserted points
        calculating the overlap of the insert atoms with the framework atoms
        p_framework, p_insert : positions of framework and inserted species
        s_framework, s_insert : symbols of framework and inserted species
        r_i : radii He
    """
    p_insert_red = []
    s_insert_red = []
    for i in range(0, len(p_insert)):
        Voverlap = 0
        for j in range(0, len(p_framework)):
            Voverlap_tmp = calc_Voverlap(p_i=p_insert[i], p_j=p_framework[j], s_j=s_framework[j], r_i=r_i)
            Voverlap += Voverlap_tmp
        if Voverlap < 0.01:
            p_insert_red.append(p_insert[i])
            s_insert_red.append(s_insert[i])
    return p_insert_red, s_insert_red


def get_M(Lx, Ly, Lz, m='low', d=2 * ATOM_data['He'][ATOM_key2idx['r']]):
    """ Determine optimal subcell sizes
        Lx, Ly, Lz ... cell lengths
        m          ... method 'low' or 'high'
        d          ... insert Atom optimal distance for insertion lattice
    """
    f = {'low': numpy.floor,
         'high': numpy.ceil}
    Mx = int(f[m](numpy.round(Lx, 6) / d))
    My = int(f[m](numpy.round(Ly, 6) / d))
    Mz = int(f[m](numpy.round(Lz, 6) / d))
    calc_error(numpy.array([Lx, Ly, Lz]), numpy.array([Mx, My, Mz]), d)
    return Mx, My, Mz


def get_opt_r(Lx, Ly, Lz, Mx, My, Mz):
    """ Determine d_opt
        Lx, Ly, Lz : lengths of the cell vectors
        Mx, My, Mz : number of subcells in x,y,z
    """
    dx = Lx / Mx
    dy = Ly / My
    dz = Lz / Mz
    d = numpy.array([dx, dy, dz])
    dopt = d.sum() / 3.
    return dopt


def calc_error(L, M, d_target):
    """ Calculate spacing error
        L = (Lx, Ly, Lz) : lengths of the cell vectors
        M = (Mx, My, Mz) : number of subcells in x,y,z
        d_target : target spacing
    """
    d_calc = L / M
    e = abs(d_calc - d_target)
    print('spacing error(|d_calc-d_target|): {}'.format(e))
    return e


def HEA(f_file, verbose=3):
    """ Helium approach
          Calculate porosity using celllist
          f_file  ... file name
          verbose ... verbosity
      """
    #print('----------------')
    #print('HEA: calculation')
    #print('----------------')
    print('f_file : {}'.format(f_file))
    t1 = time.time()
    # read the structure
    struct = read(f_file)
    p_framework = struct.get_positions()
    s_framework = struct.get_chemical_symbols()
    cell = struct.get_cell()
    Lx = numpy.linalg.norm(cell[:][0, :])
    Ly = numpy.linalg.norm(cell[:][1, :])
    Lz = numpy.linalg.norm(cell[:][2, :])
    # Guess for M
    M0 = get_M(Lx=Lx, Ly=Ly, Lz=Lz, m='high')
    # Optimize d for M0
    dopt = get_opt_r(Lx=Lx, Ly=Ly, Lz=Lz, Mx=M0[0], My=M0[1], Mz=M0[2])
    # Optimize M for dopt
    M = get_M(Lx=Lx, Ly=Ly, Lz=Lz, m='high', d=dopt)
    Mx, My, Mz = M
    # Celllist = get_subcells(cell, M, p_framework)
    # Get COMs of Celllist
    p_insert = get_COMs_celllist(cell, M)
    # Species we are inserting
    s_insert = ['He'] * len(p_insert)
    # Remove atoms an the insertion grid overlapping with the framework
    ropt = dopt / 2.
    print('r_opt:            {:10.5f} A'.format(ropt))
    p_insert_red, s_insert_red = prune_points(p_framework=p_framework, p_insert=p_insert, s_framework=s_framework,
                                              s_insert=s_insert, r_i=ropt)
    # Calculate number of atoms
    N_insert = len(p_insert)
    N_insert_red = len(p_insert_red)
    P_points = N_insert_red / N_insert
    print('N_insert:         {:10d}'.format(N_insert))
    print('N_insert_reduced: {:10d}'.format(N_insert_red))
    # Calculate volumes
    struct_subcell = Atoms(cell=cell)
    cell_subcell = struct_subcell.get_cell()
    cell_subcell[:][0, :] = cell_subcell[:][0, :] / Mx
    cell_subcell[:][1, :] = cell_subcell[:][1, :] / My
    cell_subcell[:][2, :] = cell_subcell[:][2, :] / Mz
    struct_subcell.set_cell(cell_subcell)
    V_subcell = struct_subcell.get_volume() * N_insert_red
    V_insert = 4 / 3. * numpy.pi * ropt ** 3. * N_insert_red
    V_framework = struct.get_volume()
    P_subcell = V_subcell / V_framework
    print('optimal M:        {:10d} {:10d} {:10d}'.format(Mx, My, Mz))
    # simple cubic packing
    f_sc = 0.52
    f_calc = V_insert / V_subcell
    print('------------------')
    print('packing fraction f')
    print('------------------')
    print('f_sc:             {:15.5f}'.format(f_sc))
    print('f_calc:           {:15.5f}'.format(f_calc))
    P_volume = V_insert / (V_framework * f_calc)
    print('----------------------')
    print('HEA porosity P [%]')
    print('----------------------')
    print('P_subcell:        {:15.5f}'.format(P_subcell*100.0))
    print('P_points:         {:15.5f}'.format(P_points*100.0))
    print('P_volume:         {:15.5f}'.format(P_volume*100.0))
    print('---------------')
    print(' Volume V [A^3]')
    print('---------------')
    print('V_framework:      {:15.5f}'.format(V_framework))
    print('V_subcell:        {:15.5f}'.format(V_subcell))
    print('V_He:             {:15.5f}'.format(V_insert))
    t2 = time.time()
    print('Total CPU time:   {:15.5f} s'.format(t2-t1))
    print('\n')

    if verbose > 3:
        get_ase(f_file, p_insert_red, s_insert_red)
    return [P_points*100.0, P_subcell*100.0, P_volume*100.0, V_framework, V_subcell, V_insert, N_insert, N_insert_red, f_calc, ropt]


def get_ase(f_file, p_insert_red, s_insert_red):
    """ Combine framework and insert atoms to one ase structure
        f_file : framework file name
        p_insert_red: positions of reduced inserted species
        s_insert_red: symbols of reduced inserted species
    """
    struct = read(f_file)
    insert = Atoms(s_insert_red, p_insert_red)
    struct_tot = struct + insert
    view(struct_tot)


if __name__ == '__main__':
    # Before small benchmark
    t1 = time.time()
    # f_file = ['hkust1_vesta.cif','uio66_vesta.cif','irmof10_vesta.cif','uio67_vesta.cif','dut8open_vcrelaxed_vesta_conv.cif','uio66_vesta_conv.cif','uio67_vesta_conv.cif','irmof10_vesta_conv.cif','mof5_vesta_conv.cif','hkust1_vesta_conv.cif']
    f_file = ['dut8closed_exp_vesta.cif', 'dut8open_exp_vesta.cif', 'uio66_vesta.cif', 'uio67_vesta.cif',
              'irmof10_vesta.cif', 'mof5_vesta.cif', 'hkust1_vesta.cif', 'mof210_vesta.cif']
    res = []
    for f in f_file:
        P = HEA(f_file=f, verbose=4)
        res.append(P)
    print('System P_N P_subcell P_volume')
    print('-----------------------------')
    for i in range(len(f_file)):
        print('{} : {}'.format(f_file[i], res[i]))
    t2 = time.time()
    print('-----------')
    print('Time : {} s'.format(t2 - t1))
