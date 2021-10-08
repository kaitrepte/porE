from porE.hea.HEA import HEA
from pore import porosity as p
from pore import psd 
from porE.io.ase2pore import *

def run_HEA(cif):
    '''
        Execute HEA, using a cif file

        Input
        -----
        cif ... cif file containing the structural information

        Output
        ------
        Phi_HEA ... Porosity, in %
    '''

    print('-------')
    print('Run HEA')
    print('-------')
    Phi_HEA = HEA(cif)[0]
    return Phi_HEA, True


def run_OSA(xyz):
    '''
        Execute porosity evaluation, 
        using the overlapping sphere approach (OSA)

        Input
        -----
        xyz ... structural information in a xyz format compatible with porE

        Output
        ------
        Phi       ... Porosity, in %
        density   ... density of the structure, in kg/m^3
        poreV     ... pore volume density, in cm^3/g
        V_total   ... Total unit cell volume, in A^3
        V_vdwSum  ... Sum of the volumes of the vdW spheres of all atoms, in A^3
        V_overlap ... Overlap volume, in A^3
    '''

    print('-------')
    print('Run OSA')
    print('-------')
    Phi, density, poreV, V_total, V_vdwSum, V_overlap = p.osa(xyz)
    return Phi, density, poreV, V_total, V_vdwSum, V_overlap, True


def run_GPA_fullgrid(xyz,probe_R,grid_a,grid_b,grid_c):
    '''
        Execute porosity evaluation, using the grid point apporach (GPA)
        Here: Provide an explicit size of the grid in directions a, b, c

        A probe radius for the accessible porosity needs to be provided

        Input
        -----
        xyz     ... structural information in a xyz format compatible with porE
        probe_R ... probe radius, in A
        grid_a  ... number of grid points in crystallographic direction a
        grid_b  ... number of grid points in crystallographic direction b
        grid_c  ... number of grid points in crystallographic direction c

        Output
        ------
        Phi_void   ... Void porosity, in %
        Phi_acc    ... Accessible porosity, in %
        density    ... density of the structure, in kg/m^3
        poreV_void ... pore volume density wrt void porosity, in cm^3/g
        poreV_acc  ... pore volume density wrt accessible porosity, in cm^3/g
    '''

    print('-------------------------------')
    print('Run GPA: grid_a, grid_b, grid_c')
    print('-------------------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_fullgrid(xyz,probe_R,grid_a,grid_b,grid_c)
    return Phi_void, Phi_acc, density, poreV_void, poreV_acc, True


def run_GPA_gridPerA(xyz,probe_R,grid_density):
    '''
        Execute porosity evaluation, using the grid point apporach (GPA)
        Here: Provide a grid density along all cell vectors a, b, c
              porE will determine the amount of grid points automatically

        A probe radius for the accessible porosity needs to be provided

        Input
        -----
        xyz          ... structural information in a xyz format compatible with porE
        probe_R      ... probe radius, in A
        grid_density ... grid point density along all cell vectors, in points/A

        Output
        ------
        Phi_void   ... Void porosity, in %
        Phi_acc    ... Accessible porosity, in %
        density    ... density of the structure, in kg/m^3
        poreV_void ... pore volume density wrt void porosity, in cm^3/g
        poreV_acc  ... pore volume density wrt accessible porosity, in cm^3/g
    '''

    print('---------------------')
    print('Run GPA: grid_density')
    print('---------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_gridpera(xyz,probe_R,grid_density)
    return Phi_void, Phi_acc, density, poreV_void, poreV_acc, True


def run_PSD(xyz,N_points,N_steps):
    '''
        Execute an analyis of the pore size distribution (PSD)

        Input
        -----
        xyz      ... structural information in a xyz format compatible with porE
        N_points ... Number of starting points
        N_steps  ... Number of Monte-Carlo steps

        Output
        ------
        pores         ... all pore diameters, in A
        distr         ... the corresponding ditribution, in %
        pore_pos_cart ... cartesian  coordinates of pore centers, in A
        pore_pos_frac ... fractional coordinates of pore centers
    '''

    print('-------')
    print('Run PSD')
    print('-------')
    no_pores,tmp1,tmp2,tmp3,tmp4 = psd.get_psd(xyz,N_points,N_steps)
    pores           = tmp1[0:no_pores]
    distr           = tmp2[0:no_pores]
    pore_pos_cart   = tmp3[0:no_pores]
    pore_pos_frac   = tmp4[0:no_pores]
    return pores, distr, pore_pos_cart, pore_pos_frac, True


def main(cif):
    '''
        Execute any evaluation you want
    '''
    # start from a cif file used in the HEA;
    # convert cif to porE xyz (i.e., pypore.xyz) 
    # To be used for OSA, PSD and GPA
    ase2pore(cif)
    xyz = 'pypore.xyz'

    ###################
    # parameters to use
    # for PSD: N_points, N_steps
    # Recommendation: N_points > 100, N_steps > 1000
    N_points = 200
    N_steps  = 2000
    # For GPA: probe radius, grid points or grid density (recommended)
    # Recommendation: grid_density > 5
    probe_R = 1.20
    grid_a  = grid_b = grid_c = 30
    grid_density = 5
    # To check which run was performed; for the final summary
    HEA_done = OSA_done = GPA_done = PSD_done = False

    #### Run HEA. Only cif file needed
    #Phi_HEA, HEA_done = run_HEA(cif)
    #### Run OSA. Only xyz file needed
    #Phi_OSA, density, poreV, V_total, V_vdwSum, V_overlap, OSA_done = run_OSA(xyz)
    #### Run PSD. xyz, N_points and N_steps are needed
    #pores, distr, pore_pos_cart, pore_pos_frac, PSD_done = run_PSD(xyz,N_points,N_steps)
    ### Run GPA. xyz, probe radius and grid is needed
    ### Lets provide the total number of grid points per cell vector (not really used anymore)
    #Phi_void, Phi_acc, density, poreV_void, poreV_acc = run_GPA_fullgrid(xyz,probe_R,grid_a,grid_b,grid_c)
    ### Lets provide a grid point density (Preferred method)
    Phi_void, Phi_acc, density, poreV_void, poreV_acc, GPA_done = run_GPA_gridPerA(xyz,probe_R,grid_density)

    # Print a little summary
    print('\n### porE SUMMARY ###')
    if HEA_done:
        print('Phi_HEA       = {:10.5f} %'.format(Phi_HEA))
    if OSA_done:
        print('Phi_OSA       = {:10.5f} %'.format(Phi_OSA))
    if GPA_done:
        print('Phi_GPA(void) = {:10.5f} %'.format(Phi_void))
        print('Phi_GPA(acc)  = {:10.5f} %'.format(Phi_acc))
    if PSD_done:
        form = 'Pores [A]        '+'{:10.5f}'*len(pores)
        print(form.format(*pores))
        form = 'Distribution [%] '+'{:10.2f}'*len(pores)
        print(form.format(*distr))

if __name__ == '__main__':
    # start from a cif file
    #cif = 'structures/cif/NU1000.cif'
    cif = 'nu1000.cif'
    main(cif)
